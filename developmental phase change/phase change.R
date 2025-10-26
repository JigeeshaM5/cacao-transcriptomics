install.packages("igraph")  # Install if missing
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
library(DESeq2)
library(ggplot2)
library(pheatmap)
BiocManager::install(c("DESeq2", "pheatmap"))
install.packages("ggplot2")  
install.packages("visNetwork")
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(clusterProfiler)
library(ggplot2)
library(igraph)  



#path to your file
# Define file path
raw_phase <- read.csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/for R/deseq2_OTA_vs_PTA.csv")

#read other files

expr_data_early<-read.csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/deseq2_OTA_early_vs_PTA_early.csv")
expr_data_late<-read.csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/deseq2_OTA_late_vs_PTA_late.csv")

library(ggplot2)
library(matrixStats)
library(pheatmap)



#no filtering from here July 23
library(writexl)
library(readxl)
# Step 2: Convert necessary columns to numeric
raw_phase$log2FoldChange <- as.numeric(as.character(raw_phase$log2FoldChange))
raw_phase$PAdj           <- as.numeric(as.character(raw_phase$PAdj))
raw_phase$FDR            <- as.numeric(as.character(raw_phase$FDR))
raw_phase$baseMean       <- as.numeric(as.character(raw_phase$baseMean))

#optional filtering july 25---------------------
#first inspect
# Inspect distribution of log2FoldChange and PAdj
summary(raw_phase$log2FoldChange)
summary(raw_phase$PAdj)

# See how many pass each criterion
nrow(subset(raw_phase, abs(log2FoldChange) > 1))         # LFC filter only
nrow(subset(raw_phase, PAdj <= 0.1))                     # Padj filter only
nrow(subset(raw_phase, abs(log2FoldChange) > 1 & PAdj <= 0.1))  # Combined

# View entries with significant raw p-value (< 0.01)
subset(raw_phase, PValue < 0.01)

# Check how many pass log2FC > 1 & raw P-value < 0.01
nrow(subset(raw_phase, abs(log2FoldChange) > 1 & PValue < 0.01))

# Filter using raw p-value and LFC
filtered_by_pvalue <- subset(raw_phase,
                             abs(log2FoldChange) > 1 &
                               PValue < 0.01)

# Save to Excel
library(writexl)
write_xlsx(filtered_by_pvalue,
           "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/for R/filtered_rawP.xlsx")

# Step 4: Save the filtered data as an Excel file
# First, load the writexl package (you may need to install it)
install.packages("writexl")  # Run only once if not installed
library(writexl)
library(ggplot2)

ggplot(raw_phase, aes(x = log2FoldChange, y = -log10(PValue))) +
  geom_point(alpha = 0.5) +
  geom_point(data = filtered_by_pvalue, color = "red") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)")


#or better volcano
library(ggplot2)
library(ggrepel)

# Add a column for significance
raw_phase$Significance <- "Not Significant"
raw_phase$Significance[raw_phase$log2FoldChange > 1 & raw_phase$PValue < 0.01] <- "Upregulated"
raw_phase$Significance[raw_phase$log2FoldChange < -1 & raw_phase$PValue < 0.01] <- "Downregulated"

# Add -log10(p-value)
raw_phase$negLog10P <- -log10(raw_phase$PValue)

# Volcano plot
ggplot(raw_phase, aes(x = log2FoldChange, y = negLog10P, color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  geom_text_repel(data = subset(raw_phase, PValue < 0.01 & abs(log2FoldChange) > 1),
                  aes(label = rownames(subset(raw_phase, PValue < 0.01 & abs(log2FoldChange) > 1))),
                  size = 3, max.overlaps = 50) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(p-value)",
       color = "Gene Regulation") +
  theme(legend.position = "top")



#ends here july 25----------------------------------

# Step 3: Filter out genes with log2FC == 0
filtered_phase <- raw_phase[raw_phase$log2FoldChange != 0, ]

# Step 4: Compute composite score
filtered_phase$score <- -log10(filtered_phase$PAdj + 1e-10) * abs(filtered_phase$log2FoldChange)

# Step 5: Rank and select top 298 genes
ranked_phase <- filtered_phase[!is.na(filtered_phase$score), ]
top298_genes <- ranked_phase[order(ranked_phase$score, decreasing = TRUE)[1:298], ]

# Step 6: Sort by |log2FC|, then PAdj, FDR, baseMean
top298_sorted <- top298_genes[
  order(-abs(top298_genes$log2FoldChange),
        top298_genes$PAdj,
        top298_genes$FDR,
        -top298_genes$baseMean),
]

# Step 7: Split into upregulated and downregulated
upregulated   <- top298_sorted[top298_sorted$log2FoldChange > 0, ]
downregulated <- top298_sorted[top298_sorted$log2FoldChange < 0, ]

# Step 8: Save to Excel with two sheets
output_path <- "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_DEGs_up_down.xlsx"

write_xlsx(
  list(
    Upregulated = upregulated,
    Downregulated = downregulated
  ),
  path = output_path
)

#now repeat for early
# Step 2: Convert necessary columns to numeric
expr_data_early$log2FoldChange <- as.numeric(as.character(expr_data_early$log2FoldChange))
expr_data_early$PAdj           <- as.numeric(as.character(expr_data_early$PAdj))
expr_data_early$FDR            <- as.numeric(as.character(expr_data_early$FDR))
expr_data_early$baseMean       <- as.numeric(as.character(expr_data_early$baseMean))

# Step 3: Filter out genes with log2FC == 0
filtered_phase_early <- expr_data_early[expr_data_early$log2FoldChange != 0, ]

# Step 4: Compute composite score
filtered_phase_early$score <- -log10(filtered_phase_early$PAdj + 1e-10) * abs(filtered_phase_early$log2FoldChange)

# Step 5: Rank and select top 224 genes
ranked_phase_early <- filtered_phase_early[!is.na(filtered_phase_early$score), ]
top224_genes_early <- ranked_phase_early[order(ranked_phase_early$score, decreasing = TRUE)[1:224], ]

# Step 6: Sort by |log2FC|, then PAdj, FDR, baseMean
top224_sorted_early <- top224_genes_early[
  order(-abs(top224_genes_early$log2FoldChange),
        top224_genes_early$PAdj,
        top224_genes_early$FDR,
        -top224_genes_early$baseMean),
]

# Step 7: Split into upregulated and downregulated
upregulated   <- top224_sorted_early[top224_sorted_early$log2FoldChange > 0, ]
downregulated <- top224_sorted_early[top224_sorted_early$log2FoldChange < 0, ]

# Step 8: Save to Excel with two sheets
output_path <- "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_early_DEGs_up_down.xlsx"

write_xlsx(
  list(
    Upregulated = upregulated,
    Downregulated = downregulated
  ),
  path = output_path
)


#now repeat for late
# Step 2: Convert necessary columns to numeric
expr_data_late$log2FoldChange <- as.numeric(as.character(expr_data_late$log2FoldChange))
expr_data_late$PAdj           <- as.numeric(as.character(expr_data_late$PAdj))
expr_data_late$FDR            <- as.numeric(as.character(expr_data_late$FDR))
expr_data_late$baseMean       <- as.numeric(as.character(expr_data_late$baseMean))

# Step 3: Filter out genes with log2FC == 0
filtered_phase_late <- expr_data_late[expr_data_late$log2FoldChange != 0, ]

# Step 4: Compute composite score
filtered_phase_late$score <- -log10(filtered_phase_late$PAdj + 1e-10) * abs(filtered_phase_late$log2FoldChange)

# Step 5: Rank and select top 186 genes
ranked_phase_late <- filtered_phase_late[!is.na(filtered_phase_late$score), ]
top186_genes_late <- ranked_phase_late[order(ranked_phase_late$score, decreasing = TRUE)[1:186], ]

# Step 6: Sort by |log2FC|, then PAdj, FDR, baseMean
top186_sorted_late <- top186_genes_late[
  order(-abs(top186_genes_late$log2FoldChange),
        top186_genes_late$PAdj,
        top186_genes_late$FDR,
        -top186_genes_late$baseMean),
]

# Step 7: Split into upregulated and downregulated
upregulated   <- top186_sorted_late[top186_sorted_late$log2FoldChange > 0, ]
downregulated <- top186_sorted_late[top186_sorted_late$log2FoldChange < 0, ]

# Step 8: Save to Excel with two sheets
output_path <- "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_late_DEGs_up_down.xlsx"

write_xlsx(
  list(
    Upregulated = upregulated,
    Downregulated = downregulated
  ),
  path = output_path
)

#-------------------------ignore filtering below

# Ensure relevant columns are numeric
raw_phase$log2FoldChange <- as.numeric(as.character(raw_phase$log2FoldChange))
raw_phase$FDR            <- as.numeric(as.character(raw_phase$FDR))
# Make sure PAdj is numeric
raw_phase$log2FoldChange <- as.numeric(as.character(raw_phase$log2FoldChange))
raw_phase$PAdj <- as.numeric(as.character(raw_phase$PAdj))


# OR use raw_phase$PAdj if you prefer adjusted p-values

# Remove NA values and filter
filtered_genes <- raw_phase[
  !is.na(raw_phase$log2FoldChange) &
    !is.na(raw_phase$FDR) &
    abs(raw_phase$log2FoldChange) >= lfc_cut &
    raw_phase$FDR <= fdr_cut,
]

# Filter using log2FC and PAdj
filtered_PAdj <- raw_phase[
  !is.na(raw_phase$log2FoldChange) &
    !is.na(raw_phase$PAdj) &
    abs(raw_phase$log2FoldChange) >= lfc_cut &
    raw_phase$PAdj <= padj_cut,
]


#continue on july 23 to create heatmaps and volcano plots

library(tidyverse)

# Select only expression columns (adjust as needed)
expr_cols <- grep("OTA_|PTA_", colnames(top298_sorted), value = TRUE)

# Extract and scale expression values
expression_matrix <- as.matrix(top298_sorted[, expr_cols])
row.names(expression_matrix) <- top298_sorted$gene

# Z-score normalization per gene (row-wise)
expression_scaled <- t(scale(t(expression_matrix)))

library(pheatmap)

# Sample group annotation (adjust if needed)
sample_groups <- data.frame(
  Condition = ifelse(grepl("OTA", expr_cols), "OTA", "PTA")
)
row.names(sample_groups) <- expr_cols

# Beautiful heatmap
pheatmap(
  expression_scaled,
  annotation_col = sample_groups,
  cluster_rows = TRUE,     # Cluster genes
  cluster_cols = TRUE,     # Cluster samples
  show_rownames = FALSE,   # Prevent overcrowding
  fontsize = 12,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Expression Heatmap: Top DEGs"
)

# 0. Dependencies
install.packages(c("readr","dplyr","ggplot2","ggrepel","pheatmap","biomaRt","writexl"))
library(readr); library(dplyr)
library(ggplot2); library(ggrepel)
library(pheatmap)
library(biomaRt)
library(writexl)

# 1. Load your 298-gene table
top298 <- read_xlsx(
  "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_DEGs_up_down.xlsx")

# 2. Clean & rank
top298 <- top298 %>%
  mutate(
    log2FoldChange = as.numeric(log2FoldChange),
    PAdj           = as.numeric(PAdj),
    FDR            = as.numeric(FDR),
    baseMean       = as.numeric(baseMean),
    score          = -log10(PAdj + 1e-10) * abs(log2FoldChange)
  ) %>%
  arrange(desc(score))

# 3. Extract top 20 by score
top20 <- top298 %>% slice(1:20)

# Install & load biomaRt (only run install once)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)

# Query the Plants host directly
plants_marts <- listMarts(host = "https://plants.ensembl.org")
print(plants_marts)


library(biomaRt)

#  List available marts on the Plants site
listMarts(host = "https://plants.ensembl.org")
# you should see “plants_mart” listed

# Connect using the base host only
cacao_mart <- useMart(
  biomart = "plants_mart",
  dataset = "tcacao_eg_gene",
  host    = "https://plants.ensembl.org"
)

# Double-check your dataset is available
ds <- listDatasets(cacao_mart)
if (!"tcacao_eg_gene" %in% ds$dataset) stop("Dataset not found—check mart/host.")

# Pull attributes
genes <- top298_sorted$gene
annotations <- getBM(
  mart       = cacao_mart,
  attributes = c("ensembl_gene_id","external_gene_name","description"),
  filters    = "ensembl_gene_id",
  values     = genes
)

# Merge back
top298_annotated <- merge(
  top298_sorted, annotations,
  by.x = "gene", by.y = "ensembl_gene_id",
  all.x = TRUE
)

# Connect to the Ensembl Plants mart for Theobroma cacao
library(biomaRt)


#— Prepare Expression Matrix with Gene Labels —#
expr_cols <- grep("OTA_|PTA_", colnames(top298_annotated), value = TRUE)

expr_mat <- as.matrix(top298_annotated[ , expr_cols])
# Use external_gene_name when available, otherwise fall back to Ensembl ID
row_labels <- ifelse(
  is.na(top298_annotated$external_gene_name) | top298_annotated$external_gene_name == "",
  top298_annotated$gene,
  top298_annotated$external_gene_name
)
rownames(expr_mat) <- row_labels

# Z-score normalization per gene (row)
expr_scaled <- t(scale(t(expr_mat)))


#— Step 7: Draw Heatmap —#
library(pheatmap)

# Define sample annotations
sample_groups <- data.frame(
  Condition = factor(ifelse(grepl("OTA", expr_cols), "OTA", "PTA"))
)
rownames(sample_groups) <- expr_cols
library(pheatmap)

#  Check matching names
all(colnames(expr_scaled) == rownames(sample_groups))
# Ideally returns TRUE

library(pheatmap)

# Ensure it’s a factor (not character)
sample_groups$Condition <- factor(
  sample_groups$Condition,
  levels = c("OTA","PTA")
)

# Define colors for each level
ann_colors <- list(
  Condition = c(
    OTA = "#1b9e77",
    PTA = "#d95f02"
  )
)

# (Re)order to be safe
sample_groups <- sample_groups[colnames(expr_scaled), , drop = FALSE]

# Draw the heatmap with annotation
pheatmap(
  expr_scaled,
  annotation_col    = sample_groups,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  fontsize_col      = 10,
  color             = colorRampPalette(c("navy","white","firebrick3"))(100),
  main              = "Top DEGs (Z-score normalized)"
)

library(pheatmap)

library(ComplexHeatmap)
library(circlize)

# Create column annotation
col_anno <- HeatmapAnnotation(
  Condition = sample_groups$Condition,
  col       = list(Condition = ann_colors$Condition),
  annotation_name_side = "left"
)

Heatmap(
  expr_scaled,
  name               = "Z-score",
  top_annotation     = col_anno,
  show_row_names     = FALSE,
  show_column_names  = TRUE,
  cluster_rows       = TRUE,
  cluster_columns    = TRUE,
  col                = colorRamp2(c(-2,0,4), c("navy","white","firebrick3")),
  column_title       = "Samples",
  heatmap_legend_param = list(title = "Z-score")
)


#volcano plots
library(ggplot2)
library(dplyr)

#stop july 23 at 3.45 pm

# Assuming your dataframe is already ranked and cleaned as `top298`

# Volcano plot

top298 <- top298 %>%
  filter(!is.na(PAdj), is.finite(PAdj), !is.na(log2FoldChange), is.finite(log2FoldChange))

ggplot(top298, aes(x = log2FoldChange, y = -log10(PAdj))) +
  geom_point(alpha = 0.6, size = 2.2, color = "#3b528b") +
  scale_x_continuous(limits = c(-4, 4)) +  # expand axis to make room
  scale_y_continuous(limits = c(0, 10)) +  # stretch vertical to see dispersion
  labs(
    title = "Volcano Plot: OTA vs PTA",
    x = "log₂(Fold Change)",
    y = "−log₁₀(Adjusted p-value)"
  ) +
  theme_minimal(base_size = 15)

summary(top298$log2FoldChange)
summary(-log10(top298$PAdj))

#--------------ignore filtering above



# Install and load the writexl package
install.packages("writexl")
library(writexl)

# Define the output file path
phase_file_filtered <- "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/filtered_DEGs_OTA_vs_PTA.xlsx"

# Save the filtered data frame to Excel
write_xlsx(filtered_genes, phase_file_filtered)

hist(as.numeric(raw_phase$log2FoldChange), breaks = 100,
     main = "Distribution of Log2 Fold Changes",
     xlab = "log2(FC)", col = "steelblue")
abline(v = c(-1.5, 1.5), col = "red", lwd = 2, lty = 2)

hist(as.numeric(raw_phase$FDR), breaks = 100,
     main = "Distribution of FDR values",
     xlab = "FDR", col = "forestgreen")
abline(v = 0.1, col = "red", lwd = 2, lty = 2)

plot(as.numeric(raw_phase$baseMean),
     as.numeric(raw_phase$log2FoldChange),
     log = "x", pch = 20, col = "gray",
     xlab = "BaseMean (log scale)", ylab = "log2FoldChange",
     main = "Expression vs Fold Change")
abline(h = c(-1.5, 1.5), col = "red", lwd = 2, lty = 2)

plot(log2(as.numeric(raw_phase$baseMean + 1)),
     as.numeric(raw_phase$log2FoldChange),
     pch = 20, col = "gray",
     xlab = "log2(baseMean + 1)", ylab = "log2FoldChange",
     main = "MA-style Plot")

# 1. Identify sample columns
sampleCols <- grep("^(OTA|PTA)_.*_R[1-4]$", colnames(raw_phase), value = TRUE)

# 2. Force-convert sample columns to numeric
raw_phase[sampleCols] <- lapply(raw_phase[sampleCols], function(x) as.numeric(as.character(x)))

# 3. Build expression matrix
expr_matrix <- as.matrix(raw_phase[, sampleCols])
stopifnot(is.numeric(expr_matrix))

# 4. Sample metadata
colData_phase <- data.frame(
  condition = factor(ifelse(grepl("^OTA", sampleCols), "OTA", "PTA")),
  replicate = factor(sub("^.*_R", "", sampleCols)),
  row.names = sampleCols
)

# 5. Filter top variable genes
geneVars <- rowVars(expr_matrix)
topN <- 500
topGenes <- order(geneVars, decreasing = TRUE)[1:topN]
filtered_matrix <- expr_matrix[topGenes, ]

# 6. PCA
pca_res <- prcomp(t(filtered_matrix), center = TRUE, scale. = TRUE)
varExp <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# 7. Scree plot
barplot(
  varExp[1:10] * 100,
  names.arg = paste0("PC", 1:10),
  ylab = "% Variance",
  main = paste("Scree Plot (Top", topN, "Genes)")
)

# 8. PC1 vs PC2 plot
pc_df <- data.frame(
  sample    = rownames(pca_res$x),
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  condition = colData_phase$condition
)

ggplot(pc_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 5) +
  geom_text(vjust = -1, size = 3) +
  xlab(paste0("PC1 (", round(varExp[1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(varExp[2] * 100, 1), "%)")) +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle(paste("PCA: Top", topN, "Variable Genes"))


library(ggplot2)

# Assuming 'pc_df' is your PCA result data frame
ggplot(pc_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -0.8, hjust = 0.5, size = 3, check_overlap = TRUE) +
  labs(
    title = "Principal Component Analysis of Top 500 Variable Genes",
    subtitle = "Cacao Phase Change: OTA vs PTA",
    x = "PC1 (50.8% of variance)",
    y = "PC2 (12.6% of variance)",
    color = "Condition"
  ) +
  scale_color_manual(values = c("OTA" = "#E41A1C", "PTA" = "#377EB8")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90")
  )

# 9. Extract PCA loadings (rotation matrix)
loadings <- pca_res$rotation  # dimensions: genes × PCs

# 10. Top 20 genes contributing to PC1
top_PC1_genes <- sort(abs(loadings[, 1]), decreasing = TRUE)[1:20]
top_PC1_names <- names(top_PC1_genes)

# 11. Top 20 genes contributing to PC2
top_PC2_genes <- sort(abs(loadings[, 2]), decreasing = TRUE)[1:20]
top_PC2_names <- names(top_PC2_genes)
# 12. Map gene IDs to gene symbols (if available)
PC1_gene_symbols <- raw_phase[top_PC1_names, "gene"]
PC2_gene_symbols <- raw_phase[top_PC2_names, "gene"]

# Combine into a data frame
top_PC1_df <- data.frame(
  GeneID   = top_PC1_names,
  Symbol   = PC1_gene_symbols,
  Loading  = loadings[top_PC1_names, 1]
)

top_PC2_df <- data.frame(
  GeneID   = top_PC2_names,
  Symbol   = PC2_gene_symbols,
  Loading  = loadings[top_PC2_names, 2]
)

# View results
print("Top PC1 Loading Genes:")
print(top_PC1_df)

print("Top PC2 Loading Genes:")
print(top_PC2_df)
# Combine top gene IDs from both PCs
top_genes_combined <- unique(c(top_PC1_names, top_PC2_names))

# Subset expression matrix
heatmap_matrix <- expr_matrix[top_genes_combined, ]

# Optional: scale rows (genes) for better contrast
scaled_matrix <- t(scale(t(heatmap_matrix)))
# Subset the expression matrix for these genes
heatmap_matrix <- expr_matrix[top_genes_combined, ]
# Scale each gene (row-wise) for better contrast
scaled_matrix <- t(scale(t(heatmap_matrix)))
# Sample annotation for OTA vs PTA
annotation_col <- data.frame(
  Condition = colData_phase$condition
)
rownames(annotation_col) <- rownames(colData_phase)
library(pheatmap)

pheatmap(
  scaled_matrix,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Expression Heatmap of Top PCA-Driving Genes",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)

#starting with early DEGs on July 25 2025
library(openxlsx)

# Load both sheets
deg_up <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_early_DEGs_up_down.xlsx", sheet = 1)
deg_down <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_early_DEGs_up_down.xlsx", sheet = 2)

# Combine
deg_all <- rbind(deg_up, deg_down)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

library(ggplot2)
library(ggrepel)

# Clean up just in case
deg_all$log2FoldChange <- as.numeric(deg_all$log2FoldChange)
deg_all$PValue <- as.numeric(deg_all$PValue)

ggplot(deg_all, aes(x = log2FoldChange, y = -log10(PValue))) +
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  theme_minimal() +
  labs(title = "Volcano Plot: All DEGs (No Cutoff)",
       x = "Log2 Fold Change", y = "-log10(P-value)")

top_genes <- deg_all[order(abs(deg_all$log2FoldChange), decreasing = TRUE), ][1:30, ]

ggplot(deg_all, aes(x = log2FoldChange, y = -log10(PValue))) +
  geom_point(alpha = 0.5, size = 1.5, color = "gray60") +
  geom_point(data = top_genes, aes(x = log2FoldChange, y = -log10(PValue)), color = "firebrick", size = 3) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = Inf) +
  theme_classic() +
  labs(title = "Volcano Plot: Highlighted Top 30 DEGs",
       x = "Log2 Fold Change", y = "-log10(P-value)")+ geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray70")


# Categorize regulation direction
deg_all$Category <- ifelse(deg_all$PValue < 0.05 & deg_all$log2FoldChange > 1, "Upregulated",
                           ifelse(deg_all$PValue < 0.05 & deg_all$log2FoldChange < -1, "Downregulated", "Not Significant"))

# Identify genes to label (top by abs log2FC or biological interest)
top_genes <- deg_all[order(abs(deg_all$log2FoldChange), decreasing = TRUE), ][1:30, ]







# Sort and extract top 50 by absolute log2FC
top_genes <- deg_all[order(abs(deg_all$log2FoldChange), decreasing = TRUE), ][1:50, ]





library(pheatmap)

# Extract expression columns (samples only)
sample_cols <- grep("Early_", names(deg_all))
expr_matrix <- as.matrix(deg_all[, sample_cols])
rownames(expr_matrix) <- deg_all$gene

# Z-score normalization (optional for clarity)
expr_scaled <- t(scale(t(expr_matrix)))

# Plot full heatmap
pheatmap(expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 6,
         fontsize_col = 8,
         show_rownames = FALSE,
         main = "Heatmap of All DEGs in Cacao Phase Change"
)

top_heat <- deg_all[deg_all$gene %in% top_genes$gene, ]
expr_top <- as.matrix(top_heat[, sample_cols])
rownames(expr_top) <- top_heat$gene

expr_top_scaled <- t(scale(t(expr_top)))

pheatmap(expr_top_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("green", "black", "magenta"))(100),
         fontsize_row = 7,
         show_rownames = TRUE,
         main = "Top 50 DEGs - Expression Heatmap"
)
write.xlsx(top_genes, file = "Top_50_DEGs_early.xlsx", sheetName = "TopGenes_early", rowNames = FALSE)

#try volcano
deg_up <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_early_DEGs_up_down.xlsx", sheet = 1)
deg_down <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_early_DEGs_up_down.xlsx", sheet = 2)

# Combine
deg_all <- rbind(deg_up, deg_down)
library(ggplot2)
library(ggrepel)
library(dplyr)



# Volcano plot: No FDR threshold
ggplot(deg_all, aes(x = log2FoldChange, y = -log10(FDR), label = gene)) +
  geom_point(aes(color = case_when(
    log2FoldChange > 0.1 ~ "Upregulated",
    log2FoldChange < -0.1 ~ "Downregulated",
    TRUE ~ "Neutral"
  )), size = 2) +
  
  geom_text_repel(aes(label = ifelse(abs(log2FoldChange) > 0.1, gene, "")),
                  size = 3, max.overlaps = 30) +
  
  scale_color_manual(values = c(
    "Upregulated" = "firebrick",
    "Downregulated" = "royalblue",
    "Neutral" = "gray"
  )) +
  
  theme_minimal(base_size = 14) +
  labs(x = "Log2 Fold Change",
       y = "-log10(FDR)",
       color = "Expression") +
  
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed")




#starting with late DEGs on July 25 2025
library(openxlsx)

# Load both sheets
deg_up_late <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_late_DEGs_up_down.xlsx", sheet = 1)
deg_down_late <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_late_DEGs_up_down.xlsx", sheet = 2)

# Combine
deg_all_late <- rbind(deg_up_late, deg_down_late)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)


# Clean up just in case
deg_all_late$log2FoldChange <- as.numeric(deg_all_late$log2FoldChange)
deg_all_late$PValue <- as.numeric(deg_all_late$PValue)



top_genes_late <- deg_all_late[order(abs(deg_all_late$log2FoldChange), decreasing = TRUE), ][1:50, ]

ggplot(deg_all_late, aes(x = log2FoldChange, y = -log10(PValue))) +
  geom_point(alpha = 0.5, size = 1.5, color = "gray60") +
  geom_point(data = top_genes, aes(x = log2FoldChange, y = -log10(PValue)), color = "firebrick", size = 3) +
  geom_text_repel(data = top_genes_late, aes(label = gene), size = 3, max.overlaps = Inf) +
  theme_classic() +
  labs(title = "Volcano Plot: Late Top 50 DEGs",
       x = "Log2 Fold Change", y = "-log10(P-value)")+ geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray70") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray70")









# Sort and extract top 50 by absolute log2FC
library(pheatmap)

# 🟣 Step 1: Select top 50 DEGs based on absolute log2 fold change
top_genes_late <- deg_all_late[order(abs(deg_all_late$log2FoldChange), decreasing = TRUE), ][1:50, ]

# 🟢 Step 2: Get sample expression columns (with 'Late_' prefix)
sample_cols <- grep("Late_", names(deg_all_late))

# 🔵 Step 3: Build full expression matrix and normalize (Z-score)
expr_matrix <- as.matrix(deg_all_late[, sample_cols])
rownames(expr_matrix) <- deg_all_late$gene
expr_scaled <- t(scale(t(expr_matrix)))

# 🔶 Step 4: Heatmap for all DEGs (optional overview)
pheatmap(expr_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 6,
         fontsize_col = 8,
         show_rownames = FALSE,
         main = "Heatmap of All DEGs in Late Cacao Phase Change"
)

# 🔺 Step 5: Subset expression matrix for top 50 DEGs only
top_heat_late <- deg_all_late[match(top_genes_late$gene, deg_all_late$gene), ]
expr_top_late <- as.matrix(top_heat_late[, sample_cols])
rownames(expr_top_late) <- top_heat_late$gene

# 🔻 Step 6: Normalize and plot targeted heatmap
expr_top_late_scaled <- t(scale(t(expr_top_late)))

pheatmap(expr_top_late_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("green", "black", "magenta"))(100),
         fontsize_row = 7,
         show_rownames = TRUE,
         main = "Top 50 Late DEGs - Expression Heatmap"
)
write.xlsx(top_genes_late, file = "Top_50_DEGs_late.xlsx", sheetName = "TopGenes_late", rowNames = FALSE)

#try volcano
deg_up_late <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_late_DEGs_up_down.xlsx", sheet = 1)
deg_down_late <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/top_late_DEGs_up_down.xlsx", sheet = 2)

# Combine
deg_all_late <- rbind(deg_up_late, deg_down_late)
library(ggplot2)
library(ggrepel)
library(dplyr)



# Volcano plot: No FDR threshold
ggplot(deg_all_late, aes(x = log2FoldChange, y = -log10(FDR), label = gene)) +
  geom_point(aes(color = case_when(
    log2FoldChange > 0.1 ~ "Upregulated",
    log2FoldChange < -0.1 ~ "Downregulated",
    TRUE ~ "Neutral"
  )), size = 2) +
  
  geom_text_repel(aes(label = ifelse(abs(log2FoldChange) > 0.1, gene, "")),
                  size = 3, max.overlaps = 30) +
  
  scale_color_manual(values = c(
    "Upregulated" = "firebrick",
    "Downregulated" = "royalblue",
    "Neutral" = "gray"
  )) +
  
  theme_minimal(base_size = 14) +
  labs(x = "Log2 Fold Change",
       y = "-log10(FDR)",
       color = "Expression") +
  
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed")


#try GO July 25
library(clusterProfiler)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)


library(clusterProfiler)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# 🟡 STEP 1: Load and clean GO annotation file
go_map <- read_csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/for R/Theobroma_cacao_SCA-6_chr.v1.0.summary.csv",
                   na = c("", "NA", "---NA---")) |>
  select(go_id = "GO IDs", gene_id = "SeqName") |>
  drop_na(go_id) |>
  mutate(
    gene_id = str_remove(gene_id, "\\.1$"),              # Match format: remove .1 suffix
    go_id = str_remove_all(go_id, "[FPC]:")              # Remove prefixes like F:, P:, C:
  ) |>
  separate_longer_delim(go_id, delim = "; ") |>
  mutate(
    gene_id = str_trim(gene_id),
    go_id = str_trim(go_id)
  )

# 🟢 STEP 2: Define your upregulated DEG vector manually
up_genes <- c(
  "SCA-6_Chr3v1_09040", "SCA-6_Chr3v1_08468", "SCA-6_Chr6v1_18722", "SCA-6_Chr1v1_01001",
  "SCA-6_Chr1v1_00531", "SCA-6_Chr1v1_01604", "SCA-6_Chr1v1_00398", "SCA-6_Chr1v1_00536",
  "SCA-6_Chr1v1_00373", "SCA-6_Chr1v1_01368", "SCA-6_Chr1v1_01622", "SCA-6_Chr1v1_00438",
  "SCA-6_Chr1v1_01654", "SCA-6_Chr1v1_00268", "SCA-6_Chr1v1_01649", "SCA-6_Chr1v1_00174",
  "SCA-6_Chr1v1_00666", "SCA-6_Chr1v1_00419", "SCA-6_Chr1v1_00014", "SCA-6_Chr1v1_00228",
  "SCA-6_Chr1v1_00881", "SCA-6_Chr1v1_01192", "SCA-6_Chr1v1_00278", "SCA-6_Chr1v1_00923",
  "SCA-6_Chr1v1_01066", "SCA-6_Chr1v1_01547", "SCA-6_Chr1v1_00420", "SCA-6_Chr1v1_00298",
  "SCA-6_Chr1v1_01019", "SCA-6_Chr1v1_00939", "SCA-6_Chr1v1_00964", "SCA-6_Chr1v1_01092",
  "SCA-6_Chr1v1_01088", "SCA-6_Chr1v1_01215", "SCA-6_Chr1v1_00572", "SCA-6_Chr1v1_00795",
  "SCA-6_Chr1v1_00231", "SCA-6_Chr1v1_00953", "SCA-6_Chr1v1_01330", "SCA-6_Chr1v1_01545",
  "SCA-6_Chr1v1_01408", "SCA-6_Chr1v1_00305", "SCA-6_Chr1v1_00323", "SCA-6_Chr1v1_00868",
  "SCA-6_Chr1v1_00197", "SCA-6_Chr1v1_01062", "SCA-6_Chr1v1_01221", "SCA-6_Chr1v1_00545",
  "SCA-6_Chr1v1_00025", "SCA-6_Chr1v1_00687", "SCA-6_Chr1v1_00861", "SCA-6_Chr1v1_01288",
  "SCA-6_Chr1v1_00169", "SCA-6_Chr1v1_01589", "SCA-6_Chr1v1_00093", "SCA-6_Chr1v1_00251",
  "SCA-6_Chr1v1_00325", "SCA-6_Chr1v1_01429", "SCA-6_Chr1v1_00049", "SCA-6_Chr1v1_01156",
  "SCA-6_Chr1v1_00022", "SCA-6_Chr1v1_01036", "SCA-6_Chr1v1_01239", "SCA-6_Chr1v1_00225",
  "SCA-6_Chr1v1_01068", "SCA-6_Chr1v1_01637", "SCA-6_Chr1v1_00530", "SCA-6_Chr1v1_00885",
  "SCA-6_Chr1v1_00125", "SCA-6_Chr1v1_01514", "SCA-6_Chr1v1_01530", "SCA-6_Chr1v1_01067",
  "SCA-6_Chr1v1_01597", "SCA-6_Chr1v1_00991", "SCA-6_Chr1v1_00189", "SCA-6_Chr1v1_00408",
  "SCA-6_Chr1v1_01158", "SCA-6_Chr1v1_00979", "SCA-6_Chr1v1_01327", "SCA-6_Chr1v1_00900",
  "SCA-6_Chr1v1_00166", "SCA-6_Chr1v1_00082", "SCA-6_Chr1v1_00678", "SCA-6_Chr1v1_00627",
  "SCA-6_Chr1v1_00734", "SCA-6_Chr1v1_00147"
)

# 🔵 STEP 3: Run GO enrichment
go_results_up <- enricher(
  gene = up_genes,
  TERM2GENE = go_map,
  pvalueCutoff = 0.1
)

# 🔴 STEP 4: Visualize
dotplot(go_results_up, showCategory = 100) +
  ggtitle("GO Enrichment of Upregulated Late DEGs")
library(ggplot2)

#for down degs late

library(clusterProfiler)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# 🟡 STEP 1: Load and clean GO annotation
go_map <- read_csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/for R/Theobroma_cacao_SCA-6_chr.v1.0.summary.csv",
                   na = c("", "NA", "---NA---")) |>
  select(go_id = "GO IDs", gene_id = "SeqName") |>
  drop_na(go_id) |>
  mutate(
    gene_id = str_remove(gene_id, "\\.1$"),
    go_id   = str_trim(go_id)
  ) |>
  separate_rows(go_id, sep = "; ") |>
  mutate(
    go_id = str_replace(go_id, "^[A-Z]:", ""),   # Remove P:, F:, C: prefixes
    gene_id = str_trim(gene_id),
    go_id   = str_trim(go_id)
  )



# 🧬 STEP 2: Manually define downregulated DEGs
down_genes <- c(
  "SCA-6_Chr9v1_24462", "SCA-6_Chr2v1_05915", "SCA-6_Chr1v1_00203", "SCA-6_Chr1v1_01347", "SCA-6_Chr1v1_01039", 
  "SCA-6_Chr1v1_00663", "SCA-6_Chr1v1_01567", "SCA-6_Chr1v1_00001", "SCA-6_Chr1v1_01124", "SCA-6_Chr1v1_01029", 
  "SCA-6_Chr1v1_01120", "SCA-6_Chr1v1_00114", "SCA-6_Chr1v1_00127", "SCA-6_Chr1v1_00226", "SCA-6_Chr1v1_00968", 
  "SCA-6_Chr1v1_00677", "SCA-6_Chr1v1_00444", "SCA-6_Chr1v1_01625", "SCA-6_Chr1v1_00966", "SCA-6_Chr1v1_00221", 
  "SCA-6_Chr1v1_00239", "SCA-6_Chr1v1_01149", "SCA-6_Chr1v1_00879", "SCA-6_Chr1v1_00110", "SCA-6_Chr1v1_00869", 
  "SCA-6_Chr1v1_00489", "SCA-6_Chr1v1_00422", "SCA-6_Chr1v1_01409", "SCA-6_Chr1v1_00651", "SCA-6_Chr1v1_00067", 
  "SCA-6_Chr1v1_00186", "SCA-6_Chr1v1_00367", "SCA-6_Chr1v1_00603", "SCA-6_Chr1v1_00024", "SCA-6_Chr1v1_00269", 
  "SCA-6_Chr1v1_01023", "SCA-6_Chr1v1_01437", "SCA-6_Chr1v1_00697", "SCA-6_Chr1v1_00179", "SCA-6_Chr1v1_01155", 
  "SCA-6_Chr1v1_00802", "SCA-6_Chr1v1_00158", "SCA-6_Chr1v1_00655", "SCA-6_Chr1v1_00046", "SCA-6_Chr1v1_01285", 
  "SCA-6_Chr1v1_00185", "SCA-6_Chr1v1_00759", "SCA-6_Chr1v1_00594", "SCA-6_Chr1v1_00393", "SCA-6_Chr1v1_00331", 
  "SCA-6_Chr1v1_00664", "SCA-6_Chr1v1_00113", "SCA-6_Chr1v1_01245", "SCA-6_Chr1v1_01397", "SCA-6_Chr1v1_01042", 
  "SCA-6_Chr1v1_00318", "SCA-6_Chr1v1_01582", "SCA-6_Chr1v1_00716", "SCA-6_Chr1v1_00707", "SCA-6_Chr1v1_00270", 
  "SCA-6_Chr1v1_00362", "SCA-6_Chr1v1_00755", "SCA-6_Chr1v1_01246", "SCA-6_Chr1v1_01552", "SCA-6_Chr1v1_01114", 
  "SCA-6_Chr1v1_01341", "SCA-6_Chr1v1_01147", "SCA-6_Chr1v1_00926", "SCA-6_Chr1v1_01171", "SCA-6_Chr1v1_01602", 
  "SCA-6_Chr1v1_00580", "SCA-6_Chr1v1_00495", "SCA-6_Chr1v1_01273", "SCA-6_Chr1v1_00322", "SCA-6_Chr1v1_01307", 
  "SCA-6_Chr1v1_01591", "SCA-6_Chr1v1_00299", "SCA-6_Chr1v1_01427", "SCA-6_Chr1v1_01615", "SCA-6_Chr1v1_01527", 
  "SCA-6_Chr1v1_00294", "SCA-6_Chr1v1_01332", "SCA-6_Chr1v1_01361", "SCA-6_Chr1v1_00295", "SCA-6_Chr1v1_01084", 
  "SCA-6_Chr1v1_00930", "SCA-6_Chr1v1_01151", "SCA-6_Chr1v1_00403", "SCA-6_Chr1v1_00951", "SCA-6_Chr1v1_01298", 
  "SCA-6_Chr1v1_00172", "SCA-6_Chr1v1_00810", "SCA-6_Chr1v1_01168", "SCA-6_Chr1v1_00234", "SCA-6_Chr1v1_01345", 
  "SCA-6_Chr1v1_00257", "SCA-6_Chr1v1_00656", "SCA-6_Chr1v1_00521", "SCA-6_Chr1v1_01618", "SCA-6_Chr1v1_00790"
)



# Extract top terms — even if not significant
top_terms <- go_results_down@result |>
  arrange(desc(Count)) |>
  slice_head(n = 10)

# Plot fallback
ggplot(top_terms, aes(x = reorder(Description, Count), y = Count)) +
  geom_col(fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top GO Terms in Downregulated DEGs",
    x = "GO Term Description",
    y = "Gene Count"
  )


#rerun deseq july 25
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
# Set your working directory
setwd("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/")



# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Set working directory
setwd("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/")

# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(EnhancedVolcano)

# Read the count data

countdata <- read.csv("deseq_ready_counts.csv", row.names=1, check.names=FALSE)
# The column names in countdata must match the sample IDs below, in order


# Create sample information
sample_names <- colnames(count_data)
condition <- ifelse(grepl("OTA", sample_names), "OTA", "PTA")
col_data <- data.frame(row.names = sample_names, condition = factor(condition))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)

# Pre-filtering to remove low count genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds)
resOrdered <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(resOrdered), file = "DESeq2_results.csv")

# PCA plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()


# Heatmap of top 40 DEGs
top_genes <- head(order(res$padj), 40)
pheatmap(assay(vsd)[top_genes, ],
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         annotation_col = col_data,
         main = "Top 40 Differentially Expressed Genes")


# Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Volcano Plot',
                subtitle = 'OTA vs PTA',
                caption = 'DESeq2 Analysis',
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)
# Create a custom color scheme for up/downregulated genes
keyvals <- ifelse(res$log2FoldChange < -1 & res$pvalue < 0.05, 'blue',
                  ifelse(res$log2FoldChange > 1 & res$pvalue < 0.05, 'red', 'grey'))
names(keyvals)[keyvals == 'red'] <- 'Upregulated'
names(keyvals)[keyvals == 'blue'] <- 'Downregulated'
names(keyvals)[keyvals == 'grey'] <- 'Not Significant'

# Volcano plot with smaller labels and custom colors
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Volcano Plot',
                subtitle = 'OTA vs PTA',
                caption = 'DESeq2 Analysis',
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                labSize = 2.5,  # smaller font size for gene labels
                colCustom = keyvals
)

# Load necessary library
library(readr)
library(dplyr)


# Read the cleaned DESeq2 results file
df <- read.csv("DESeq2_results.csv")

# Filter and count genes based on thresholds
n_padj_05 <- df %>% filter(padj < 0.05) %>% nrow()
n_padj_10 <- df %>% filter(padj < 0.1) %>% nrow()
n_padj_05_lfc_15 <- df %>% filter(padj < 0.05, abs(log2FoldChange) > 1.5) %>% nrow()
n_padj_10_lfc_15 <- df %>% filter(padj < 0.1, abs(log2FoldChange) > 1.5) %>% nrow()

# Print results
cat("Number of genes with padj < 0.05:", n_padj_05, "\n")
cat("Number of genes with padj < 0.1:", n_padj_10, "\n")
cat("Number of genes with padj < 0.05 and abs(log2FoldChange) > 1.5:", n_padj_05_lfc_15, "\n")
cat("Number of genes with padj < 0.1 and abs(log2FoldChange) > 1.5:", n_padj_10_lfc_15, "\n")

#rerun July 28 part 2 from original file

# Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

countdata <- read.csv("deseq_ready_counts.csv", row.names=1, check.names=FALSE)
# The column names in countdata must match the sample IDs below, in order
sample_ids <- c(
  "MA-OTA-Sca6-S1 001", "MA-OTA-Sca6-S1 002", "MA-OTA-Sca6-S1 003",
  "MA-OTA-Sca6-S2 001", "MA-OTA-Sca6-S2 002", "MA-OTA-Sca6-S2 003",
  "MA-OTA-Sca6-S3 001", "MA-OTA-Sca6-S3 002", "MA-OTA-Sca6-S3 003",
  "MA-OTA-Sca6-S4 001", "MA-OTA-Sca6-S4 002", "MA-OTA-Sca6-S4 003",
  "MA-PTA-Sca6-S1 001", "MA-PTA-Sca6-S1 002", "MA-PTA-Sca6-S1 003",
  "MA-PTA-Sca6-S2 001", "MA-PTA-Sca6-S2 002", "MA-PTA-Sca6-S2 003",
  "MA-PTA-Sca6-S3 001", "MA-PTA-Sca6-S3 002", "MA-PTA-Sca6-S3 003",
  "MA-PTA-Sca6-S4 001", "MA-PTA-Sca6-S4 002", "MA-PTA-Sca6-S4 003"
)
sample_ids <- colnames(countdata)
Group <- rep(LETTERS[1:8], each=3)
Type <- rep(c(rep("Ortho", 12), rep("Plagio", 12)))
Stage <- rep(rep(1:4, each=3), 2)

coldata <- data.frame(
  Group = Group,
  Type = Type,
  Stage = Stage,
  row.names = sample_ids  # set rownames explicitly to match your count columns
)


all(colnames(countdata) == rownames(coldata))  # Should return TRUE






dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~Type)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Type", "Plagio", "Ortho"))
resOrdered <- res[order(res$pvalue), ]
write.csv(as.data.frame(resOrdered), file="DEG_results.csv")

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="Type", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

topGenes <- head(order(res$padj), 50)
mat <- assay(vsd)[topGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=coldata)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)
# to parse
res <- results(dds, contrast = c("Type", "Plagio", "Ortho"))
# Remove NA padj values and filter for significant ones
sig_res <- res[!is.na(res$padj) & res$padj < 0.05, ]
write.csv(as.data.frame(sig_res), file = "significant_genes_padj_0.05.csv")

# Example: reading a TSV file named "counts.tsv" in your working directory
countdata_tsv <- read.delim("Theobroma_cacao_SCA-6_chr.v1.0.Meristem_Atlas.fractional.counts.tsv", row.names=1, check.names=FALSE)

# Use your actual colnames from the count matrix as sample IDs
sample_ids <- colnames(countdata_tsv)

# Assign groups (A to H, 3 replicates each)
Group <- rep(LETTERS[1:8], each=3)

# Assign Type based on group: A-D = Ortho, E-H = Plagio
Type <- ifelse(Group %in% c("A", "B", "C", "D"), "Ortho", "Plagio")

# Assign Stage 1 to 4 in triplicates for each group
Stage <- rep(rep(1:4, each=3), 2)  # Because you have 8 groups total (4 stages x 2 types)

# Create metadata dataframe with rownames as sample IDs
coldata <- data.frame(
  Group = Group,
  Type = factor(Type),   # factor is recommended for DESeq2
  Stage = factor(Stage),
  row.names = sample_ids
)
all(rownames(coldata) == colnames(countdata_tsv))  # Should return TRUE
countdata_int <- round(countdata_tsv)

# Check if all values are non-negative integers
if(all(countdata_int >= 0) && all(countdata_int == floor(countdata_int))) {
  dds <- DESeqDataSetFromMatrix(countData = countdata_int,
                                colData = coldata,
                                design = ~ Type)
} else {
  stop("Count data contains negative or non-integer values even after rounding")
}

dds <- DESeqDataSetFromMatrix(countData = countdata_int,
                              colData = coldata,
                              design = ~ Type)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Type", "Plagio", "Ortho"))
write.csv(as.data.frame(res), file = "DESeq2_results_Plagio_vs_Ortho.csv")

# Filter for genes with padj < 0.05 and remove NA padj values
sig_genes <- res[!is.na(res$padj) & res$padj < 0.05, ]

# Convert to data frame for saving
sig_genes_df <- as.data.frame(sig_genes)
head(sig_genes_df[order(sig_genes_df$padj), ])


# Write filtered results to CSV
write.csv(sig_genes_df, file = "DESeq2_significant_genes_padj_0.05.csv", row.names = TRUE)

#to compare with apical bud file

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

library(readxl)

# Read the "Results" sheet (sheet 3) from your Excel file
original_res <- read_excel("PhaseChange_ApicalBudDESeq.xlsx", sheet = "Results")

# Inspect the first few rows and column names to confirm structure
head(original_res)
colnames(original_res)
library(readxl)

# 1. Read the 'Results' sheet again (already done)
original_res <- read_excel("PhaseChange_ApicalBudDESeq.xlsx", sheet = "Results")

# 2. Convert to data.frame and check the first few rows
original_res <- as.data.frame(original_res)
head(original_res)

# 3. Rename the gene ID column from ...1 to "GeneID"
colnames(original_res)[1] <- "GeneID"

# 4. Set row names as GeneID and optionally remove the GeneID column
rownames(original_res) <- original_res$GeneID
original_res$GeneID <- NULL

# 5. Examine column names again to select columns related to DESeq2 results only
# For your comparison, focus on these columns only:
columns_of_interest <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

# 6. Subset original_res to only these columns for cleaner comparison
original_res_sub <- original_res[, columns_of_interest]

# 7. Prepare your DESeq2 results similarly if not done already
my_res_df <- as.data.frame(res)  # 'res' is your DESeq2 results object

# 8. Identify common genes between your results and original
common_genes <- intersect(rownames(my_res_df), rownames(original_res_sub))

# 9. Subset both datasets by common genes
my_res_sub <- my_res_df[common_genes, ]
original_res_sub <- original_res_sub[common_genes, ]

# 10. Rename columns for clarity
colnames(my_res_sub) <- paste0("My_", colnames(my_res_sub))
colnames(original_res_sub) <- paste0("Orig_", colnames(original_res_sub))

# 11. Combine side-by-side
comparison_df <- cbind(my_res_sub, original_res_sub)

# 12. Calculate differences for key numerical columns to highlight discrepancies
comparison_df$log2FC_diff <- comparison_df$My_log2FoldChange - comparison_df$Orig_log2FoldChange
comparison_df$padj_diff <- comparison_df$My_padj - comparison_df$Orig_padj

# View genes with largest discrepancies in log2 fold change
head(comparison_df[order(abs(comparison_df$log2FC_diff), decreasing = TRUE), ])

# 13. Optionally save comparison to CSV for detailed review
write.csv(comparison_df, file = "DESeq2_vs_Original_Results_Comparison.csv", row.names = TRUE)

#rerun
dds <- DESeqDataSetFromMatrix(countData = countdata_int,
                              colData = coldata,
                              design = ~ Stage + Type)

dds <- DESeq(dds)
res <- results(dds, contrast = c("Type", "Plagio", "Ortho"))
write.csv(as.data.frame(res), file = "DESeq2_results_Plagio_vs_Ortho_2.csv")
library(DESeq2)
library(ggplot2)

# Transformation (choose vst or rlog, vst is faster)
vsd <- vst(dds, blind=FALSE)

# Extract PCA data for plotting
pcaData <- plotPCA(vsd, intgroup="Type", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot with ggplot2
ggplot(pcaData, aes(PC1, PC2, color=Type)) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  theme_minimal()
library(pheatmap)

# Choose top 50 genes by smallest padj
topGenes <- head(order(res$padj), 50)

# Get variance stabilized counts for those genes
mat <- assay(vsd)[topGenes, ]

# Center by row (gene)
mat <- mat - rowMeans(mat)

# Annotation for samples (use colData from dds)
annotation_col <- as.data.frame(colData(dds)[, c("Type", "Stage")])

# Plot heatmap
pheatmap(mat, annotation_col = annotation_col)



# Append regulation status factor column to the res object (your DESeq2 result data frame)
threshold_log2FC <- 1
threshold_padj <- 0.05

res$regulation_status <- "Not Significant"
res$regulation_status[res$log2FoldChange > threshold_log2FC & res$padj < threshold_padj] <- "Upregulated"
res$regulation_status[res$log2FoldChange < -threshold_log2FC & res$padj < threshold_padj] <- "Downregulated"
res$regulation_status <- factor(res$regulation_status, levels = c("Upregulated", "Downregulated", "Not Significant"))

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = rownames(res)[res$regulation_status != "Not Significant"],
                labSize = 3,
                pointSize = 2.5,
                pCutoff = threshold_padj,
                FCcutoff = threshold_log2FC,
                colCustom = c('Not Significant' = 'grey30', 'Upregulated' = 'red', 'Downregulated' = 'blue'),
                legendPosition = 'right',
                drawConnectors = TRUE,
                max.overlaps = 30
)


#deseq rerun July 29
# Set working directory
setwd("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/")

# Load libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)

# Load count data
counts <- read.csv("deseq_ready_counts.csv", row.names = 1, check.names = FALSE)

# Remove the outlier sample: MA-PTA-Sca6-S2_R1
counts <- counts[, !grepl("MA-PTA-Sca6-S2_R1", colnames(counts))]

# Create sample metadata (colData)
sample_names <- colnames(counts)
condition <- ifelse(grepl("MA-OTA", sample_names), "Ortho", "Plagio")
stage <- sub(".*Sca6-(S[1-4])_.*", "\\1", sample_names)
replicate <- sub(".*_R([1-3])", "\\1", sample_names)

colData <- data.frame(
  row.names = sample_names,
  condition = factor(condition, levels = c("Ortho", "Plagio")),
  stage = factor(stage),
  replicate = factor(replicate)
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Get results for Plagio vs Ortho
res <- results(dds, contrast = c("condition", "Plagio", "Ortho"))
write.csv(as.data.frame(res), "DEG_Plagio_vs_Ortho_results_july29.csv")

# Filter significant genes
res_sig_0.2 <- res[which(res$padj < 0.3 & abs(res$log2FoldChange) >= 1), ]
res_sig_0.2 <- res_sig[order(res_sig$padj), ]

# Save results
write.csv(as.data.frame(res_sig), "DEG_Plagio_vs_Ortho_0.2cutoff.csv")

res_sig <- res[which(res$padj < 0.1 & abs(res$log2FoldChange) >= 1), ]
res_sig <- res_sig[order(res_sig$padj), ]

# Save results
write.csv(as.data.frame(res_sig), "DEG_Plagio_vs_Ortho.csv")

# Heatmap of top 50 DEGs
vsd <- vst(dds, blind = FALSE)
top_genes <- head(rownames(res_sig), 50)
pheatmap(assay(vsd)[top_genes, ], cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, annotation_col = colData)

# Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.1,
                FCcutoff = 1,
                title = 'Plagiotropic vs Orthotropic Apical Buds',
                subtitle = 'DESeq2 results',
                caption = 'padj < 0.1 & |log2FC| ≥ 1')


#using the file from mg July 29
# Set working directory
setwd("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/Biodiversity/patrick rerun/")

# Load required libraries
library(readxl)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(DESeq2)

# Load DESeq2 results
deg <- read_excel("DESeq2_results_Plagio_vs_Ortho_20.2 padj cutoff.xlsx")
colnames(deg)[1] <- "GeneID"

EnhancedVolcano(deg,
                lab = deg$GeneID,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.2,
                FCcutoff = 1,
                title = 'Plagiotropic vs Orthotropic Apical Buds',
                subtitle = 'Significant DEGs (padj < 0.2)',
                caption = 'No additional filtering applied',
                pointSize = 3.0,
                labSize = 3.5)


