# Code under work for this project as of October 2025

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




# Define file path
raw_phase <- read.csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/deseq2_OTA_vs_PTA.csv")

#read other files

expr_data_early<-read.csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/deseq2_OTA_early_vs_PTA_early.csv")
expr_data_late<-read.csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/deseq2_OTA_late_vs_PTA_late.csv")

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
           "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/filtered_rawP.xlsx")

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
output_path <- "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/top_DEGs_up_down.xlsx"

write_xlsx(
  list(
    Upregulated = upregulated,
    Downregulated = downregulated
  ),
  path = output_path
)









 


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


#continued on july 23 to create heatmaps and volcano plots

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
  "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/top_DEGs_up_down.xlsx")

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


#Draw Heatmap —#
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


# PCA analysis

# PC1 vs PC2 plot
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

# Extract PCA loadings (rotation matrix)
loadings <- pca_res$rotation  # dimensions: genes × PCs

# Top 20 genes contributing to PC1
top_PC1_genes <- sort(abs(loadings[, 1]), decreasing = TRUE)[1:20]
top_PC1_names <- names(top_PC1_genes)

# Top 20 genes contributing to PC2
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
deg_up <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/top_early_DEGs_up_down.xlsx", sheet = 1)
deg_down <- read.xlsx("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/top_early_DEGs_up_down.xlsx", sheet = 2)

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


# GO  enrichment July 25
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

# STEP 1: Load and clean GO annotation file
go_map <- read_csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/Theobroma_cacao_SCA-6_chr.v1.0.summary.csv",
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

# STEP 2: Define your upregulated DEG vector manually
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

# STEP 3: Run GO enrichment
go_results_up <- enricher(
  gene = up_genes,
  TERM2GENE = go_map,
  pvalueCutoff = 0.1
)

# STEP 4: Visualize
dotplot(go_results_up, showCategory = 100) +
  ggtitle("GO Enrichment of Upregulated Late DEGs")
library(ggplot2)

# repeat for down degs late





