install.packages("igraph")  # Install if missing
install.packages("visNetwork")
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(clusterProfiler)
library(ggplot2)
library(igraph)  
# Load the expression matrix (contains Cacao ID, Arabidopsis ID, gene names, and expression values)
expr_data <- read_excel("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2025/Cacao_Embryo_GeneExpression.xlsx")# Create matrix of expression values only
data_matrix <- expr_data[, 4:ncol(expr_data)] |> as.matrix()
rownames(data_matrix) <- expr_data$`Name/Description`


# Z-score normalization per gene (row)
data_scaled <- t(scale(t(data_matrix)))
cacao_lookup <- expr_data |>
  select(`Name/Description`, `Cacao ID`) |>
  mutate(`Cacao ID` = str_remove(`Cacao ID`, "\\.1$"))

go_map <- read_csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/for R/Theobroma_cacao_SCA-6_chr.v1.0.summary.csv",
                   na = c("", "NA", "---NA---")) |>
  select(go_id = "GO IDs", gene_id = "SeqName") |>
  drop_na(go_id) |>
  mutate(gene_id = str_remove(gene_id, "\\.1$"),
         go_id = str_remove_all(go_id, "[FPC]:")) |>
  separate_longer_delim(go_id, delim = "; ")

# Match gene names in expression matrix to Cacao IDs
target_genes <- cacao_lookup |>
  filter(`Name/Description` %in% rownames(data_scaled)) |>
  pull(`Cacao ID`)

go_results <- enricher(gene = target_genes,
                       TERM2GENE = go_map,
                       pvalueCutoff = 0.05)

# Check result table
head(go_results)

if (!is.null(go_results) && nrow(go_results@result) > 0) {
  barplot(go_results, showCategory = 20, title = "GO Enrichment")
} else {
  message("No enriched GO terms found — consider checking input gene IDs or adjusting cutoff.")
}

#by functions
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Load custom GO mapping file
go_map <- read_csv("C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/for R/Theobroma_cacao_SCA-6_chr.v1.0.summary.csv",
                   na = c("", "NA", "---NA---")) %>%
  select(go_id = "GO IDs", gene_id = "SeqName") %>%
  drop_na(go_id) %>%
  mutate(gene_id = str_remove(gene_id, "\\.1$")) %>%
  separate_longer_delim(go_id, delim = "; ") %>%
  mutate(category = case_when(
    str_starts(go_id, "P:") ~ "BP",
    str_starts(go_id, "F:") ~ "MF",
    str_starts(go_id, "C:") ~ "CC"
  )) %>%
  mutate(go_id = str_remove(go_id, "^[PFC]:"))  # Remove category prefix
library(clusterProfiler)

# Biological Process (BP)
go_bp <- enricher(gene = target_genes, TERM2GENE = go_map %>% filter(category == "BP"), pvalueCutoff = 0.05)

# Molecular Function (MF)
go_mf <- enricher(gene = target_genes, TERM2GENE = go_map %>% filter(category == "MF"), pvalueCutoff = 0.05)

# Cellular Component (CC)
go_cc <- enricher(gene = target_genes, TERM2GENE = go_map %>% filter(category == "CC"), pvalueCutoff = 0.05)
# Biological Process (BP)
barplot(go_bp, showCategory = 15, title = "GO Enrichment (Biological Processes)")

# Molecular Function (MF)
barplot(go_mf, showCategory = 15, title = "GO Enrichment (Molecular Functions)")

# Cellular Component (CC)
barplot(go_cc, showCategory = 15, title = "GO Enrichment (Cellular Components)")

# Compute gene expression correlation
cor_matrix <- cor(t(data_matrix))

# Construct network graph
adj_matrix <- cor_matrix > 0.8  # Adjust threshold if necessary
adj_matrix[is.na(adj_matrix)] <- FALSE  # Remove NA values
gene_network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

# Visualize network structure
# Apply force-directed layout
layout <- layout_with_fr(gene_network)

plot(gene_network, layout = layout, 
     vertex.label = rownames(cor_matrix), 
     vertex.size = 8, 
     vertex.label.cex = 0.6, # Reduce font size
     edge.width = 0.5,       # Make edges thinner
     edge.color = "gray50")  # Use lighter edge colors for visibility

#alt layout
plot(gene_network, layout = layout_with_kk(gene_network),  # Kamada-Kawai layout
     vertex.label = rownames(cor_matrix), 
     vertex.size = 8, 
     vertex.label.cex = 0.5, # Smaller labels to prevent clutter
     vertex.label.color = "black", 
     edge.width = 0.5)

library(visNetwork)

visNetwork(toVisNetworkData(gene_network)$nodes, toVisNetworkData(gene_network)$edges) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
library(visNetwork)



layout <- layout_with_kk(gene_network)


#now looking into a comparative analysis
colnames(expr_data)
# Extract numerical count data (skip metadata columns)
count_data <- expr_data[, c("CCN51 Immature Fruit Embryo", 
                            "CCN51 Developing Fruit Embryo", 
                            "CCN51 Mature Embryo", 
                            "CCN51 Germinating Seed Shoot")]

count_data <- as.data.frame(count_data)  # Convert tibble to data frame
rownames(count_data) <- expr_data$`Cacao ID`  # Now assign row names properly

meta$Stage <- factor(meta$Stage, levels = c("Immature", "Developing", "Mature", "Germinating"))

# Make "Immature" the reference level
meta$Stage <- relevel(meta$Stage, ref = "Immature")

# Create design matrix
design <- model.matrix(~ Stage, data = meta)

# Check column names
colnames(design)
fit <- glmFit(dge, design)

# Perform pairwise comparisons between embryo stages
res_imm_dev <- glmLRT(fit, contrast = c(-1, 1, 0, 0))  # Immature vs Developing
res_dev_mat <- glmLRT(fit, contrast = c(0, -1, 1, 0))  # Developing vs Mature
res_mat_germ <- glmLRT(fit, contrast = c(0, 0, -1, 1)) # Mature vs Germinating
sig_imm_dev <- topTags(res_imm_dev, n = Inf)$table
sig_dev_mat <- topTags(res_dev_mat, n = Inf)$table
sig_mat_germ <- topTags(res_mat_germ, n = Inf)$table


# Filter significant genes with adjusted p-value < 0.05
sig_imm_dev <- sig_imm_dev[sig_imm_dev$FDR < 0.05, ]
sig_dev_mat <- sig_dev_mat[sig_dev_mat$FDR < 0.05, ]
sig_mat_germ <- sig_mat_germ[sig_mat_germ$FDR < 0.05, ]

# Move gene descriptions from rownames to a new column
sig_imm_dev$GeneName <- rownames(sig_imm_dev)
sig_dev_mat$GeneName <- rownames(sig_dev_mat)
sig_mat_germ$GeneName <- rownames(sig_mat_germ)

library(openxlsx)

# Create a named list of data frames
sig_gene_sheets <- list(
  Immature_to_Developing = sig_imm_dev,
  Developing_to_Mature = sig_dev_mat,
  Mature_to_Germinating = sig_mat_germ
)

# Save to an Excel workbook
write.xlsx(sig_gene_sheets, file = "Significant_Genes_by_Transition.xlsx", rowNames = FALSE)



#to plot them
library(ggplot2)
library(dplyr)

#using LFC
# Create a combined data frame for all transitions
df <- data.frame(
  GeneName = c(
    # Immature → Developing
    "2S seed storage protein 5", "Sorting nexin 1", "Lysine-specific demethylase REF6", 
    "Glucosamine inositolphosphorylceramide transferase 1", "Peroxisome biogenesis protein 16", 
    "Ureide permease 2", "Cellulose synthase A catalytic subunit 3", "BURP domain protein RD22", 
    "3-ketoacyl-CoA thiolase 2, peroxisomal", "EM6; ABA-responsive LEA",
    
    # Developing → Mature
    "Vicilin-like seed storage protein", "2S seed storage protein 5", "3-ketoacyl-CoA thiolase 2, peroxisomal",
    "Ureide permease 2", "Polygalacturonase 1 beta-like protein 3", "Sec1 family domain-containing protein MIP3",
    "Peroxisome biogenesis protein 16", "BURP domain protein RD22", "Lysine-specific demethylase REF6",
    
    # Mature → Germinating
    "2S seed storage protein 5", "Sorting nexin 1", "Lysine-specific demethylase REF6", 
    "Glucosamine inositolphosphorylceramide transferase 1", "Peroxisome biogenesis protein 16", 
    "Ureide permease 2", "Cellulose synthase A catalytic subunit 3", "BURP domain protein RD22", 
    "3-ketoacyl-CoA thiolase 2, peroxisomal", "EM6; ABA-responsive LEA"
  ),
  logFC = c(
    # Immature → Developing
    3.7308, -3.5697, -3.0723, -3.9851, -2.7487, 5.6560, -1.9619, 3.4852, 1.8394, 3.0765,
    
    # Developing → Mature
    14.5177, 10.0769, -2.5972, -6.2728, -2.9938, -3.1879, 2.5836, -3.7463, 1.9346,
    
    # Mature → Germinating
    3.7308, -3.5697, -3.0723, -3.9851, -2.7487, 5.6560, -1.9619, 3.4852, 1.8394, 3.0765
  ),
  Transition = c(
    rep("Immature → Developing", 10),
    rep("Developing → Mature", 9),
    rep("Mature → Germinating", 10)
  )
)

# Create a unique identifier for each gene-transition pair
df$GeneTransition <- paste(df$GeneName, df$Transition, sep = " | ")

# Convert GeneTransition into a factor with proper ordering
df$GeneTransition <- factor(df$GeneTransition, levels = df$GeneTransition[order(df$logFC)])

ggplot(df, aes(x = GeneTransition, y = logFC, fill = Transition)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = TRUE) +
  coord_flip() +  # Flip for better readability
  labs(title = "Differential Expression Across Embryo Transitions",
       x = "Gene Name | Transition", y = "Log Fold Change (logFC)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors for transitions

#using LR
library(ggplot2)
library(dplyr)

# Create a combined data frame for all transitions
df_2 <- data.frame(
  GeneName = c(
    # Immature → Developing
    "2S seed storage protein 5", "Sorting nexin 1", "Lysine-specific demethylase REF6", 
    "Glucosamine inositolphosphorylceramide transferase 1", "Peroxisome biogenesis protein 16", 
    "Ureide permease 2", "Cellulose synthase A catalytic subunit 3", "BURP domain protein RD22", 
    "3-ketoacyl-CoA thiolase 2, peroxisomal", "EM6; ABA-responsive LEA",
    
    # Developing → Mature
    "Vicilin-like seed storage protein", "2S seed storage protein 5", "3-ketoacyl-CoA thiolase 2, peroxisomal",
    "Ureide permease 2", "Polygalacturonase 1 beta-like protein 3", "Sec1 family domain-containing protein MIP3",
    "Peroxisome biogenesis protein 16", "BURP domain protein RD22", "Lysine-specific demethylase REF6",
    
    # Mature → Germinating
    "2S seed storage protein 5", "Sorting nexin 1", "Lysine-specific demethylase REF6", 
    "Glucosamine inositolphosphorylceramide transferase 1", "Peroxisome biogenesis protein 16", 
    "Ureide permease 2", "Cellulose synthase A catalytic subunit 3", "BURP domain protein RD22", 
    "3-ketoacyl-CoA thiolase 2, peroxisomal", "EM6; ABA-responsive LEA"
  ),
  LR = c(
    # Immature → Developing
    26.50, 21.77, 18.62, 17.65, 12.04, 9.23, 8.07, 7.76, 7.22, 5.63,
    
    # Developing → Mature
    96.80, 38.73, 13.70, 11.33, 11.24, 11.00, 9.13, 8.08, 8.00,
    
    # Mature → Germinating
    26.50, 21.77, 18.62, 17.65, 12.04, 9.23, 8.07, 7.76, 7.22, 5.63
  ),
  Transition = c(
    rep("Immature → Developing", 10),
    rep("Developing → Mature", 9),
    rep("Mature → Germinating", 10)
  )
)

# Create a unique identifier for each gene-transition pair
df_2$GeneTransition <- paste(df_2$GeneName, df_2$Transition, sep = " | ")

# Convert GeneTransition into a factor with proper ordering
df_2$GeneTransition <- factor(df_2$GeneTransition, levels = df_2$GeneTransition[order(df_2$LR, decreasing = TRUE)])
ggplot(df_2, aes(x = GeneTransition, y = LR, fill = Transition)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = TRUE) +
  coord_flip() +  # Flip for better readability
  labs(title = "Likelihood Ratio (LR) Across Embryo Transitions",
       x = "Gene Name | Transition", y = "Likelihood Ratio (LR)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors for transitions


library(ggplot2)

#Volcano Plot for Immature → Developing
library(ggplot2)
library(ggrepel)

ggplot(sig_imm_dev, aes(x = logFC, y = -log10(FDR), label = GeneID)) +
  geom_point(aes(color = case_when(
    FDR < 0.05 & logFC > 1 ~ "Upregulated",
    FDR < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )), size = 2) +
  geom_text_repel(aes(label = ifelse(FDR < 0.05, GeneID, "")), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Immature → Developing",
       x = "Log2 Fold Change",
       y = "-log10(FDR)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

#Developing → Mature
ggplot(sig_dev_mat, aes(x = logFC, y = -log10(FDR), label = GeneID)) +
  geom_point(aes(color = case_when(
    FDR < 0.05 & logFC > 1 ~ "Upregulated",
    FDR < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )), size = 2) +
  geom_text_repel(aes(label = ifelse(FDR < 0.05, GeneID, "")), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Developing → Mature",
       x = "Log2 Fold Change",
       y = "-log10(FDR)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")

#Mature → Germinating
ggplot(sig_mat_germ, aes(x = logFC, y = -log10(FDR), label = GeneID)) +
  geom_point(aes(color = case_when(
    FDR < 0.05 & logFC > 1 ~ "Upregulated",
    FDR < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )), size = 2) +
  geom_text_repel(aes(label = ifelse(FDR < 0.05, GeneID, "")), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Mature → Germinating",
       x = "Log2 Fold Change",
       y = "-log10(FDR)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")


library(pheatmap)

# Standardize expression data to Z-scores
expr_matrix_scaled <- t(scale(t(expr_matrix_filtered)))

# Create heatmap
pheatmap(expr_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 6,  # Adjust font size for readability
         fontsize_col = 10,
         show_rownames = FALSE,  # Hide row names if too cluttered
         main = "Differential Gene Expression Heatmap")





#GO enrichment
library(readr)

# Load annotation file
file_path <- "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/for R/Theobroma_cacao_SCA-6_chr.v1.0.summary.csv"
go_custom <- read.csv(file_path, check.names = FALSE)

# Check structure
str(go_custom)
head(go_custom)

library(tidyr)

# Extract cacao gene IDs and GO terms
go_map <- go_custom[, c("SeqName", "GO IDs")]

# Rename columns for compatibility
colnames(go_map) <- c("CacaoGeneID", "GO")

# Split multiple GO terms into separate rows
go_map <- separate_rows(go_map, GO, sep = ";")

# Clean leading spaces and prefixes
go_map$GO <- trimws(go_map$GO)  # Remove spaces
go_map$GO <- gsub("^[PFC]:", "", go_map$GO)  # Remove category prefixes

# Check structure
head(go_map)

go_map <- as.data.frame(go_map)

# Define gene sets for each transition
genes_imm_dev <- c(
  "SCA-6_Chr4v1_12901.1", "SCA-6_Chr7v1_19601.1", "SCA-6_Chr4v1_12655.1",
  "SCA-6_Chr2v1_07230.1", "SCA-6_Chr2v1_06511.1", "SCA-6_Chr1v1_01433.1",
  "SCA-6_Chr1v1_03538.1", "SCA-6_Chr9v1_22983.1", "SCA-6_Chr3v1_08609.1"
)

genes_dev_mat <- c(
  "SCA-6_Chr7v1_19601.1", "SCA-6_Chr5v1_14510.1", "SCA-6_Chr3v1_08609.1",
  "SCA-6_Chr10v1_27396.1", "SCA-6_Chr1v1_03538.1", "SCA-6_Chr2v1_07230.1",
  "SCA-6_Chr1v1_03934.1", "SCA-6_Chr9v1_22983.1", "SCA-6_Chr4v1_12655.1",
  "SCA-6_Chr10v1_26675.1"
)

genes_mat_germ <- c(
  "SCA-6_Chr4v1_12901.1", "SCA-6_Chr2v1_06511.1", "SCA-6_Chr1v1_03934.1",
  "SCA-6_Chr1v1_01433.1", "SCA-6_Chr5v1_14510.1", "SCA-6_Chr9v1_23744.1",
  "SCA-6_Chr4v1_13592.1", "SCA-6_Chr1v1_02584.1", "SCA-6_Chr1v1_00141.1",
  "SCA-6_Chr3v1_08609.1", "SCA-6_Chr3v1_08281.1", "SCA-6_Chr2v1_07230.1",
  "SCA-6_Chr8v1_21218.1", "SCA-6_Chr3v1_08238.1", "SCA-6_Chr9v1_22983.1",
  "SCA-6_Chr4v1_12580.1", "SCA-6_Chr10v1_27396.1", "SCA-6_Chr1v1_01786.1",
  "SCA-6_Chr1v1_03538.1", "SCA-6_Chr5v1_13788.1", "SCA-6_Chr4v1_12120.1",
  "SCA-6_Chr4v1_12655.1"
)

library(clusterProfiler)

genes_imm_dev <- as.character(genes_imm_dev)
genes_dev_mat <- as.character(genes_dev_mat)
genes_mat_germ <- as.character(genes_mat_germ)

library(clusterProfiler)

library(clusterProfiler)

# Ensure TERM2GENE format: two columns (Gene ID, GO Term)
term2gene <- go_map[, c("GO", "CacaoGeneID")]

# Check structure
head(term2gene)

# Run GO enrichment for each transition
go_imm_dev <- enricher(genes_imm_dev, TERM2GENE = term2gene, pvalueCutoff = 0.05)
go_dev_mat <- enricher(genes_dev_mat, TERM2GENE = term2gene, pvalueCutoff = 0.05)
go_mat_germ <- enricher(genes_mat_germ, TERM2GENE = term2gene, pvalueCutoff = 0.05)

# Visualize results
dotplot(go_imm_dev, title = "GO Enrichment (Immature → Developing)")
dotplot(go_dev_mat, title = "GO Enrichment (Developing → Mature)")
dotplot(go_mat_germ, title = "GO Enrichment (Mature → Germinating)")


#combined plot
library(dplyr)

# Convert enrichment results to data frames
df_imm_dev <- as.data.frame(go_imm_dev)
df_dev_mat <- as.data.frame(go_dev_mat)
df_mat_germ <- as.data.frame(go_mat_germ)
# Create workbook
wb <- createWorkbook()

# Add sheets
addWorksheet(wb, "Immature to Developing")
addWorksheet(wb, "Developing to Mature")
addWorksheet(wb, "Mature to Germinating")

# Write data
writeData(wb, "Immature to Developing", df_imm_dev)
writeData(wb, "Developing to Mature", df_dev_mat)
writeData(wb, "Mature to Germinating", df_mat_germ)

# Save file
saveWorkbook(wb, file = "GO_Enrichment_Results.xlsx", overwrite = TRUE)

# Add transition labels
df_imm_dev$Transition <- "Immature → Developing"
df_dev_mat$Transition <- "Developing → Mature"
df_mat_germ$Transition <- "Mature → Germinating"

# Combine into one dataset
df_combined <- bind_rows(df_imm_dev, df_dev_mat, df_mat_germ)

library(ggplot2)

ggplot(df_combined, aes(x = Description, y = -log10(p.adjust), color = Transition)) +
  geom_point(size = 3) +
  coord_flip() +
  facet_wrap(~ Transition, scales = "free_y") +
  labs(title = "GO Enrichment Across Embryo Transitions",
       x = "GO Term", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors for transitions


library(ggplot2)

ggplot(df_combined, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust), fill = Transition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), show.legend = TRUE) +
  coord_flip() +
  labs(title = "GO Enrichment Across Embryo Transitions",
       x = "GO Term", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors for transitions


ggplot(df_combined, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust), fill = Transition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), show.legend = TRUE) +
  coord_flip() +
  labs(title = "GO Enrichment Across Embryo Transitions",
       x = "GO Term", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),  # Reduce font size for terms
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors


ggplot(df_combined, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust), fill = Transition)) +
  geom_bar(stat = "identity") +
  facet_grid(Transition ~ ., scales = "free_y", space = "free_y") +  # Creates separate panels with free spacing
  coord_flip() +
  labs(title = "GO Enrichment Across Embryo Transitions",
       x = "GO Term", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors


ggplot(df_combined, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust), fill = Transition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1.2), show.legend = TRUE) +  # Increase bar spacing
  coord_flip() +
  labs(title = "GO Enrichment Across Embryo Transitions",
       x = "GO Term", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),  # Further reduce font size
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c"))  # Custom colors


# Prepare GO enrichment results
go_data <- data.frame(Gene = df_combined$Description,
                      Transition = df_combined$Transition,
                      LogPval = -log10(df_combined$p.adjust))

# Convert to matrix format
library(reshape2)
go_matrix <- acast(go_data, Gene ~ Transition, value.var = "LogPval")
go_matrix[is.na(go_matrix)] <- min(go_matrix, na.rm = TRUE) / 2  # Replace NA with half of min value
go_matrix_scaled <- apply(go_matrix, 2, function(x) (x - min(x)) / (max(x) - min(x)))
# Create heatmap
library(pheatmap)

pheatmap(go_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "orange", "darkred"))(100),
         fontsize_row = 4,  # Reduce font size for GO terms
         fontsize_col = 8,  # Slightly reduce column labels
         main = "GO Enrichment Heatmap")


#co-expression analysis
library(DESeq2)

# Apply variance-stabilizing transformation
vsd <- varianceStabilizingTransformation(dds)

# Extract normalized expression matrix
expr_matrix <- assay(vsd)

library(factoextra)
dist_matrix <- dist(t(expr_matrix_filtered), method = "euclidean")
hc_res <- hclust(dist_matrix, method = "ward.D2")

# Plot dendrogram
plot(hc_res, main = "Hierarchical Clustering of Expression Data")

