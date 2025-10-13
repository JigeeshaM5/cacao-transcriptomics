# run_gsea.R
library(clusterProfiler)
library(org.At.tair.db)  # Replace with appropriate annotation package for cacao
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read DESeq2 results
deg <- read_csv(input_file)

# Prepare ranked gene list
gene_list <- deg %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange)) %>%
  select(Gene_ID, log2FoldChange)

ranks <- gene_list$log2FoldChange
names(ranks) <- gene_list$Gene_ID

# Run GSEA (example using GO Biological Process)
gsea_res <- gseGO(
  geneList = ranks,
  OrgDb = org.At.tair.db,  # Replace with cacao annotation
  ont = "BP",
  keyType = "SYMBOL",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Save results
write.csv(as.data.frame(gsea_res), output_file, row.names = FALSE)
