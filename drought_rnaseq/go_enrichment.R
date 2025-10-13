# run_go_enrichment.R
library(clusterProfiler)
library(org.At.tair.db)  # Replace with cacao annotation
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

deg <- read_csv(input_file)

sig_genes <- deg %>%
  filter(padj < 0.05) %>%
  pull(Gene_ID)

ego <- enrichGO(
  gene = sig_genes,
  OrgDb = org.At.tair.db,
  keyType = \"SYMBOL\",
  ont = \"BP\",
  pAdjustMethod = \"BH\",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.csv(as.data.frame(ego), output_file, row.names = FALSE)
