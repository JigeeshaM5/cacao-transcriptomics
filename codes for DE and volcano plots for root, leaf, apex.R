#code for the DE and volcano plot Feb 22, 2025
library(PCAtools)
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(UpSetR)  # for upset plots
library(VennDiagram)  # for venn diagrams
library(dplyr)
library(grid)
library(gridExtra)
library(tidyverse)
library(WGCNA)
library(clusterProfiler)
library(EnhancedVolcano)
library(fgsea)
library(ComplexHeatmap)
library(circlize)

#data path
expr_data_path <- "C:/Users/Jigeesha Mukherjee/Desktop/Jigeesha/work 2023/job applications/guiltinan/For R/counts.tsv"
#-------------------------------------------------------------------------------

df <- read.table("counts.tsv", header = TRUE, sep = "\t")
data <- read.delim(expr_data_path, row.names = 1, check.names = FALSE)
count_data <- round(data)

#-------------------------------------------------------------------------------

metadata <- data.frame(
  sample = colnames(count_data),
  treatment = gsub("DD-(.)T.-T.*", "\\1", colnames(count_data)),
  tissue = gsub("DD-.T(.)-T.*", "\\1", colnames(count_data)),
  time = gsub(".*-T(\\d)_.*", "\\1", colnames(count_data)),
  replicate = gsub(".*_R(\\d)", "\\1", colnames(count_data))
)
metadata$treatment<-factor(metadata$treatment,levels = c('W','D'))

time_mapping <- c("1" = 5, "2" = 9, "3" = 13, "4" = 17, 
                  "5" = 21, "6" = 25, "7" = 29)
metadata$hours <- as.numeric(time_mapping[metadata$time])
metadata$hours<-factor(metadata$hours)

root_metadata<- metadata %>% filter(tissue =="R")
leaf_metadata<- metadata %>% filter(tissue =="L")
apex_metadata<- metadata %>% filter(tissue =="A")

#for root
# Create a new factor combining treatment and time
root_metadata$group <- factor(paste(root_metadata$treatment, root_metadata$hours, sep="_"))
# Create new DESeq object with modified design
dds_root_new <- DESeqDataSetFromMatrix(
  countData = root_count_data,
  colData = root_metadata,
  design = ~ group
)

# Run DESeq
dds_root_new <- DESeq(dds_root_new)

# Then for each time point:
res_5h_new <- results(dds_root_new, contrast=c("group", "D_5", "W_5"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_9h_new <- results(dds_root_new, contrast=c("group", "D_9", "W_9"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_13h_new <- results(dds_root_new, contrast=c("group", "D_13", "W_13"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_17h_new <- results(dds_root_new, contrast=c("group", "D_17", "W_17"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_21h_new <- results(dds_root_new, contrast=c("group", "D_21", "W_21"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_25h_new <- results(dds_root_new, contrast=c("group", "D_25", "W_25"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_29h_new <- results(dds_root_new, contrast=c("group", "D_29", "W_29"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

#-------------------------------------------------------------------------------
# Combine top 100 genes for each time point
all_top_genes <- rbind(
  cbind(top_5h, Timepoint = "5h"),
  cbind(top_9h, Timepoint = "9h"),
  cbind(top_13h, Timepoint = "13h"),
  cbind(top_17h, Timepoint = "17h"),
  cbind(top_21h, Timepoint = "21h"),
  cbind(top_25h, Timepoint = "25h"),
  cbind(top_29h, Timepoint = "29h")
)

# Separate into upregulated (lfc > 0) and downregulated (lfc < 0)
upregulated_genes <- all_top_genes[all_top_genes$log2FoldChange > 0, ]
downregulated_genes <- all_top_genes[all_top_genes$log2FoldChange < 0, ]

# Count occurrences of upregulated genes across time points
up_gene_counts <- table(upregulated_genes$Gene_ID)

# Keep genes that appear in at least 2 time points
overlapping_up_genes <- names(up_gene_counts[up_gene_counts >= 2])

# Subset and sort upregulated genes
overlap_up_results <- upregulated_genes[upregulated_genes$Gene_ID %in% overlapping_up_genes, ]
overlap_up_results <- overlap_up_results[order(overlap_up_results$padj, -overlap_up_results$baseMean), ]

# Select top 50 upregulated genes
top_50_up <- head(overlap_up_results, 50)

# Save to Excel
write.xlsx(top_50_up, file = "DESeq2_Top50_Overlapping_Upregulated.xlsx", row.names = FALSE)
#-------------------------------------------------------------------------------
# Count occurrences of downregulated genes across time points
down_gene_counts <- table(downregulated_genes$Gene_ID)

# Keep genes that appear in at least 2 time points
overlapping_down_genes <- names(down_gene_counts[down_gene_counts >= 2])

# Subset and sort downregulated genes
overlap_down_results <- downregulated_genes[downregulated_genes$Gene_ID %in% overlapping_down_genes, ]
overlap_down_results <- overlap_down_results[order(overlap_down_results$padj, -overlap_down_results$baseMean), ]

# Select top 50 downregulated genes
top_50_down <- head(overlap_down_results, 50)

# Save to Excel
write.xlsx(top_50_down, file = "DESeq2_Top50_Overlapping_Downregulated.xlsx", row.names = FALSE)
#-------------------------------------------------------------------------------


# Add gene names as a column to each result dataframe
res_5h_new$Gene_ID <- rownames(res_5h_new)
res_9h_new$Gene_ID <- rownames(res_9h_new)
res_13h_new$Gene_ID <- rownames(res_13h_new)
res_17h_new$Gene_ID <- rownames(res_17h_new)
res_21h_new$Gene_ID <- rownames(res_21h_new)
res_25h_new$Gene_ID <- rownames(res_25h_new)
res_29h_new$Gene_ID <- rownames(res_29h_new)

# Function to get top 100 significant genes sorted by padj, then arrange by baseMean
get_top_100 <- function(df) {
  df <- df[order(df$padj), ] # Sort by padj
  top_100 <- head(df, 100)  # Select top 100
  top_100 <- top_100[order(-top_100$baseMean), ] # Rearrange by baseMean in descending order
  return(top_100)
}

# Apply the function to each result dataframe
top_5h <- get_top_100(res_5h_new)
top_9h <- get_top_100(res_9h_new)
top_13h <- get_top_100(res_13h_new)
top_17h <- get_top_100(res_17h_new)
top_21h <- get_top_100(res_21h_new)
top_25h <- get_top_100(res_25h_new)
top_29h <- get_top_100(res_29h_new)

# Save all to an Excel file
library(openxlsx)

write.xlsx(list(
  "5h" = top_5h,
  "9h" = top_9h,
  "13h" = top_13h,
  "17h" = top_17h,
  "21h" = top_21h,
  "25h" = top_25h,
  "29h" = top_29h
), file = "DESeq2_Top100_Genes.xlsx", row.names = FALSE)

# Combine all time points into one data frame and include a "Timepoint" column
combined_results <- rbind(
  cbind(top_5h, Timepoint = "5h"),
  cbind(top_9h, Timepoint = "9h"),
  cbind(top_13h, Timepoint = "13h"),
  cbind(top_17h, Timepoint = "17h"),
  cbind(top_21h, Timepoint = "21h"),
  cbind(top_25h, Timepoint = "25h"),
  cbind(top_29h, Timepoint = "29h")
)

# Sort by padj (ascending) and then by baseMean (descending)
combined_results <- combined_results[order(combined_results$padj, -combined_results$baseMean), ]

# Select the top 100 genes across all timepoints
top_100_combined <- head(combined_results, 100)

# Save to an Excel file
write.xlsx(top_100_combined, file = "DESeq2_Top100_Combined.xlsx", row.names = FALSE)

# Optional: Output the data to check the structure in the console
head(top_100_combined)


#-------------------------------------------------------------------------------
#volcano plot function
create_volcano_plot <- function(res_obj, title) {
  # Convert to data frame and add diffexpressed column
  df <- as.data.frame(res_obj) %>%
    mutate(
      diffexpressed = case_when(
        padj < 0.05 & log2FoldChange > 1.5 & baseMean > 30 ~ "UP",
        padj < 0.05 & log2FoldChange < -1.5 & baseMean > 30 ~ "DOWN",
        TRUE ~ "NO"
      )
    )
  
  # Create plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = diffexpressed), alpha = 0.5) +
    scale_color_manual(values = c("DOWN" = "pink", "NO" = "black", "UP" = "turquoise")) +
    geom_vline(xintercept = c(-1.5, 1.5),linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
    theme_minimal() +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    )
  
  # Print summary
  summary <- df %>%
    group_by(diffexpressed) %>%
    summarise(count = n()) %>%
    filter(diffexpressed %in% c("UP", "DOWN"))
  print(summary)
  
  return(p)
}

# Use for each time point
p5 <- create_volcano_plot(res_5h_new, "Root - 5 hours")
p9 <- create_volcano_plot(res_9h_new, "Root - 9 hours")
p13 <- create_volcano_plot(res_13h_new, "Root - 13 hours")
p17 <- create_volcano_plot(res_17h_new, "Root - 17 hours")
p21 <- create_volcano_plot(res_21h_new, "Root - 21 hours")
p25 <- create_volcano_plot(res_25h_new, "Root - 25 hours")
p29 <- create_volcano_plot(res_29h_new, "Root - 29 hours")

# If you want to display all plots together, you can use patchwork
library(patchwork)
(p5 + p9 + p13) / (p17 + p21 + p25) / (p29 + plot_spacer() + plot_spacer())

#for root Feb 22, 2025
# Use for each time point
p5 <- create_volcano_plot(res_5h_new, "Root - 5 hours")
p9 <- create_volcano_plot(res_9h_new, "Root - 9 hours")
p13 <- create_volcano_plot(res_13h_new, "Root - 13 hours")
p17 <- create_volcano_plot(res_17h_new, "Root - 17 hours")
p21 <- create_volcano_plot(res_21h_new, "Root - 21 hours")
p25 <- create_volcano_plot(res_25h_new, "Root - 25 hours")
p29 <- create_volcano_plot(res_29h_new, "Root - 29 hours")

# Combine plots using patchwork
library(patchwork)
combined_plot <- (p5 + p9) / 
  (p13 + p17) / 
  (p21 + p25) / 
  (plot_spacer() + p29 + plot_spacer())  # Center p29 in the final row

# Add a title to the combined plot
final_plot <- combined_plot +
  plot_annotation(
    title = "Volcano plots corresponding to seven time points for root DEGs",
    theme = theme(
      plot.title = element_text(
        size = 16,                  # Increase title font size
        face = "bold",              # Make title bold
        hjust = 0.5,                # Center the title horizontally
        vjust = 1,                  # Adjust vertical position
        margin = margin(t = 10, b = 20)  # Add margin above and below the title
      ),
      plot.title.position = "plot"  # Place the title at the bottom
    )
  )

# Display the final plot
final_plot
#-------------------------------------------------------------------------------


#for leaf
dds_leaf <- DESeqDataSetFromMatrix(
  countData = leaf_count_data,
  colData = leaf_metadata,
  design = ~ hours + treatment + hours:treatment 
)
dds_leaf <- DESeq(dds_leaf)

resultsNames(dds_leaf)

leaf_metadata$group <- factor(paste(leaf_metadata$treatment, leaf_metadata$hours, sep="_"))

# Create new DESeq object with modified design
dds_leaf_new <- DESeqDataSetFromMatrix(
  countData = leaf_count_data,
  colData = leaf_metadata,
  design = ~ group
)

# Run DESeq
dds_leaf_new <- DESeq(dds_leaf_new)

# Then for each time point:
res_5h_leaf <- results(dds_leaf_new, contrast=c("group", "D_5", "W_5"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_9h_leaf <- results(dds_leaf_new, contrast=c("group", "D_9", "W_9"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_13h_leaf <- results(dds_leaf_new, contrast=c("group", "D_13", "W_13"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_17h_leaf <- results(dds_leaf_new, contrast=c("group", "D_17", "W_17"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_21h_leaf <- results(dds_leaf_new, contrast=c("group", "D_21", "W_21"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_25h_leaf <- results(dds_leaf_new, contrast=c("group", "D_25", "W_25"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_29h_leaf <- results(dds_leaf_new, contrast=c("group", "D_29", "W_29"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

#example use
p5_leaf <- create_volcano_plot(res_5h_leaf, "Volcano Plot - 5 hours")
p9_leaf <- create_volcano_plot(res_9h_leaf, "Volcano Plot - 9 hours")
p17_leaf <- create_volcano_plot(res_17h_leaf, "Volcano Plot - 17 hours")
p21_leaf <- create_volcano_plot(res_21h_leaf, "Volcano Plot - 21 hours")
p25_leaf <- create_volcano_plot(res_25h_leaf, "Volcano Plot - 25 hours")
p29_leaf <- create_volcano_plot(res_29h_leaf, "Volcano Plot - 29 hours")
p13_leaf <- create_volcano_plot(res_13h_leaf, "Volcano Plot - 13 hours")
#all together
library(patchwork)
(p5_leaf + p9_leaf + p13_leaf) / (p17_leaf + p21_leaf + p25_leaf) / (p29_leaf + plot_spacer() + plot_spacer())

# Combine plots using patchwork
library(patchwork)
combined_plot_leaf <- (p5_leaf + p9_leaf) / 
  (p13_leaf + p17_leaf) / 
  (p21_leaf + p25_leaf) / 
  (plot_spacer() + p29_leaf + plot_spacer())  # Center p29 in the final row

# Add a title to the combined plot
final_plot_leaf <- combined_plot_leaf +
  plot_annotation(
    title = "Volcano plots corresponding to seven time points for leaf DEGs",
    theme = theme(
      plot.title = element_text(
        size = 16,                  # Increase title font size
        face = "bold",              # Make title bold
        hjust = 0.5,                # Center the title horizontally
        vjust = 1,                  # Adjust vertical position
        margin = margin(t = 10, b = 20)  # Add margin above and below the title
      ),
      plot.title.position = "plot"  # Place the title at the bottom
    )
  )

# Display the final plot
final_plot_leaf

#-------------------------------------------------------------------------------
#for apex
dds_apex <- DESeqDataSetFromMatrix(
  countData = apex_count_data,
  colData = apex_metadata,
  design = ~ hours + treatment + hours:treatment 
)
dds_apex <- DESeq(dds_apex)

resultsNames(dds_apex)

apex_metadata$group <- factor(paste(apex_metadata$treatment, apex_metadata$hours, sep="_"))

# Create new DESeq object with modified design
dds_apex_new <- DESeqDataSetFromMatrix(
  countData = apex_count_data,
  colData = apex_metadata,
  design = ~ group
)

#jan 27
# Create a new factor combining treatment and time
apex_metadata$group <- factor(paste(apex_metadata$treatment, apex_metadata$hours, sep="_"))
# Run DESeq
dds_apex_new <- DESeq(dds_apex_new)

# Then for each time point:
res_5h_apex <- results(dds_apex_new, contrast=c("group", "D_5", "W_5"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_9h_apex <- results(dds_apex_new, contrast=c("group", "D_9", "W_9"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_13h_apex <- results(dds_apex_new, contrast=c("group", "D_13", "W_13"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_17h_apex <- results(dds_apex_new, contrast=c("group", "D_17", "W_17"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_21h_apex <- results(dds_apex_new, contrast=c("group", "D_21", "W_21"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_25h_apex <- results(dds_apex_new, contrast=c("group", "D_25", "W_25"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)

res_29h_apex <- results(dds_apex_new, contrast=c("group", "D_29", "W_29"), alpha=0.05) %>%
  subset(padj < 0.05 & abs(log2FoldChange) > 1.5 & baseMean > 30)


# Combine top 100 genes for each time point
all_top_genes <- rbind(
  cbind(top_5h, Timepoint = "5h"),
  cbind(top_9h, Timepoint = "9h"),
  cbind(top_13h, Timepoint = "13h"),
  cbind(top_17h, Timepoint = "17h"),
  cbind(top_21h, Timepoint = "21h"),
  cbind(top_25h, Timepoint = "25h"),
  cbind(top_29h, Timepoint = "29h")
)

# Add gene names as a column to each result dataframe
res_5h_apex$gene_id <- rownames(res_5h_apex)
res_9h_apex$gene_id <- rownames(res_9h_apex)
res_13h_apex$gene_id <- rownames(res_13h_apex)
res_17h_apex$gene_id <- rownames(res_17h_apex)
res_21h_apex$gene_id <- rownames(res_21h_apex)
res_25h_apex$gene_id <- rownames(res_25h_apex)
res_29h_apex$gene_id <- rownames(res_29h_apex)



# Apply the function to each result dataframe for apex
top_5h_apex <- get_top_100(res_5h_apex)
top_9h_apex <- get_top_100(res_9h_apex)
top_13h_apex <- get_top_100(res_13h_apex)
top_17h_apex <- get_top_100(res_17h_apex)
top_21h_apex <- get_top_100(res_21h_apex)
top_25h_apex <- get_top_100(res_25h_apex)
top_29h_apex <- get_top_100(res_29h_apex)

library(openxlsx)

write.xlsx(list(
  "5h" = top_5h_apex,
  "9h" = top_9h_apex,
  "13h" = top_13h_apex,
  "17h" = top_17h_apex,
  "21h" = top_21h_apex,
  "25h" = top_25h_apex,
  "29h" = top_29h_apex
), file = "DESeq2_Top100_Genes_apex.xlsx")

# Combine all time points into one data frame and include a "Timepoint" column
combined_results_apex <- rbind(
  cbind(top_5h_apex, Timepoint = "5h"),
  cbind(top_9h_apex, Timepoint = "9h"),
  cbind(top_13h_apex, Timepoint = "13h"),
  cbind(top_17h_apex, Timepoint = "17h"),
  cbind(top_21h_apex, Timepoint = "21h"),
  cbind(top_25h_apex, Timepoint = "25h"),
  cbind(top_29h_apex, Timepoint = "29h")
)

# Sort by padj (ascending) and then by baseMean (descending) for apex
combined_results_apex_sorted <- combined_results_apex[order(combined_results_apex$padj, -combined_results_apex$baseMean), ]
# Select the top 100 genes across all timepoints for apex
top_100_combined_apex <- head(combined_results_apex_sorted, 100)

# Save to an Excel file
write.xlsx(top_100_combined_apex, file = "DESeq2_Top100_Combined_apex.xlsx", row.names = FALSE)


# Combine top 100 genes for each time point for apex
all_top_genes_apex <- rbind(
  cbind(top_5h_apex, Timepoint = "5h"),
  cbind(top_9h_apex, Timepoint = "9h"),
  cbind(top_13h_apex, Timepoint = "13h"),
  cbind(top_17h_apex, Timepoint = "17h"),
  cbind(top_21h_apex, Timepoint = "21h"),
  cbind(top_25h_apex, Timepoint = "25h"),
  cbind(top_29h_apex, Timepoint = "29h")
)

all_top_genes_apex_tibble <- as_tibble(all_top_genes_apex, rownames = "gene_id")

# Separate into upregulated (lfc > 0) and downregulated (lfc < 0)
upregulated_genes <- all_top_genes[all_top_genes$log2FoldChange > 0, ]
downregulated_genes <- all_top_genes[all_top_genes$log2FoldChange < 0, ]

# Separate into upregulated (lfc > 0) and downregulated (lfc < 0) for apex
upregulated_genes_apex <- all_top_genes_apex_tibble[all_top_genes_apex_tibble$log2FoldChange > 0, ]
downregulated_genes_apex <- all_top_genes_apex_tibble[all_top_genes_apex_tibble$log2FoldChange < 0, ]
#------------------

# Count occurrences of upregulated genes across time points
library(readxl)

up_gene_counts_apex <- read_excel("Upregulated_apex.xlsx")
gene_count_table_apex <- table(up_gene_counts_apex$gene_id)



# Keep genes that appear in at least 2 time points
overlapping_up_genes_apex <- names(gene_count_table_apex[gene_count_table_apex >= 2])

# Keep genes that appear in at least 2 time points for apex
overlapping_up_genes_apex <- names(up_gene_counts_apex[up_gene_counts_apex >= 2])

# Subset and sort upregulated genes
overlap_up_results <- upregulated_genes[upregulated_genes$Gene_ID %in% overlapping_up_genes, ]
overlap_up_results <- overlap_up_results[order(overlap_up_results$padj, -overlap_up_results$baseMean), ]

# Subset and sort upregulated genes for apex
overlap_up_results_apex <- upregulated_genes_apex[upregulated_genes_apex$gene_id %in% overlapping_up_genes_apex, ]
overlap_up_results_apex <- overlap_up_results_apex[order(overlap_up_results_apex$padj, -overlap_up_results_apex$baseMean), ]

# Select top 50 upregulated genes
top_50_up <- head(overlap_up_results, 50)

# Select top 50 upregulated genes for apex
top_50_up_apex <- head(overlap_up_results_apex, 50)

# Save to Excel
write.xlsx(top_50_up, file = "DESeq2_Top50_Overlapping_Upregulated.xlsx", row.names = FALSE)

# Save to Excel for apex
write.xlsx(top_50_up_apex, file = "DESeq2_Top50_Overlapping_Upregulated_apex.xlsx", row.names = FALSE)
write.xlsx(upregulated_genes_apex, file = "Upregulated_apex.xlsx")
write.xlsx(downregulated_genes_apex, file = "Downregulated_apex.xlsx")

# Count occurrences of downregulated genes across time points
down_gene_counts_apex <- table(downregulated_genes_apex$gene_id)

# Keep genes that appear in at least 2 time points
overlapping_down_genes_apex <- names(down_gene_counts_apex[down_gene_counts_apex >= 2])

# Subset and sort downregulated genes
overlap_down_results_apex <- downregulated_genes[downregulated_genes$gene_id %in% overlapping_down_genes_apex, ]
overlap_down_results <- overlap_down_results[order(overlap_down_results$padj, -overlap_down_results$baseMean), ]

# Select top 50 downregulated genes
top_50_down <- head(overlap_down_results, 50)

# Save to Excel
write.xlsx(top_50_down, file = "DESeq2_Top50_Overlapping_Downregulated.xlsx", row.names = FALSE)

#example use#example useTRUE
p5_apex <- create_volcano_plot(res_5h_apex, "Volcano Plot - 5 hours")
p9_apex <- create_volcano_plot(res_9h_apex, "Volcano Plot - 9 hours")
p17_apex <- create_volcano_plot(res_17h_apex, "Volcano Plot - 17 hours")
p21_apex <- create_volcano_plot(res_21h_apex, "Volcano Plot - 21 hours")
p25_apex <- create_volcano_plot(res_25h_apex, "Volcano Plot - 25 hours")
p29_apex <- create_volcano_plot(res_29h_apex, "Volcano Plot - 29 hours")
p13_apex <- create_volcano_plot(res_13h_apex, "Volcano Plot - 13 hours")
#all together
library(patchwork)
(p5_apex + p9_apex + p13_apex) / (p17_apex + p21_apex + p25_apex) / (p29_apex + plot_spacer() + plot_spacer())

# Combine plots using patchwork
library(patchwork)
combined_plot_apex <- (p5_apex + p9_apex) / 
  (p13_apex + p17_apex) / 
  (p21_apex + p25_apex) / 
  (plot_spacer() + p29_apex + plot_spacer())  # Center p29 in the final row

# Add a title to the combined plot
final_plot_apex <- combined_plot_apex +
  plot_annotation(
    title = "Volcano plots corresponding to seven time points for apex DEGs",
    theme = theme(
      plot.title = element_text(
        size = 16,                  # Increase title font size
        face = "bold",              # Make title bold
        hjust = 0.5,                # Center the title horizontally
        vjust = 1,                  # Adjust vertical position
        margin = margin(t = 10, b = 20)  # Add margin above and below the title
      ),
      plot.title.position = "plot"  # Place the title at the bottom
    )
  )

# Display the final plot
final_plot_apex
