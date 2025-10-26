# First, let's combine all DESeq2 results into one dataframe
library(tidyverse)

# Function to convert DESeq results to dataframe with timepoint info
convert_res_to_df <- function(res, timepoint) {
  as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(timepoint = timepoint)
}

# Combine all results
leaf_all_results <- bind_rows(
  convert_res_to_df(res_5h_leaf, "5h"),
  convert_res_to_df(res_9h_leaf, "9h"),
  convert_res_to_df(res_13h_leaf, "13h"),
  convert_res_to_df(res_17h_leaf, "17h"),
  convert_res_to_df(res_21h_leaf, "21h"),
  convert_res_to_df(res_25h_leaf, "25h"),
  convert_res_to_df(res_29h_leaf, "29h")
) %>%
  mutate(diffexpressed = case_when(
    log2FoldChange >= 1.5 & padj < 0.05 ~ "UP",
    log2FoldChange <= -1.5 & padj < 0.05 ~ "DOWN",
    TRUE ~ "NO"
  ))

# Function to perform GO enrichment for specific timepoint and direction
perform_go_enrichment <- function(data, timepoint, direction, gene_set, go_terms) {
  # Filter genes for specific timepoint and direction
  signif_genes <- data %>%
    filter(timepoint == !!timepoint, 
           diffexpressed == !!direction) %>%
    arrange(padj, abs(log2FoldChange)) %>%
    pull(gene_id)
  
  # Perform enrichment
  enrichment_results <- enricher(gene = signif_genes, 
                                TERM2GENE = gene_set) %>%
    as_tibble() %>%
    mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  
  # Add GO term information and filter
  enrichment_results %>%
    left_join(go_terms, join_by(ID == ids)) %>%
    filter(category == 'biological process')
}

# Function to plot GO enrichment results
plot_go_enrichment <- function(enrichment_data, title) {
  enrichment_data %>%
    slice(1:20) %>%
    ggplot(showCategory = 50,
           aes(richFactor,
               fct_reorder(go_names, richFactor))) +
    geom_segment(aes(xend = 0, yend = go_names)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
      colours = c("#f7ca64", "#46bac2", "#7e62a3"),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    theme_dose(12) +
    xlab("Rich Factor") +
    labs(y = NULL) +
    ggtitle(title)
}

# Example usage for one timepoint:
# For upregulated genes at 5h
go_up_5h <- perform_go_enrichment(
  leaf_all_results, 
  timepoint = "5h", 
  direction = "UP",
  gene_set = gene_set,
  go_terms = go_terms
)

# Plot the results
plot_go_enrichment(go_up_5h, "Biological Processes - 5h Upregulated Genes")

# To analyze all timepoints and both directions:
timepoints <- unique(leaf_all_results$timepoint)
directions <- c("UP", "DOWN")

# Create list to store all results
all_enrichment_results <- list()

for(tp in timepoints) {
  for(dir in directions) {
    result_name <- paste0(tp, "_", dir)
    all_enrichment_results[[result_name]] <- perform_go_enrichment(
      leaf_all_results,
      timepoint = tp,
      direction = dir,
      gene_set = gene_set,
      go_terms = go_terms
    )
  }
}

#-------------------------------------------------------------------------------

#to visualize and save
library(patchwork)
library(writexl)

# Create plots for all timepoints and directions
enrichment_plots <- list()
enrichment_tables <- list()

for(tp in timepoints) {
  for(dir in directions) {
    result_name <- paste0(tp, "_", dir)
    
    # Create plot
    plot_title <- paste("Biological Processes -", tp, dir, "regulated Genes")
    enrichment_plots[[result_name]] <- all_enrichment_results[[result_name]] %>%
      slice(1:20) %>%
      ggplot(showCategory = 50,
             aes(richFactor,
                 fct_reorder(go_names, richFactor))) +
      geom_segment(aes(xend = 0, yend = go_names)) +
      geom_point(aes(color = p.adjust, size = Count)) +
      scale_color_gradientn(
        colours = c("#f7ca64", "#46bac2", "#7e62a3"),
        trans = "log10",
        guide = guide_colorbar(reverse = TRUE, order = 1)
      ) +
      scale_size_continuous(range = c(2, 10)) +
      theme_dose(12) +
      xlab("Rich Factor") +
      labs(y = NULL) +
      ggtitle(plot_title)
    
    # Store table with additional information
    enrichment_tables[[result_name]] <- all_enrichment_results[[result_name]] %>%
      select(ID, go_names, category, GeneRatio, BgRatio, 
             pvalue, p.adjust, qvalue, Count, richFactor, geneID)
  }
}

# Combine all plots using patchwork
# First, combine UP and DOWN for each timepoint
combined_plots <- list()
for(tp in timepoints) {
  up_plot <- enrichment_plots[[paste0(tp, "_UP")]]
  down_plot <- enrichment_plots[[paste0(tp, "_DOWN")]]
  combined_plots[[tp]] <- up_plot + down_plot +
    plot_layout(guides = "collect") +
    plot_annotation(title = paste("Timepoint:", tp))
}

# Save the combined plot with adjusted dimensions
ggsave("all_timepoints_GO_enrichment.pdf", final_plot, 
       width = 15, height = length(timepoints) * 6,  # Reduced height per timepoint
       units = "in", limitsize = FALSE)  # Added limitsize parameter

# Save individual plots for each timepoint with adjusted dimensions
for(tp in timepoints) {
  plot_name <- paste0("GO_enrichment_", tp, ".pdf")
  ggsave(plot_name, combined_plots[[tp]], 
         width = 15, height = 6,  # Adjusted height
         units = "in")
}

# Alternative approach: Split into multiple pages if still too large
# Calculate number of timepoints per page (3 timepoints per page)
timepoints_per_page <- 3
num_pages <- ceiling(length(timepoints) / timepoints_per_page)

# Create and save multiple pages
for(page in 1:num_pages) {
  # Get timepoints for this page
  start_idx <- (page - 1) * timepoints_per_page + 1
  end_idx <- min(page * timepoints_per_page, length(timepoints))
  page_timepoints <- timepoints[start_idx:end_idx]
  
  # Combine plots for this page
  page_plots <- wrap_plots(combined_plots[page_timepoints], ncol = 1)
  
  # Save this page
  ggsave(
    paste0("GO_enrichment_page_", page, ".pdf"),
    page_plots,
    width = 15,
    height = length(page_timepoints) * 6,
    units = "in"
  )
}



# Save tables to Excel - one sheet per timepoint and direction
# Create a list of sheets for the Excel file
excel_sheets <- list()
for(tp in timepoints) {
  for(dir in directions) {
    sheet_name <- paste0(tp, "_", dir)
    excel_sheets[[sheet_name]] <- enrichment_tables[[sheet_name]]
  }
}

# Save to Excel file
write_xlsx(excel_sheets, "GO_enrichment_results.xlsx")

# Also create a summary table
summary_table <- tibble(
  Timepoint = character(),
  Direction = character(),
  Total_DEGs = numeric(),
  Significant_GO_Terms = numeric(),
  Top_GO_Terms = character()
)

for(tp in timepoints) {
  for(dir in directions) {
    result_name <- paste0(tp, "_", dir)
    
    # Get top 5 GO terms
    top_terms <- enrichment_tables[[result_name]] %>%
      slice(1:5) %>%
      pull(go_names) %>%
      paste(collapse = "; ")
    
    # Add to summary table
    summary_table <- summary_table %>%
      add_row(
        Timepoint = tp,
        Direction = dir,
        Total_DEGs = length(unique(unlist(strsplit(enrichment_tables[[result_name]]$geneID, "/")))),
        Significant_GO_Terms = nrow(enrichment_tables[[result_name]]),
        Top_GO_Terms = top_terms
      )
  }
}

# Add summary table to Excel file
write_xlsx(list(
  Summary = summary_table,
  `All Results` = bind_rows(enrichment_tables, .id = "Timepoint_Direction")
), "GO_enrichment_summary.xlsx")

# Print summary to console
print(summary_table)

# Return the plot object for display
final_plot



#-------------------------------------------------------------------------------
# Get all significant genes across all time points
leaf_all_signif <- leaf_all_results %>%
  filter(diffexpressed %in% c("UP", "DOWN")) %>%
  distinct(gene_id) %>%
  pull(gene_id)

# Split by up/down regulation across all time points
leaf_up_all <- leaf_all_results %>%
  filter(diffexpressed == "UP") %>%
  distinct(gene_id) %>%
  pull(gene_id)

leaf_down_all <- leaf_all_results %>%
  filter(diffexpressed == "DOWN") %>%
  distinct(gene_id) %>%
  pull(gene_id)

# Perform enrichment for all combined genes
leaf_all_enrichment <- enricher(gene = leaf_all_signif, 
                                TERM2GENE = gene_set) %>%
  as_tibble() %>%
  mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

# Perform separate enrichment for up and down regulated genes
leaf_up_enrichment <- enricher(gene = leaf_up_all, 
                               TERM2GENE = gene_set) %>%
  as_tibble() %>%
  mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

leaf_down_enrichment <- enricher(gene = leaf_down_all, 
                                 TERM2GENE = gene_set) %>%
  as_tibble() %>%
  mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

# Create plots for each analysis
# 1. All significant genes plot
plot_all <- leaf_all_enrichment %>%
  left_join(go_terms, join_by(ID == ids)) %>%
  filter(category == 'biological process') %>%
  slice(1:20) %>%
  ggplot(showCategory = 50,
         aes(richFactor,
             fct_reorder(go_names, richFactor))) +
  geom_segment(aes(xend = 0, yend = go_names)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(
    colours = c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10",
    guide = guide_colorbar(reverse = TRUE, order = 1)
  ) +
  scale_size_continuous(range = c(2, 10)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  labs(y = NULL) +
  ggtitle("Biological Processes - All DEGs Combined")

# 2. Upregulated genes plot
plot_up <- leaf_up_enrichment %>%
  left_join(go_terms, join_by(ID == ids)) %>%
  filter(category == 'biological process') %>%
  slice(1:20) %>%
  ggplot(showCategory = 50,
         aes(richFactor,
             fct_reorder(go_names, richFactor))) +
  geom_segment(aes(xend = 0, yend = go_names)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(
    colours = c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10",
    guide = guide_colorbar(reverse = TRUE, order = 1)
  ) +
  scale_size_continuous(range = c(2, 10)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  labs(y = NULL) +
  ggtitle("Biological Processes - Upregulated DEGs")

# 3. Downregulated genes plot
plot_down <- leaf_down_enrichment %>%
  left_join(go_terms, join_by(ID == ids)) %>%
  filter(category == 'biological process') %>%
  slice(1:20) %>%
  ggplot(showCategory = 50,
         aes(richFactor,
             fct_reorder(go_names, richFactor))) +
  geom_segment(aes(xend = 0, yend = go_names)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(
    colours = c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10",
    guide = guide_colorbar(reverse = TRUE, order = 1)
  ) +
  scale_size_continuous(range = c(2, 10)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  labs(y = NULL) +
  ggtitle("Biological Processes - Downregulated DEGs")

# Arrange all plots together using patchwork
library(patchwork)
combined_plots <- plot_all / plot_up / plot_down +
  plot_layout(heights = c(1, 1, 1))

# Display the combined plot
combined_plots

# Optional: Save the plots
ggsave("leaf_combined_GO_enrichment.pdf", combined_plots, 
       width = 12, height = 20, units = "in")

# Get summary statistics
summary_stats <- tibble(
  Category = c("All DEGs", "Upregulated", "Downregulated"),
  Gene_Count = c(length(leaf_all_signif), 
                 length(leaf_up_all), 
                 length(leaf_down_all)),
  Enriched_GO_Terms = c(nrow(leaf_all_enrichment),
                        nrow(leaf_up_enrichment),
                        nrow(leaf_down_enrichment))
)

print(summary_stats)
