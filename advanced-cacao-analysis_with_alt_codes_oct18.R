library(ggplot2)
library(dplyr)
library(DESeq2)



res <- results(dds, contrast = list('time_point29.partRoot.treatmentDrought', 
                                    'time_point29.partRoot.treatmentWatered'), alpha = 0.05)

# Convert to tibble and filter
df.res <- as_tibble(res, rownames = 'gene_id')
df.top <- df.res %>% 
  filter(baseMean > 1) %>% 
  mutate(diffexpressed = case_when(log2FoldChange > 1.5 & padj < 0.05 ~ 'UP',
                                   log2FoldChange < -1.5 & padj < 0.05 ~ 'DOWN',
                                   TRUE ~ 'NO'))

# Add labels for UP and DOWN genes
df.top <- df.top %>%
  mutate(lab = if_else(diffexpressed %in% c('UP', 'DOWN'), gene_id, NA))

# Volcano plot
ggplot(df.top, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpressed)) +
  geom_point() +
  geom_text(aes(label = lab), size = 2, vjust = -1) +
  theme_minimal() +
  scale_color_manual(values = c('pink', 'black', 'turquoise')) +
  geom_vline(xintercept = c(-1.5, 1.5), colour = 'black') +  # Adjust intercept to match the cutoff (1.5)
  geom_hline(yintercept = -log10(0.05), colour = 'black') +
  theme(text = element_text(size = 20))

# Install VennDiagram and UpSetR packages
# Load required library
library(VennDiagram)

# Load required library
library(VennDiagram)

# Example gene sets (replace with actual gene lists for each condition)
drought_5am_down <- paste0("Gene", 1:676)
drought_9am_down <- paste0("Gene", 1:25)
drought_24hpt_down <- paste0("Gene", 1:1)  # 1 down-regulated gene at 24 hrs P.T.

# Create a Venn diagram for DOWN-regulated genes at different times
venn.plot_down <- venn.diagram(
  x = list(
    `5 AM` = drought_5am_down,
    `9 AM` = drought_9am_down,
    `24h PT` = drought_24hpt_down
  ),
  filename = NULL,
  category.names = c("5 AM", "9 AM", "24h PT"),
  fill = c("turquoise", "pink", "lightgreen"),
  alpha = 0.5,
  cat.cex = 1.0,
  cat.pos = c(-20, 20, 180),
  
)

# Display the Venn diagram
grid.draw(venn.plot_down)

# Load required library
library(VennDiagram)

# Example gene sets for DOWN-regulated genes at different times (replace with actual lists)
drought_5am_down <- paste0("Gene", 1:651)
drought_9am_down <- paste0("Gene", 1:24)
drought_24hpt_down <- "Gene1"

# Create a Venn diagram for DOWN-regulated genes with better aesthetics
venn.plot_down <- venn.diagram(
  x = list(
    `5 AM` = drought_5am_down,
    `9 AM` = drought_9am_down,
    `24 h PT` = drought_24hpt_down
  ),
  filename = NULL,  # We want to display it in R, not save it as a file
  category.names = c("5 AM", "9 AM", "24 h PT"),
  fill = c("turquoise", "lightpink", "lightgreen"),  # Softer color palette
  alpha = 0.7,  # Set transparency
  cex = 1.0,  # Larger font size for numbers
  cat.cex = 1.0,  # Larger font size for category labels
  cat.pos = c(-20, 20, 180),  # Adjust positions to avoid overlap
  cat.dist = c(0.06, 0.06, 0.1),  # Adjust distances of labels from the circles
  fontfamily = "sans",  # Clean sans-serif font
  lwd = 2,  # Thicker lines for circles
  
  

# Display the Venn diagram
grid.draw(venn.plot_down)

# Load necessary library
library(UpSetR)

# Create a data frame with correct sizes
upset_data <- data.frame(
  `5 AM (Up)` = c(rep(1, 598), rep(0, 0), rep(0, 676), rep(0, 0), rep(0, 0), rep(0, 0), rep(0, 0)),
  `5 AM (Down)` = c(rep(0, 598), rep(1, 676), rep(0, 0), rep(0, 0), rep(0, 0), rep(0, 0), rep(0, 0)),
  `9 AM (Up)` = c(rep(0, 598), rep(0, 676), rep(1, 2), rep(0, 0), rep(0, 0), rep(0, 0), rep(0, 0)),
  `9 AM (Down)` = c(rep(0, 598), rep(0, 676), rep(0, 2), rep(1, 25), rep(0, 0), rep(0, 0), rep(0, 0)),
  `24h PT (Down)` = c(rep(0, 598), rep(0, 676), rep(0, 2), rep(0, 25), rep(1, 1))
)

# Load necessary library
library(UpSetR)

# Total number of genes
total_genes <- 1300  # Update this if your total gene count is different

# Create a data frame with the correct number of genes
upset_data <- data.frame(
  `5 AM (Up)` = c(rep(1, 598), rep(0, total_genes - 598)),  # 598 up-regulated, rest non-up
  `5 AM (Down)` = c(rep(0, 598), rep(1, 676), rep(0, total_genes - 598 - 676)),  # 676 down-regulated, rest non-down
  `9 AM (Up)` = c(rep(0, total_genes - 2), rep(1, 2)),  # 2 up-regulated at 9 AM
  `9 AM (Down)` = c(rep(0, total_genes - 25), rep(1, 25)),  # 25 down-regulated at 9 AM
  `24h PT (Down)` = c(rep(0, total_genes - 1), 1)  # 1 down-regulated at 24h PT
)

# Check the data frame structure
str(upset_data)

# Create UpSet plot with modified sets argument
upset(upset_data, sets = c("X5.AM..Up.", "X5.AM..Down.", "X9.AM..Up.", "X9.AM..Down.", "X24h.PT..Down."),
      main.bar.color = "turquoise",
      sets.bar.color = "pink",
      order.by = "freq",
      mainbar.y.label = "Gene Overlap",
      sets.x.label = "Gene Set Size")

# Load required library
library(VennDiagram)

# Example gene sets (replace with actual gene lists for each condition)
drought_5am_up <- paste0("Gene", 1:598)
drought_5am_down <- paste0("Gene", 1:676)

drought_9am_up <- paste0("Gene", 1:2)
drought_9am_down <- paste0("Gene", 1:25)

drought_24hpt_down <- paste0("Gene", 1:1) # 1 down-regulated gene

# Create a Venn diagram for UP-regulated genes at different times
venn.plot <- venn.diagram(
  x = list(
    `5 AM` = drought_5am_up,
    `9 AM` = drought_9am_up
  ),
  filename = NULL,
  category.names = c("5 AM", "9 AM"),
  fill = c("turquoise", "pink"),
  alpha = 0.5,
  cat.cex = 0.6,
  cat.pos = 0.6,
  cat.dist = 0.06
)

# Display the Venn diagram
grid.draw(venn.plot)


# Load necessary library
library(UpSetR)

# Total number of genes
total_genes <- 2000  # Adjust based on the sum of UP and DOWN gene counts provided for all time points

# Create a data frame with the correct number of genes for UP and DOWN categories at different time points
upset_data <- data.frame(
  `5 AM (Up)` = c(rep(1, 1308), rep(0, total_genes - 1308)),  # 1308 up-regulated at 5 AM
  `9 AM (Up)` = c(rep(0, 1308), rep(1, 3), rep(0, total_genes - 1308 - 3)),  # 3 up-regulated at 9 AM
  `24h PT (Up)` = c(rep(0, 1308 + 3), rep(1, 22), rep(0, total_genes - 1308 - 3 - 22)),  # 22 up-regulated at 24h PT
  `21h PT (Up)` = c(rep(0, 1308 + 3 + 22), rep(1, 25), rep(0, total_genes - 1308 - 3 - 22 - 25)),  # 25 up-regulated at 21h PT
  `17h PT (Up)` = c(rep(0, 1308 + 3 + 22 + 25), rep(1, 35), rep(0, total_genes - 1308 - 3 - 22 - 25 - 35)),  # 35 up-regulated at 17h PT
  `13h PT (Up)` = c(rep(0, 1308 + 3 + 22 + 25 + 35), rep(1, 5), rep(0, total_genes - 1308 - 3 - 22 - 25 - 35 - 5)),  # 5 up-regulated at 13h PT
  
  `5 AM (Down)` = c(rep(1, 1048), rep(0, total_genes - 1048)),  # 1048 down-regulated at 5 AM
  `9 AM (Down)` = c(rep(0, 1048), rep(1, 16), rep(0, total_genes - 1048 - 16)),  # 16 down-regulated at 9 AM
  `24h PT (Down)` = c(rep(0, 1048 + 16), rep(1, 21), rep(0, total_genes - 1048 - 16 - 21)),  # 21 down-regulated at 24h PT
  `21h PT (Down)` = c(rep(0, 1048 + 16 + 21), rep(1, 48), rep(0, total_genes - 1048 - 16 - 21 - 48)),  # 48 down-regulated at 21h PT
  `17h PT (Down)` = c(rep(0, 1048 + 16 + 21 + 48), rep(1, 35), rep(0, total_genes - 1048 - 16 - 21 - 48 - 35)),  # 35 down-regulated at 17h PT
  `13h PT (Down)` = c(rep(0, 1048 + 16 + 21 + 48 + 35), rep(1, 12), rep(0, total_genes - 1048 - 16 - 21 - 48 - 35 - 12))  # 12 down-regulated at 13h PT
)

# Check the structure of the data frame
str(upset_data)

# Create UpSet plot with the modified sets argument and the gene count data
upset(
  upset_data, 
  sets = c("5 AM (Up)", "9 AM (Up)", "24h PT (Up)", "21h PT (Up)", "17h PT (Up)", "13h PT (Up)", 
           "5 AM (Down)", "9 AM (Down)", "24h PT (Down)", "21h PT (Down)", "17h PT (Down)", "13h PT (Down)"),
  main.bar.color = "turquoise",          # Color for the main intersection bars
  sets.bar.color = "pink",               # Color for the set size bars
  order.by = "freq",                     # Order the intersections by frequency
  mainbar.y.label = "Gene Overlap",      # Y-axis label for the main bar plot
  sets.x.label = "Gene Set Size"         # X-axis label for the set size plot
)

#
# Load necessary library
library(UpSetR)

# Calculate total number of genes
total_genes <- 1308 + 3 + 22 + 25 + 35 + 5 + 1048 + 16 + 21 + 48 + 63 + 12

# Create a data frame with the correct number of genes for UP and DOWN categories at different time points
upset_data <- data.frame(
  `5 AM (Up)` = c(rep(1, 1308), rep(0, total_genes - 1308)),
  `9 AM (Up)` = c(rep(1, 3), rep(0, total_genes - 3)),
  `24h PT (Up)` = c(rep(1, 22), rep(0, total_genes - 22)),
  `21h PT (Up)` = c(rep(1, 25), rep(0, total_genes - 25)),
  `17h PT (Up)` = c(rep(1, 35), rep(0, total_genes - 35)),
  `13h PT (Up)` = c(rep(1, 5), rep(0, total_genes - 5)),
  `5 AM (Down)` = c(rep(1, 1048), rep(0, total_genes - 1048)),
  `9 AM (Down)` = c(rep(1, 16), rep(0, total_genes - 16)),
  `24h PT (Down)` = c(rep(1, 21), rep(0, total_genes - 21)),
  `21h PT (Down)` = c(rep(1, 48), rep(0, total_genes - 48)),
  `17h PT (Down)` = c(rep(1, 63), rep(0, total_genes - 63)),
  `13h PT (Down)` = c(rep(1, 12), rep(0, total_genes - 12))
)

# Check the structure of the data frame
str(upset_data)

# Create UpSet plot
upset(
  upset_data, 
  sets = c("5 AM (Up)", "9 AM (Up)", "24h PT (Up)", "21h PT (Up)", "17h PT (Up)", "13h PT (Up)", 
           "5 AM (Down)", "9 AM (Down)", "24h PT (Down)", "21h PT (Down)", "17h PT (Down)", "13h PT (Down)"),
  main.bar.color = "turquoise",
  sets.bar.color = "pink",
  order.by = "freq",
  mainbar.y.label = "Gene Overlap",
  sets.x.label = "Gene Set Size"
)


#
> # Define counts for upregulated and downregulated genes
  > up_counts <- c(`5 AM` = 1308, `9 AM` = 3, `17 h PT` = 1, `24 h PT` = 22, `21 h PT` = 25,`13 h PT` = 5)
> down_counts <- c(`5 AM` = 1048, `9 AM` = 16, `24 h PT` = 21, `21 h PT` = 48, `17 h PT` = 63, `13 h PT` = 12)
> # Calculate the total number of samples
  > total_samples <- 1300  # Adjust this based on your total number of genes analyzed (up + down)
> # Create UpSet data
  > upset_data <- data.frame(
    +   `5 AM (Up)` = c(rep(1, up_counts["5 AM"]), rep(0, total_samples - up_counts["5 AM"])),
    +   `9 AM (Up)` = c(rep(0, up_counts["5 AM"]), rep(1, up_counts["9 AM"]), rep(0, total_samples - (up_counts["5 AM"] + up_counts["9 AM"]))),
    +   `17 h PT (Up)` = c(rep(0, up_counts["5 AM"] + up_counts["9 AM"]), rep(1, up_counts["17 h PT"]), rep(0, total_samples - (up_counts["5 AM"] + up_counts["9 AM"] + up_counts["17 h PT"]))),
    +   `5 AM (Down)` = c(rep(1, down_counts["5 AM"]), rep(0, total_samples - down_counts["5 AM"])),
    +   `9 AM (Down)` = c(rep(0, down_counts["5 AM"]), rep(1, down_counts["9 AM"]), rep(0, total_samples - (down_counts["5 AM"] + down_counts["9 AM"]))),
    +   `24 h PT (Down)` = c(rep(0, down_counts["5 AM"] + down_counts["9 AM"]), rep(1, down_counts["24 h PT"]), rep(0, total_samples - (down_counts["5 AM"] + down_counts["9 AM"] + down_counts["24 h PT"]))),
    +   `21 h PT (Down)` = c(rep(0, down_counts["5 AM"] + down_counts["9 AM"] + down_counts["24 h PT"]), rep(1, down_counts["21 h PT"]), rep(0, total_samples - (down_counts["5 AM"] + down_counts["9 AM"] + down_counts["24 h PT"] + down_counts["21 h PT"]))),
    +   `17 h PT (Down)` = c(rep(0, down_counts["5 AM"] + down_counts["9 AM"] + down_counts["24 h PT"] + down_counts["21 h PT"]), rep(1, down_counts["17 h PT"]), rep(0, total_samples - (down_counts["5 AM"] + down_counts["9 AM"] + down_counts["24 h PT"] + down_counts["21 h PT"] + down_counts["17 h PT"]))),
    +   `13 h PT (Down)` = c(rep(0, down_counts["5 AM"] + down_counts["9 AM"] + down_counts["24 h PT"] + down_counts["21 h PT"] + down_counts["17 h PT"]), rep(1, down_counts["13 h PT"]), rep(0, total_samples - (down_counts["5 AM"] + down_counts["9 AM"] + down_counts["24 h PT"] + down_counts["21 h PT"] + down_counts["17 h PT"] + down_counts["13 h PT"])))
    + )
> # Print the structure of the upset data for verification
  > str(upset_data)
'data.frame':	1300 obs. of  9 variables:
  $ X5.AM..Up.     : num  1 1 1 1 1 1 1 1 1 1 ...
$ X9.AM..Up.     : num  0 0 0 0 0 0 0 0 0 0 ...
$ X17.h.PT..Up.  : num  0 0 0 0 0 0 0 0 0 0 ...
$ X5.AM..Down.   : num  1 1 1 1 1 1 1 1 1 1 ...
$ X9.AM..Down.   : num  0 0 0 0 0 0 0 0 0 0 ...
$ X24.h.PT..Down.: num  0 0 0 0 0 0 0 0 0 0 ...
$ X21.h.PT..Down.: num  0 0 0 0 0 0 0 0 0 0 ...
$ X17.h.PT..Down.: num  0 0 0 0 0 0 0 0 0 0 ...
$ X13.h.PT..Down.: num  0 0 0 0 0 0 0 0 0 0 ...
> # Create the UpSet plot
  > upset(upset_data, sets = colnames(upset_data), order.by = "freq")


#
library(UpSetR)

# Calculate total number of genes
total_genes <- 1715 + 36 + 34 + 40 + 65 + 45 + 11 + 1210 + 14 + 38 + 14 + 31 + 28 + 4

# Create a data frame with the correct number of genes for UP and DOWN categories at different time points
upset_data <- data.frame(
  `5 AM (Up)` = c(rep(1, 1210), rep(0, total_genes - 1210)),
  `9 AM (Up)` = c(rep(1, 14), rep(0, total_genes - 14)),
  `25h PT (Up)` = c(rep(1, 28), rep(0, total_genes - 28)),
  `29h PT (Up)` = c(rep(1, 4), rep(0, total_genes - 4)),
  `21h PT (Up)` = c(rep(1, 38), rep(0, total_genes - 38)),
  `17h PT (Up)` = c(rep(1, 31), rep(0, total_genes - 31)),
  `13h PT (Up)` = c(rep(1, 14), rep(0, total_genes - 14)),
  `5 AM (Down)` = c(rep(1, 1715), rep(0, total_genes - 1715)),
  `9 AM (Down)` = c(rep(1, 40), rep(0, total_genes - 40)),
  `25h PT (Down)` = c(rep(1, 45), rep(0, total_genes - 45)),
  `29h PT (Down)` = c(rep(1, 11), rep(0, total_genes - 11)),
  `21h PT (Down)` = c(rep(1, 34), rep(0, total_genes - 34)),
  `17h PT (Down)` = c(rep(1, 65), rep(0, total_genes - 65)),
  `13h PT (Down)` = c(rep(1, 36), rep(0, total_genes - 36))
)

# Check the structure of the data frame
str(upset_data)

# Create UpSet plot
upset(
  upset_data, 
  sets = c("5 AM (Up)", "9 AM (Up)", "25h PT (Up)", "29h PT (Up)", "21h PT (Up)", "17h PT (Up)", "13h PT (Up)", 
           "5 AM (Down)", "9 AM (Down)", "25h PT (Down)", "29h PT (Up)", "21h PT (Down)", "17h PT (Down)", "13h PT (Down)"),
  main.bar.color = "turquoise",
  sets.bar.color = "pink",
  order.by = "freq",
  mainbar.y.label = "Gene Overlap",
  sets.x.label = "Gene Set Size"
)

#
# Load necessary library
library(UpSetR)

# Calculate total number of genes
total_genes <- 1308 + 3 + 22 + 25 + 35 + 5 + 1048 + 16 + 21 + 48 + 63 + 12

# Create a data frame with the correct number of genes for UP and DOWN categories at different time points
upset_data <- data.frame(
  `X5_AM_Up` = c(rep(1, 1308), rep(0, total_genes - 1308)),  # Replace spaces and parentheses with underscores
  `X9_AM_Up` = c(rep(1, 3), rep(0, total_genes - 3)),
  `X24h_PT_Up` = c(rep(1, 22), rep(0, total_genes - 22)),
  `X21h_PT_Up` = c(rep(1, 25), rep(0, total_genes - 25)),
  `X17h_PT_Up` = c(rep(1, 35), rep(0, total_genes - 35)),
  `X13h_PT_Up` = c(rep(1, 5), rep(0, total_genes - 5)),
  `X5_AM_Down` = c(rep(1, 1048), rep(0, total_genes - 1048)),
  `X9_AM_Down` = c(rep(1, 16), rep(0, total_genes - 16)),
  `X24h_PT_Down` = c(rep(1, 21), rep(0, total_genes - 21)),
  `X21h_PT_Down` = c(rep(1, 48), rep(0, total_genes - 48)),
  `X17h_PT_Down` = c(rep(1, 63), rep(0, total_genes - 63)),
  `X13h_PT_Down` = c(rep(1, 12), rep(0, total_genes - 12))
)

# Check the structure of the data frame
str(upset_data)

# Create UpSet plot
upset(
  upset_data, 
  sets = c("X5_AM_Up", "X9_AM_Up", "X24h_PT_Up", "X21h_PT_Up", "X17h_PT_Up", "X13h_PT_Up", 
           "X5_AM_Down", "X9_AM_Down", "X24h_PT_Down", "X21h_PT_Down", "X17h_PT_Down", "X13h_PT_Down"),
  main.bar.color = "turquoise",
  sets.bar.color = "pink",
  order.by = "freq",
  mainbar.y.label = "Gene Overlap",
  sets.x.label = "Gene Set Size"
)

#barplot for apices
# Install and load required packages if not already installed
# install.packages("ggplot2")
library(ggplot2)
library(tidyr)

# Create the data frame with counts for each time point
data <- data.frame(
  Time = c("5 am", "9 am", "13 hrs PT", "17 hrs PT", "21 hrs PT", "25 hrs PT", "29 hrs PT"),
  UP = c(1210, 14, 14, 31, 38, 28, 4),   # UP gene counts
  DOWN = c(1715, 40, 36, 65, 34, 45, 11) # DOWN gene counts
)

# Reshape the data to long format for ggplot
long_data <- data %>%
  pivot_longer(cols = c(UP, DOWN), names_to = "Regulation", values_to = "Count")

# Create the bar plot
ggplot(long_data, aes(x = Time, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gene Expression Counts in Apices Under Drought",
       x = "Time Points",
       y = "Gene Counts") +
  theme_minimal() +
  scale_fill_manual(values = c("UP" = "turquoise", "DOWN" = "lightpink")) +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))





