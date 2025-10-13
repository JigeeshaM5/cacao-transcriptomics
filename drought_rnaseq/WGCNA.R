# run_wgcna.R
library(WGCNA)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

expr_data <- read_csv(input_file)
datExpr <- as.data.frame(t(expr_data[,-1]))
names(datExpr) <- expr_data$Gene_ID

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
net <- blockwiseModules(datExpr, power = sft$powerEstimate, TOMType = \"unsigned\", minModuleSize = 30)

module_df <- data.frame(Gene = names(net$colors), Module = net$colors)
write.csv(module_df, output_file, row.names = FALSE)
