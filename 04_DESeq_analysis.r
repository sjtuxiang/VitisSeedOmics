# 1. Check and load required packages
packages <- c("getopt", "DESeq2", "ggplot2", "pheatmap", "RColorBrewer", "grid")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("ERROR: Required package '", pkg, "' is not installed. Install it with install.packages('", pkg, "') or BiocManager::install('", pkg, "').", sep=""))
  }
}

# 2. Define command-line options
spec <- matrix(c(
  'help', 'h', 0, "logical", "Display help information",
  'input', 'i', 1, "character", "Input count matrix file (required)",
  'group', 'g', 1, "character", "Sample group file (required; first column: sample ID, second column: group)",
  'fpkm', 'k', 1, "character", "Input FPKM/TPM matrix (required; used for the correlation plot)",
  'ref', 'r', 1, "character", "Reference/control group name (required; e.g., Control)",
  'fdr', 'f', 1, "double", "FDR threshold [default: 0.05]",
  'fc', 'c', 1, "double", "Fold-change threshold [default: 2]",
  'outdir', 'o', 1, "character", "Output directory [default: current directory]",
  'prefix', 'p', 1, "character", "Output filename prefix [default: Result]"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

# Print usage information
print_usage <- function() {
  cat(getopt(spec, usage = TRUE))
  cat("\nUsage example:\n")
  cat("Rscript DESeq2_OneStep.R -i all_gene_count.tsv -g different_group.txt -k all_gene_tpm.tsv -r LL-CK -p LL-Pv_vs_LL-CK\n")
  q(status = 1)
}

if (!is.null(opt$help)) print_usage()
if (is.null(opt$input) || is.null(opt$group) || is.null(opt$ref) || is.null(opt$fpkm)) print_usage()

# Set default values
if (is.null(opt$fdr)) opt$fdr <- 0.05
if (is.null(opt$fc)) opt$fc <- 2
if (is.null(opt$outdir)) opt$outdir <- getwd()
if (is.null(opt$prefix)) opt$prefix <- "Result"

# Create the output directory
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# ==============================================================================
# 3. Plotting functions (integrated from the original deg_plot_funcs.r)
# ==============================================================================

# Custom colors
mycol <- c("#E41A1C", "#999999", "#377EB8") # Up (red), Normal (gray), Down (blue)

# Draw an MA plot
plot_MA_gg <- function(res_df, fdr_cut, fc_cut, title) {
  p <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = Regulation)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = mycol[1], "Normal" = mycol[2], "Down" = mycol[3])) +
    geom_hline(yintercept = c(log2(fc_cut), -log2(fc_cut)), linetype = "dashed", color = "black") +
    labs(x = "log10(BaseMean + 1)", y = "log2 Fold Change", title = title) +
    theme_bw() + theme(panel.grid = element_blank())
  return(p)
}

# Draw a volcano plot
plot_Volcano_gg <- function(res_df, fdr_cut, fc_cut, title) {
  # Replace zero adjusted P-values to avoid infinite -log10 values
  res_df$padj[res_df$padj == 0] <- min(res_df$padj[res_df$padj > 0]) * 0.1
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = mycol[1], "Normal" = mycol[2], "Down" = mycol[3])) +
    geom_vline(xintercept = c(log2(fc_cut), -log2(fc_cut)), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "gray") +
    labs(x = "log2 Fold Change", y = "-log10(FDR)", title = title) +
    theme_bw() + theme(panel.grid = element_blank())
  return(p)
}

# Draw a correlation plot
plot_corr_gg <- function(fpkm_data, group_info, out_file) {
  # Extract mean FPKM values for the two groups
  groups <- unique(group_info$Group)
  g1 <- groups[1]
  g2 <- groups[2]
  
  samples_g1 <- group_info$ID[group_info$Group == g1]
  samples_g2 <- group_info$ID[group_info$Group == g2]
  
  # Calculate group means
  mean_g1 <- rowMeans(fpkm_data[, samples_g1, drop=FALSE])
  mean_g2 <- rowMeans(fpkm_data[, samples_g2, drop=FALSE])
  
  df_cor <- data.frame(x = log10(mean_g1 + 1), y = log10(mean_g2 + 1))
  pearson_r <- round(cor(mean_g1, mean_g2), 4)
  
  p <- ggplot(df_cor, aes(x = x, y = y)) +
    geom_point(size = 0.5, color = "blue", alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype="dashed") +
    annotate("text", x = min(df_cor$x), y = max(df_cor$y), 
             label = paste0("Pearson r = ", pearson_r), hjust=0, vjust=1, size=5) +
    labs(x = paste0("log10(", g1, " FPKM + 1)"), 
         y = paste0("log10(", g2, " FPKM + 1)"), 
         title = "Group Correlation") +
    theme_bw()
  
  ggsave(paste0(out_file, ".png"), p, width = 6, height = 6, dpi = 300)
  ggsave(paste0(out_file, ".pdf"), p, width = 6, height = 6)
}

# ==============================================================================
# 4. Data processing and main analysis
# ==============================================================================

cat(">>> Reading input data...\n")

# Read sample group information
sample_info <- read.table(opt$group, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(sample_info) <- c("ID", "Group") # Standardize column names

# Read the count matrix
counts_raw <- read.table(opt$input, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")

# Read the FPKM/TPM matrix for result annotation and correlation analysis
fpkm_raw <- read.table(opt$fpkm, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")

# Verify that sample names match across input files
sample_ids <- as.character(sample_info$ID)
count_cols <- colnames(counts_raw)

# Check whether all samples are present in the count matrix
missing_samples <- setdiff(sample_ids, count_cols)
if (length(missing_samples) > 0) {
  cat("\n================ ERROR ================\n")
  cat("Sample IDs from the group file were not found in the count matrix.\n")
  cat("Missing sample IDs: ", paste(missing_samples, collapse = ", "), "\n")
  cat("Columns in the count matrix: ", paste(head(count_cols), collapse = ", "), "...\n")
  cat("Check the file delimiter and whether special characters in sample names were modified (for example, '-' converted to '.').\n")
  cat("=======================================\n")
  stop("Sample names do not match. Analysis terminated.")
}

# Retain samples listed in the group file and match their order
counts_data <- counts_raw[, sample_ids]
fpkm_data <- fpkm_raw[, sample_ids]

# Round count values because DESeq2 requires integer counts
counts_data <- round(counts_data)

# Filter genes with very low total counts
keep <- rowSums(counts_data) >= 2
counts_data <- counts_data[keep, ]
fpkm_data <- fpkm_data[keep, ] # Apply the same filter to the FPKM/TPM matrix

cat(">>> Data validation passed. Samples:", ncol(counts_data), " Genes:", nrow(counts_data), "\n")

# Construct the DESeq2 dataset
colData <- data.frame(row.names = sample_info$ID, condition = factor(sample_info$Group))

# Verify that the reference group exists
if (!opt$ref %in% levels(colData$condition)) {
  stop(paste("ERROR: The specified reference group (-r) '", opt$ref, "' is not present in the group file. Available groups: ", paste(levels(colData$condition), collapse=","), sep=""))
}

# Set the reference level
colData$condition <- relevel(colData$condition, ref = opt$ref)

cat(">>> Running DESeq2 analysis...\n")
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ condition)
dds <- DESeq(dds)

# Extract differential expression results
res <- results(dds)
res_df <- as.data.frame(res)

# Assign regulation categories
log2fc_cut <- log2(opt$fc)
res_df$Regulation <- "Normal"
res_df$Regulation[which(res_df$padj < opt$fdr & res_df$log2FoldChange > log2fc_cut)] <- "Up"
res_df$Regulation[which(res_df$padj < opt$fdr & res_df$log2FoldChange < -log2fc_cut)] <- "Down"
res_df$Regulation <- factor(res_df$Regulation, levels = c("Up", "Down", "Normal"))

# Build the final table containing counts, FPKM/TPM values, and statistics
# Columns: ID | Counts... | FPKM/TPM... | log2FC | pvalue | padj | Regulation
final_table <- cbind(GeneID = rownames(res_df),
                     counts_data[rownames(res_df), ],
                     fpkm_data[rownames(res_df), ],
                     res_df[, c("log2FoldChange", "pvalue", "padj", "Regulation")])

# Rename columns to distinguish count and FPKM/TPM values
new_colnames <- c("GeneID", 
                  paste0(colnames(counts_data), "_Count"),
                  paste0(colnames(fpkm_data), "_FPKM/TPM"),
                  "log2FoldChange", "pvalue", "padj", "Regulation")
colnames(final_table) <- new_colnames

# ==============================================================================
# 5. Output results
# ==============================================================================
out_prefix <- paste0(opt$outdir, "/", opt$prefix)

cat(">>> Writing result files...\n")

# 1. Write the complete results table
write.table(final_table, file = paste0(out_prefix, ".all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 2. Write the differentially expressed gene table
deg_table <- final_table[final_table$Regulation != "Normal", ]
write.table(deg_table, file = paste0(out_prefix, ".DEGs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Count upregulated and downregulated genes
n_up <- sum(final_table$Regulation == "Up", na.rm = TRUE)
n_down <- sum(final_table$Regulation == "Down", na.rm = TRUE)
cat(paste0("Differential expression analysis completed: ", n_up, " upregulated genes and ", n_down, " downregulated genes.\n"))

if ((n_up + n_down) == 0) {
  warning("No genes passed the specified differential expression thresholds; some plots may be affected.")
}

# 3. Generate plots
cat(">>> Generating plots...\n")

# Volcano plot
p_vol <- plot_Volcano_gg(res_df, opt$fdr, opt$fc, paste0("Volcano Plot: ", opt$prefix))
ggsave(paste0(out_prefix, ".Volcano.png"), p_vol, width = 6, height = 5, dpi = 300)
ggsave(paste0(out_prefix, ".Volcano.pdf"), p_vol, width = 6, height = 5)

# MA plot
p_ma <- plot_MA_gg(res_df, opt$fdr, opt$fc, paste0("MA Plot: ", opt$prefix))
ggsave(paste0(out_prefix, ".MA.png"), p_ma, width = 6, height = 6, dpi = 300)
ggsave(paste0(out_prefix, ".MA.pdf"), p_ma, width = 6, height = 6)

# Correlation plot based on group means
plot_corr_gg(fpkm_data, sample_info, paste0(out_prefix, ".Correlation"))

cat(">>> Analysis completed successfully. Results were saved to:", opt$outdir, "\n")
