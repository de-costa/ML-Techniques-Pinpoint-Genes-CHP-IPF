# ============================================================
# STEP 01: DATA ACQUISITION AND PROCESSING
# Dataset: GSE150910
# File: GSE150910_gene-level_count_file.csv.gz
# Tools: GEOquery, DESeq2, limma
# ============================================================

# ── SECTION 1: LOAD LIBRARIES ──────────────────────────────

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "DESeq2", "limma"),
                     update = FALSE)

install.packages(c("dplyr", "ggplot2", "pheatmap"),
                 dependencies = TRUE)

# Load libraries
library(GEOquery)   # Download GEO metadata
library(DESeq2)     # Normalization and VST
library(limma)      # Batch correction
library(dplyr)      # Data manipulation
library(ggplot2)    # PCA visualization
library(pheatmap)   # Heatmap visualization

cat("All libraries loaded successfully!\n")

# ── SECTION 2: LOAD COUNT MATRIX ───────────────────────────

cat("\nLoading gene-level count matrix...\n")

# Load the gene-level count file
# This file already has genes as rows and samples as columns
# No aggregation needed — cleaner than transcript-level file
counts <- read.csv(
  "../../DATASET/GSE150910_gene-level_count_file.csv",
  row.names = 1,    # First column = gene names as row names
  check.names = FALSE)  # Keep original sample names

# Check dimensions
cat("Count matrix dimensions:\n")
cat("  Genes:   ", nrow(counts), "\n")
cat("  Samples: ", ncol(counts), "\n")

# Preview first few rows and columns
cat("\nFirst 5 genes, first 5 samples:\n")
print(counts[1:5, 1:5])

# ── SECTION 3: DOWNLOAD METADATA ───────────────────────────

cat("\nDownloading metadata from GEO...\n")

# Download series matrix (contains sample information)
gse <- getGEO("GSE150910",
              GSEMatrix  = TRUE,
              AnnotGPL   = FALSE)

# Extract metadata
metadata <- pData(gse[[1]])

cat("Metadata dimensions:", dim(metadata), "\n")

# Check diagnosis groups
cat("\nDiagnosis groups:\n")
print(table(metadata$`diagnosis:ch1`))

# ── SECTION 4: CLEAN METADATA ──────────────────────────────

cat("\nCleaning metadata...\n")

# Create clean metadata table with relevant columns
meta_clean <- data.frame(
  sample_id = metadata$title,
  diagnosis = metadata$`diagnosis:ch1`,
  batch     = metadata$`batch:ch1`,
  sex       = metadata$`Sex:ch1`,
  age       = metadata$`age:ch1`,
  row.names = metadata$title
)

# Check clean metadata
cat("Clean metadata preview:\n")
print(head(meta_clean))

cat("\nDiagnosis distribution:\n")
print(table(meta_clean$diagnosis))

cat("\nBatch groups:\n")
print(table(meta_clean$batch))

# ── SECTION 5: MATCH SAMPLES ───────────────────────────────

cat("\nMatching samples between count matrix and metadata...\n")

# Check if all count matrix samples exist in metadata
cat("All samples found in metadata:",
    all(colnames(counts) %in% rownames(meta_clean)), "\n")

# Reorder metadata to match count matrix column order
meta_clean <- meta_clean[colnames(counts), ]

# Verify perfect alignment
cat("Perfect alignment:",
    all(colnames(counts) == rownames(meta_clean)), "\n")

# ── SECTION 6: FILTER LOW EXPRESSION GENES ─────────────────

cat("\nFiltering low expression genes...\n")

# Round to integers for DESeq2
# Gene-level file may still have decimal values
counts_int <- round(counts)

# Filter: keep genes with at least 10 counts
# in at least 10% of samples
min_samples <- round(0.1 * ncol(counts_int))

cat("Minimum samples threshold:", min_samples, "\n")

keep <- rowSums(counts_int >= 10) >= min_samples
counts_filtered <- counts_int[keep, ]

cat("\nFiltering results:\n")
cat("  Genes before filtering:", nrow(counts_int), "\n")
cat("  Genes after filtering: ", nrow(counts_filtered), "\n")
cat("  Genes removed:         ",
    nrow(counts_int) - nrow(counts_filtered), "\n")

# ── SECTION 7: DESeq2 NORMALIZATION ────────────────────────

cat("\nRunning DESeq2 normalization...\n")

# Create DESeq2 object
# design = ~ diagnosis tells DESeq2 about our groups
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = meta_clean,
  design    = ~ diagnosis)

# Estimate size factors (median-of-ratios normalization)
dds <- estimateSizeFactors(dds)

cat("Size factors (first 10):\n")
print(round(sizeFactors(dds)[1:10], 3))

cat("\nSize factor range:\n")
cat("  Min:", round(min(sizeFactors(dds)), 3), "\n")
cat("  Max:", round(max(sizeFactors(dds)), 3), "\n")

# ── SECTION 8: VST TRANSFORMATION ──────────────────────────

cat("\nApplying Variance Stabilizing Transformation...\n")

# blind=TRUE: don't use group info during transformation
# ensures unbiased quality control
vsd <- vst(dds, blind = TRUE)

# Extract VST matrix
vst_matrix <- assay(vsd)

cat("VST matrix dimensions:\n")
cat("  Genes:   ", nrow(vst_matrix), "\n")
cat("  Samples: ", ncol(vst_matrix), "\n")
cat("  Range:   ", round(min(vst_matrix), 2),
    "to", round(max(vst_matrix), 2), "\n")

# ── SECTION 9: BATCH EFFECT CORRECTION ─────────────────────

cat("\nApplying batch effect correction...\n")

# Convert batch to factor
batch <- as.factor(meta_clean$batch)

cat("Batch groups:\n")
print(table(batch))

# Remove batch effects using limma
# design matrix protects biological signal
vst_corrected <- removeBatchEffect(
  vst_matrix,
  batch  = batch,
  design = model.matrix(~ diagnosis, data = meta_clean))

cat("\nBatch correction complete!\n")
cat("Corrected matrix range:\n")
cat("  Min:", round(min(vst_corrected), 2), "\n")
cat("  Max:", round(max(vst_corrected), 2), "\n")
cat("  NA values:", any(is.na(vst_corrected)), "\n")

# ── SECTION 10: PCA VISUALIZATION ──────────────────────────

cat("\nGenerating PCA plot...\n")

# Run PCA on corrected matrix
# Transpose: PCA expects samples as rows
pca_result <- prcomp(t(vst_corrected), scale. = FALSE)

# Calculate variance explained
pca_var <- round(
  100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Create PCA dataframe
pca_df <- data.frame(
  PC1       = pca_result$x[, 1],
  PC2       = pca_result$x[, 2],
  diagnosis = meta_clean$diagnosis,
  batch     = meta_clean$batch)

# Plot PCA
pca_plot <- ggplot(pca_df,
                   aes(PC1, PC2, color = diagnosis)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "PCA — GSE150910 After Batch Correction",
    x     = paste0("PC1 (", pca_var[1], "% variance)"),
    y     = paste0("PC2 (", pca_var[2], "% variance)")) +
  scale_color_manual(
    values = c("chp"     = "#E74C3C",
               "control" = "#2ECC71",
               "ipf"     = "#3498DB")) +
  theme_bw() +
  theme(legend.title = element_text(face = "bold"))

# Save PCA plot
ggsave(
  "D:/Research/implementation/step01/PCA_plot.png",
  plot   = pca_plot,
  width  = 8,
  height = 6,
  dpi    = 300)

cat("PCA plot saved!\n")

# ── SECTION 11: SAMPLE DISTANCE HEATMAP ────────────────────

cat("\nGenerating sample distance heatmap...\n")

# Calculate sample distances
sample_dists <- dist(t(vst_corrected))
dist_matrix  <- as.matrix(sample_dists)

# Create annotation for heatmap
annotation <- data.frame(
  Diagnosis = meta_clean$diagnosis,
  row.names = rownames(meta_clean))

# Define colors
ann_colors <- list(
  Diagnosis = c(chp     = "#E74C3C",
                control = "#2ECC71",
                ipf     = "#3498DB"))

# Plot and save heatmap
pheatmap(
  dist_matrix,
  annotation_col    = annotation,
  annotation_row    = annotation,
  annotation_colors = ann_colors,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  main              = "Sample Distance Heatmap — GSE150910",
  filename          = "D:/Research/implementation/step01/heatmap.png",
  width             = 10,
  height            = 8)

cat("Heatmap saved!\n")

# ── SECTION 12: SAVE ALL OUTPUT FILES ──────────────────────

cat("\nSaving all output files...\n")

# Save corrected VST matrix (main file for Python)
write.csv(
  vst_corrected,
  "D:/Research/implementation/step01/vst_corrected.csv")

# Save clean metadata
write.csv(
  meta_clean,
  "D:/Research/implementation/step01/meta_final.csv",
  row.names = TRUE)

# Save filtered counts
write.csv(
  counts_filtered,
  "D:/Research/implementation/step01/counts_filtered.csv")

# ── SECTION 13: FINAL SUMMARY ──────────────────────────────

cat("\n")
cat("=" , strrep("=", 45), "\n")
cat("STEP 01 COMPLETE — SUMMARY\n")
cat("=" , strrep("=", 45), "\n")
cat("Dataset:              GSE150910\n")
cat("File used:            gene-level count file\n")
cat("Total samples:        ", ncol(counts), "\n")
cat("Diagnosis groups:\n")
print(table(meta_clean$diagnosis))
cat("Genes before filter:  ", nrow(counts_int), "\n")
cat("Genes after filter:   ", nrow(counts_filtered), "\n")
cat("Genes removed:        ",
    nrow(counts_int) - nrow(counts_filtered), "\n")
cat("VST range:            ",
    round(min(vst_matrix), 2), "to",
    round(max(vst_matrix), 2), "\n")
cat("Corrected range:      ",
    round(min(vst_corrected), 2), "to",
    round(max(vst_corrected), 2), "\n")
cat("\nFiles saved:\n")
cat("  - vst_corrected.csv\n")
cat("  - meta_final.csv\n")
cat("  - counts_filtered.csv\n")
cat("  - PCA_plot.png\n")
cat("  - heatmap.png\n")
cat("=" , strrep("=", 45), "\n")
cat("Ready for Step 02 — Feature Selection!\n")