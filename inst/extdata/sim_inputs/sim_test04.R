# Load required library
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("Please install DESeq2 to run this script: BiocManager::install('DESeq2')")
}
library(DESeq2)

# --- Simulation Rationale ---
# This script generates a synthetic RNA-seq dataset using the
# DESeq2::makeExampleDESeqDataSet function. It is intended to create a
# standard test case with known differential expression signals, which can be
# used to validate the nexodiff pipeline.
#
# Key features:
# - 10,000 genes and 6 samples.
# - Two groups (A and B), 3 replicates each.
# - A single batch.
# - Non-zero beta standard deviation (betaSD=1) to introduce true
#   log2 fold changes between groups.
# - Randomly generated size factors to simulate library size differences.
# - The true beta values (log2 fold changes) and size factors are exported
#   for later validation.
# - A random seed is used for reproducibility.

# Set seed for reproducibility
set.seed(123)

# Define output directory for the "test04" dataset
output_dir <- "inst/extdata/sim_inputs/test04/config"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Simulation Parameters ---
n_genes <- 10000
n_samples_per_batch <- 6
beta_sd <- 1 # Standard deviation for non-intercept betas (log2 fold changes)

# --- Generate Data for Batch 1 ---
# Generate random size factors for batch 1
size_factors1 <- runif(n_samples_per_batch, min = 0.8, max = 1.2)

# Generate DESeqDataSet for batch 1
dds1 <- makeExampleDESeqDataSet(
  n = n_genes,
  m = n_samples_per_batch,
  betaSD = beta_sd,
  sizeFactors = size_factors1
)

# --- Generate Data for Batch 2 ---
# Generate random size factors for batch 2
size_factors2 <- runif(n_samples_per_batch, min = 0.9, max = 1.3)

# Generate DESeqDataSet for batch 2
dds2 <- makeExampleDESeqDataSet(
  n = n_genes,
  m = n_samples_per_batch,
  betaSD = beta_sd,
  sizeFactors = size_factors2
)


# --- Assemble and write data to CSV files ---

# Define sample names for clarity, including batch identifiers
sample_names_b1 <- paste0("b1_", rep(c("A", "B"), each = n_samples_per_batch / 2), 1:(n_samples_per_batch / 2))
sample_names_b2 <- paste0("b2_", rep(c("C", "D"), each = n_samples_per_batch / 2), 1:(n_samples_per_batch / 2))
all_sample_names <- c(sample_names_b1, sample_names_b2)

# 1. Create and write tx_raw.csv (raw counts from both batches)
counts_df1 <- as.data.frame(counts(dds1))
colnames(counts_df1) <- sample_names_b1
counts_df2 <- as.data.frame(counts(dds2))
colnames(counts_df2) <- sample_names_b2
counts_df <- cbind(counts_df1, counts_df2)
write.table(
  counts_df,
  file.path(output_dir, "tx_raw.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
)

# 2. Create and write design.csv
design <- data.frame(
  sample = all_sample_names,
  batch = c(rep("b1", n_samples_per_batch), rep("b2", n_samples_per_batch)),
  b_label = c(rep("Batch 1", n_samples_per_batch), rep("Batch 2", n_samples_per_batch)),
  group = c(as.character(dds1$condition), as.character(dds2$condition)),
  g_label = c(rep(c("Group A", "Group B"), each = n_samples_per_batch / 2), rep(c("Group C", "Group D"), each = n_samples_per_batch / 2)),
  ctrl = c(rep(c(TRUE, FALSE), each = n_samples_per_batch / 2), rep(c(TRUE, FALSE), each = n_samples_per_batch / 2))
)
write.table(
  design,
  file.path(output_dir, "design.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
)

# 3. Create and write tx_len.csv (dummy transcript lengths)
tx_len <- as.data.frame(matrix(1000, nrow = n_genes, ncol = length(all_sample_names)))
colnames(tx_len) <- all_sample_names
write.table(
  tx_len,
  file.path(output_dir, "tx_len.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
)

# 4. Create and write annotation.csv (dummy annotation)
annotation <- data.frame(
  gene_id = 1:n_genes,
  transcript_id = paste0("t", 1:n_genes),
  biotype = "mRNA",
  tax_id = 9606,
  species = "human"
)

write.table(
  annotation,
  file.path(output_dir, "annotation.csv"),
  sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE
)

id_mapping <- data.frame(
  "Entry" = paste0("UN", 1:n_genes),
  "Status" = rep("reviewed", n_genes),
  "Protein names" = paste0("name for ", 1:n_genes),
  "Gene names" = paste0("SYM", 1:n_genes),
  "Cross-reference (GeneID)" = paste(1:n_genes, ";", sep = ""),
  "Annotation" = rep("5 out of 5", n_genes),
  check.names = FALSE
)

write.table(
  id_mapping,
  file.path(output_dir, "idmapping.tab"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
)

# 5. Create and write size_factors.csv
# Stored as a single row with samples as columns for consistency
all_size_factors <- c(size_factors1, size_factors2)
sf_df <- data.frame(t(all_size_factors))
colnames(sf_df) <- all_sample_names
write.table(
  sf_df,
  file.path(output_dir, "size_factors.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
)

# 6. Create and write true_betas.csv (true log2 fold changes for both batches)
betas_df <- data.frame(
  gene_id = paste0("g", 1:n_genes),
  true_beta_b1 = mcols(dds1)$trueBeta,
  true_beta_b2 = mcols(dds2)$trueBeta
)
write.table(
  betas_df,
  file.path(output_dir, "true_betas.csv"),
  sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
)

print(paste("Generated simulation files for test04 in:", file.path(getwd(), output_dir)))