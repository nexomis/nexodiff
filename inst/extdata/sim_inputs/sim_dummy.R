# --- Simulation Rationale ---
# This script generates a synthetic RNA-seq dataset designed to test the
# robustness of inter-sample normalization methods.
#
# Key features:
# - A 'block' of genes with varying expression levels is used as a base.
# - Control groups (A, D) represent a stable baseline.
# - Test groups (B, C, E, F) introduce compositional biases:
#   - Some samples have a subset of genes that are strongly up-regulated (e.g., B2, E).
#   - Some samples have a subset of genes that are down-regulated (e.g., B3, F).
#
# The design is intentionally BALANCED: the number of samples with up-regulation
# bias is matched by samples with down-regulation bias. This is crucial for
# testing normalization methods like median-of-ratios (used in DESeq2, edgeR)
# which assume that for `ref_type="all"`, most genes are not changing, or that
# up- and down-regulation are symmetric across the experiment. An unbalanced
# design would skew the reference and lead to incorrect normalization factors.

# Create a directory for the simulation output
output_dir <- "inst/extdata/sim_inputs/test03/config"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Base block for generating counts. Represents genes with a wide dynamic range.
block <- c(
  0:10 * 8,
  (0:10) * 80,
  (0:10) * 800
)

# --- Expected Normalization Behavior Notes ---
# These comments outline the expected outcomes for different normalization strategies.
# They are for documentation and not used in the script directly.

# --- Generate sample counts based on the initial script ---

# Batch 1, Group A (3 reps, control)
# These represent the baseline control group.
# A1 and A2 have different library sizes, which should be corrected by TPM.
A1_counts <- c(block, block, block) / 2
A2_counts <- c(block, block, block) * 2
A3_counts <- c(block, block, block)

# Batch 1, Group B (3 reps)
# This group tests various compositional biases.
# median ratio = 1 (balanced up- and down-regulation within the sample)
B1_counts <- c(0.5 * block, block, 1.5 * block)
# median ratio > 1 (biased towards up-regulation)
B2_counts <- c(block, block, 1.5 * block)
# median ratio < 1 (biased towards down-regulation)
# The 0.625 factor is calculated to balance the up-regulation in other samples.
# x = 1/ (1 - (1/3.5) - (1/3))
B3_counts <- c(block, block, 0.625 * block)

# Batch 1, Group C (3 reps)
# Similar to Group B but with different library sizes.
# median ratio = 1
C1_counts <- c(1.5 * block, block, 0.5 * block) * 2
# median ratio > 1
C2_counts <- c(block, block, block * 1.5) * 2
# median ratio < 1 (CHANGED TO BE BALANCED)
# This sample is biased towards down-regulation to help balance the overall design.
C3_counts <- c(block, block, block * 0.625) / 2

# Batch 2, Group D (3 reps, control)
# A second, stable control group in a different batch.
# NOTE: Replicates are identical as per the original script's structure
D1_counts <- c(block, block, block)
D2_counts <- c(block, block, block)
D3_counts <- c(block, block, block)

# Batch 2, Group E (3 reps)
# This group is strongly biased towards up-regulation.
# NOTE: Replicates are identical as per the original script's structure
E1_counts <- c(block, block, 1.5 * block)
E2_counts <- c(block, block, 1.5 * block)
E3_counts <- c(block, block, 1.5 * block)

# Batch 2, Group F (3 reps, NEW GROUP TO BALANCE THE DESIGN)
# This group is strongly biased towards down-regulation. It is added to
# counteract the up-regulation bias from Group E, making the overall
# experimental design symmetric.
F1_counts <- c(block, block, 0.625 * block)
F2_counts <- c(block, block, 0.625 * block)
F3_counts <- c(block, block, 0.625 * block)


# --- Assemble and write data to CSV files ---

# 1. Create and write tx_raw.csv
tx_raw <- data.frame(
  A1 = A1_counts, A2 = A2_counts, A3 = A3_counts,
  B1 = B1_counts, B2 = B2_counts, B3 = B3_counts,
  C1 = C1_counts, C2 = C2_counts, C3 = C3_counts,
  D1 = D1_counts, D2 = D2_counts, D3 = D3_counts,
  E1 = E1_counts, E2 = E2_counts, E3 = E3_counts,
  F1 = F1_counts, F2 = F2_counts, F3 = F3_counts
)
write.table(tx_raw, file.path(output_dir, "tx_raw.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 2. Create and write design.csv
design <- data.frame(
  sample = colnames(tx_raw),
  batch = rep(c("b1", "b2"), times = c(9, 9)),
  group = rep(c("A", "B", "C", "D", "E", "F"), each = 3),
  ctrl = rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE), each = 3)
)
write.table(design, file.path(output_dir, "design.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 3. Create and write tx_len.csv
num_transcripts <- nrow(tx_raw)
num_samples <- ncol(tx_raw)
tx_len <- as.data.frame(matrix(1000, nrow = num_transcripts, ncol = num_samples))
colnames(tx_len) <- colnames(tx_raw)
write.table(tx_len, file.path(output_dir, "tx_len.csv"), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 4. Create and write annotation.csv
annotation <- data.frame(
  gene_id = paste0("g", 1:num_transcripts),
  transcript_id = paste0("t", 1:num_transcripts),
  biotype = "mRNA",
  tax_id = 9606,
  species = "human"
)
write.table(annotation, file.path(output_dir, "annotation.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

print(paste("Generated simulation files in:", file.path(getwd(), output_dir)))

# --- Expected Ratios ---
# This data frame defines the expected BIOLOGICAL ratios for each sample,
# which should be recoverable after proper normalization.
r <- rep(1, 33)
rc <- c(r, r, r)
req <- c(0.5 * r, r, 1.5 * r)
req2 <- c(1.5 * r, r, 0.5 * r)
rup <- c(r, r, 1.5 * r)
rdown <- c(r, r, 0.625 * r)

exp_ratios <- data.frame(
  A1 = rc, A2 = rc, A3 = rc,
  B1 = req, B2 = rup, B3 = rdown,
  C1 = req2, C2 = rup, C3 = rdown,
  D1 = rc, D2 = rc, D3 = rc,
  E1 = rup, E2 = rup, E3 = rup,
  F1 = rdown, F2 = rdown, F3 = rdown
)

write.table(
  exp_ratios, file.path(output_dir, "exp_ratios.csv"), sep = ",",
  row.names = FALSE, col.names = TRUE, quote = FALSE
)

