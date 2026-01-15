# ==============================================================================
# EP-means clustering - stage 1
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

# I process each centile as a task
array_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID', unset = 1))
# I use CPUs to run iterations in parallel
n_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', unset = 1))

# Set-up
k_candidates = 2:10 # K values to try
feature_split = 0.5 # Proportion of CpGs to use for tuning at each iteration
n_iterations = 10   # Number of tuning iterations 

data_descfile <- "~/MPSR/data/batch_corrected/mega_ComBat.desc"
metadata_file <- "~/MPSR/DNAm-clustering/metadata/CpG_metadata.rds"
output_dir <- "~/MPSR/metadata/clusters"

# Meta-data matrix
meta_bins <- readRDS(metadata_file)
centiles <- sort(unique(meta_bins[,'centile_bins'])) # 100 in this version

centile <- centiles[array_id]

centile_subset <- meta_bins[meta_bins[,'centile_bins'] == centile, ]

# Attach (pre-processed) methylation data matrix[rows = CpGs, cols = samples]
data <- bigmemory::attach.big.matrix(data_descfile, lockfile = TRUE)

# Select relevant cpgs (from the centile)
df <- data[rownames(centile_subset), ]

message('=====================================================================')
message("Centile: ", centile)
message("\nData median: ")
print(summary(as.numeric(centile_subset[,'median_mega'])))
message("\nData range: ")
print(summary(as.numeric(centile_subset[,'range_mega'])))

message("\nSubset data size: ")
print(dim(df))

# drop CpGs with any missing values (TMP)
df <- df[complete.cases(df),]

message("\nActual data size (drop NA): ")
print(dim(df))
message('=====================================================================')

# Pre-compute CpG ECDFs 
n_cpg = nrow(df)
ecdf_list <- lapply(seq_len(n_cpg), function(cpg_idx) stats::ecdf(df[cpg_idx, ]))
names(ecdf_list) <- rownames(df)

rm(df, data, meta_bins, centiles, centile_subset)

# ==============================================================================
source("2.0.epmeans_helpers.R")

k_stability <- stability_selection(ecdf_list,
                                   k_values = k_candidates,
                                   beta = feature_split,
                                   max_iter = n_iterations,
                                   seed = 3108,
                                   n_cores = n_cores)

message('=====================================================================')
message("Stability estimates: ")
print(k_stability)
message('=====================================================================')

k_hat <- as.integer(names(which.max(k_stability)))

# Handle edge case with more than 1 k value (parsimony)
if (length(k_hat) > 1) k_hat <- min(k_hat)

# ==============================================================================
saveRDS(t(k_stability),
        file = file.path(output_dir, paste0('c',centile,'_tuning.rds')))

final_cluster <- fit_epmeas(ecdf_list, k = k_hat, centile = centile, 
                            output_folder = output_dir)

