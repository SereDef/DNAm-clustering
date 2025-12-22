# ==============================================================================
# EP-means clustering - stage 2
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

# I process each centile as a task
# array_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID', unset = 1))
# I use CPUs to run iterations in parallel
n_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', unset = 1))

# Set-up
k_candidates = 2:20 # K values to try
feature_split = 0.25 # Proportion of CpGs to use for tuning at each iteration
n_iterations = 10   # Number of tuning iterations 

input_dir <- "~/MPSR/metadata/clusters"

ecdf_list <- readRDS(file.path(input_dir, 'phase1_centroids.rds'))

# ==============================================================================
source("2.0.epmeans_helpers.R")

k_stability <- stability_selection(ecdf_list,
                                   k_values = k_candidates,
                                   beta = feature_split,
                                   max_iter = n_iterations,
                                   seed = 3108,
                                   n_cores = n_cores)

message("Stability estimates: ")
print(k_stability)

k_hat <- as.integer(names(which.max(k_stability)))

# Handle edge case with more than 1 k value (parsimony)
if (length(k_hat) > 1) k_hat <- min(k_hat)

# ==============================================================================
saveRDS(k_stability, file.path(input_dir, 'phase2_tuning.rds'))

final_cluster <- fit_epmeas(ecdf_list, k = k_hat, centile = 'phase2', output_folder = input_dir)

