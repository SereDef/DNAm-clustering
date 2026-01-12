# ==============================================================================
# Post-processing: plot clustering results 
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

data_desc <- "~/MPSR/data/batch_corrected/mega_ComBat.desc"
input_dir <- "~/MPSR/metadata/clusters"
output_dir <- "~/MPSR/metadata/plots"

source("3.0.postprocess_helpers.R")

# Inspect tuning procedure -----------------------------------------------------

phase1_tuning <- readRDS(file.path(input_dir, 'phase1_tuning.rds'))
phase2_tuning <- readRDS(file.path(input_dir, 'phase2_tuning.rds'))

message("Number of clusters, phase 1:")
table(phase1_tuning$final_k)

message("Number of clusters, phase 2:")
phase2_k <- names(which.max(phase2_tuning))
cat(phase2_k)

tuning_plot(phase1_tuning)
# phase2_tuning

# Inspect phase-1 clusters -----------------------------------------------------

phase1_centroids <- readRDS(file.path(input_dir, 'phase1_centroids.rds'))

# Plot all phase 1 centroids (grouped by range centiles) 
centroid_plot(get_empirical_pdf(phase1_centroids), # turn list of ECDF into empirical PDF
              output_file = file.path(output_dir, 'phase1_centroids.pdf')) 

# Plot all CPG data, colored by phase-1 cluster
cpg_data <- bigmemory::attach.big.matrix(data_desc, lockfile = TRUE)

phase1_clusters <- readRDS(file.path(input_dir, 'phase1_clusters.rds'))

cluster_data <- phase1_clusters |> 
  tidyr::separate(cluster, into = c("range_centile", "median_centile", "cluster"), 
                  sep  = "\\.")

clusters_plot(cpg_data, cluster_data,
              output_file = file.path(output_dir, "phase1_cpg_by_cluster.pdf"),
              fixed_x_range = FALSE)
  
# Inspect phase-2 clusters -----------------------------------------------------

phase2_centroids <- readRDS(file.path(input_dir, 'phase2_centroids.rds'))

centroid_plot(get_empirical_pdf(phase2_centroids, type='phase2'),
              output_file = file.path(output_dir, 'phase2_centroids.pdf'))

# plot(table(cl$p2_cluster), col = '#048503', ylab = 'CpG counts')
