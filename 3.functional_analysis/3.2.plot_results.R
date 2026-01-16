# ==============================================================================
# Post-processing: plot clustering results 
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

data_desc <- "~/MPSR/data/batch_corrected/mega_ComBat.desc"
metadata_dir <- "~/MPSR/DNAm-clustering/metadata"
input_dir <- "~/MPSR/metadata/clusters"
output_dir <- "~/MPSR/DNAm-clustering/plots"

dir.create(output_dir, showWarnings = FALSE)

source("3.0.postprocess_helpers.R")

# Inspect tuning procedure -----------------------------------------------------

stage1_tuning <- readRDS(file.path(input_dir, 'stage1_tuning.rds'))
stage2_tuning <- readRDS(file.path(input_dir, 'stage2_tuning.rds'))

message("Number of clusters, stage 1:")
table(stage1_tuning$final_k)

message("Number of clusters, stage 2:")
stage2_k <- names(which.max(stage2_tuning))
cat(stage2_k)

tuning_plot(stage1_tuning)
# stage2_tuning

# Inspect stage-1 and 2 centroids ----------------------------------------------

stage1_centroids <- readRDS(file.path(input_dir, 'stage1_centroids.rds'))

# Plot all stage 1 centroids (grouped by range centiles) 
centroid_plot(get_empirical_pdf(stage1_centroids, smoothing = 2),
              output_file = file.path(output_dir, 'stage1_centroids.pdf')) 

stage2_centroids <- readRDS(file.path(input_dir, 'stage2_centroids.rds'))

centroid_plot(get_empirical_pdf(stage2_centroids, type='stage2', smoothing = 3),
              output_file = file.path(output_dir, 'stage2_centroids.pdf'))

# Inspect stage-1 and 2 clusters -----------------------------------------------

cpg_data <- bigmemory::attach.big.matrix(data_desc, lockfile = TRUE)

cluster_data <- readRDS(file.path(metadata_dir, 'CpG_metadata.rds'))

clusters_plot(cpg_data, cluster_data, cluster_var = 'p1_cluster', 
              random_subset = 0.5,
              output_file = file.path(output_dir, "stage1_cpg_by_cluster.pdf"),
              fixed_x_range = FALSE)

clusters_plot(cpg_data, cluster_data, cluster_var = 'p2_cluster', 
              random_subset = 0.5,
              output_file = file.path(output_dir, "stage2_cpg_by_cluster.pdf"),
              fixed_x_range = TRUE)

# EWAS functional analyses -----------------------------------------------------

sumstat_data <- readRDS(file.path(metadata_dir, 'summstats_prenatalrisk.rds'))

dset <- merge(sumstat_data, cluster_data, by = "cpg", all.x = TRUE)

rm(cluster_data, sumstat_data)

ewases <- gsub('pvalue_', '', grep('^pvalue_', names(dset), value = TRUE))

ewas_plot(dset, ewases, cluster_var = 'p2_cluster', 
          thresh_gnmwide = 1e-7, 
          output_file = file.path(output_dir, 'EWAS_cluster_repr.pdf'))

