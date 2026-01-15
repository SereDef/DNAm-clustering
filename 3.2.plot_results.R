# ==============================================================================
# Post-processing: plot clustering results 
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

data_desc <- "~/MPSR/data/batch_corrected/mega_ComBat.desc"
input_dir <- "~/MPSR/metadata/clusters"
output_dir <- "~/MPSR/metadata/plots"

dir.create(output_dir, showWarnings = FALSE)

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

# Inspect phase-1 and 2 centroids ----------------------------------------------

phase1_centroids <- readRDS(file.path(input_dir, 'phase1_centroids.rds'))

# Plot all phase 1 centroids (grouped by range centiles) 
centroid_plot(get_empirical_pdf(phase1_centroids), # turn list of ECDF into empirical PDF
              output_file = file.path(output_dir, 'phase1_centroids.pdf')) 

phase2_centroids <- readRDS(file.path(input_dir, 'phase2_centroids.rds'))

centroid_plot(get_empirical_pdf(phase2_centroids, type='phase2'),
              output_file = file.path(output_dir, 'phase2_centroids.pdf'))

# Inspect phase-1 and 2 clusters -----------------------------------------------

cpg_data <- bigmemory::attach.big.matrix(data_desc, lockfile = TRUE)

cluster_data <- readRDS('metadata.rds')

clusters_plot(cpg_data, cluster_data, cluster_var = 'p1_cluster', 
              output_file = file.path(output_dir, "phase1_cpg_by_cluster.pdf"),
              fixed_x_range = FALSE)

clusters_plot(cpg_data, cluster_data, cluster_var = 'p2_cluster', 
              output_file = file.path(output_dir, "phase2_cpg_by_cluster.pdf"),
              fixed_x_range = TRUE)

# EWAS functional analyses -----------------------------------------------------

sumstat_data <- readRDS('summstats_prenatalrisk.rds')

dset <- merge(sumstat_data, cluster_data, by = "cpg", all.x = TRUE)

rm(cluster_data, sumstat_data)

ewases <- gsub('pvalue_', '', grep('^pvalue_', names(dset), value = TRUE))

ewas_plot(dset, ewases, cluster_var = 'p2_cluster', 
          thresh_gnmwide = 1e-7, 
          output_file = file.path(output_dir, 'EWAS_cluster_repr.pdf'))

