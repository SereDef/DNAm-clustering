# ==============================================================================
# Post-processing: merge clustering results from stage 1 and 2
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

input_dir <- "~/MPSR/metadata/clusters"
output_dir <- "~/MPSR/metadata/clusters"

phase2_tuning <- readRDS(file.path(input_dir, 'phase2_tuning.rds'))

# Number of clusters, phase 2
phase2_k <- names(which.max(phase2_tuning))

# Save final metadata file =====================================================

# Cluster assignments 
phase1_clusters <- readRDS(file.path(input_dir, 'phase1_clusters.rds'))

# Meta-data matrix
meta <- readRDS(file.path(dirname(input_dir), 'metadata_mega_ComBat_binned.rds'))

metadata <- as.data.frame(meta) |> tibble::rownames_to_column(var = "cpg")

metadata <- merge(metadata, phase1_clusters, by = "cpg")

# Phase 2 (macro) cluster assignment
phase2_clusters <- readRDS(file.path(input_dir, 
                                     paste0('phase2_k', phase2_k, '.rds'))) |>
  as.data.frame() |>
  dplyr::rename(cluster = V1, p2_cluster = V2)

metadata <- merge(metadata, phase2_clusters, by = "cluster") |>
  dplyr::select(cpg, # centile_bin = centile_bins, 
                range_cpg = range_mega, median_cpg = median_mega,
                p1_cluster = cluster, p2_cluster)

saveRDS(metadata, file.path(output_dir, 'metadata.rds'))

