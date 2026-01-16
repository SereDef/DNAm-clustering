# ==============================================================================
# Post-processing: merge clustering results from stage 1 and 2
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

input_dir <- "~/MPSR/metadata/clusters"
output_file <- "~/MPSR/DNAm-clustering/metadata/CpG_metadata.rds"

stage2_tuning <- readRDS(file.path(input_dir, 'stage2_tuning.rds'))

# Number of clusters, stage 2
stage2_k <- names(which.max(stage2_tuning))

# Cluster assignments 
stage1_clusters <- readRDS(file.path(input_dir, 'stage1_clusters.rds'))

# Meta-data matrix
metadata <- readRDS(output_file)|> 
  as.data.frame() |> 
  tibble::rownames_to_column(var = "cpg")

metadata <- merge(metadata, stage1_clusters, by = "cpg")

# stage 2 (macro) cluster assignment
stage2_clusters <- readRDS(file.path(input_dir, 
                                     paste0('stage2_k', stage2_k, '.rds'))) |>
  as.data.frame() |>
  dplyr::rename(cluster = V1, p2_cluster = V2)

metadata <- merge(metadata, stage2_clusters, by = "cluster") |>
  dplyr::select(
    cpg,
    p1_cluster = cluster,
    p2_cluster,
    matches("^(mean|median|range|variance|skewness|kurtosis)_")
  )

# Save final metadata file =====================================================

# Overwrite file 
saveRDS(metadata, output_file)
