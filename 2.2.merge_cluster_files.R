# ==============================================================================
# EP-means clustering - merge cluster files from stage 1 (prep for stage 2)
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

input_dir <- "~/MPSR/metadata/clusters"
output_dir <- "~/MPSR/metadata/clusters"

source("2.0.epmeans_helpers.R")
# library(dplyr)
# require(readr)
# require(purrr)
# require(tidyr)

# Identify the files containing the clusters assigned to each CpG in phase 1 ---
phase1_cluster_files <- list.files(
  path   = input_dir,
  pattern = "^c.*_k[1-9][0-9]?\\.rds$",
  full.names = TRUE
)

# Combine them into one
phase1_clusters <- purrr::map_dfr(
  phase1_cluster_files,
  ~ readr::read_rds(.x) |>               # read one at the time
    as.data.frame()  |>                  # convert matrix to data.frame
    dplyr::mutate(file = basename(.x)),  # record file name
  .id = NULL) |>  
  # add centile to cluster name 
  dplyr::mutate(centile = sub("^c(.*)_k[0-9]+\\.rds$", "\\1", file),
                cluster = paste(centile, V2, sep = '.')) |>
  dplyr::select(cpg = V1, cluster)

# Save single file
saveRDS(phase1_clusters,
        file.path(input_dir, 'phase1_clusters.rds'))

# Remove individual files
file.remove(phase1_cluster_files)

# Identify the tuning files, and combine them ----------------------------------

phase1_tune_files <- list.files(
  path   = input_dir,
  pattern = "^c.*_tuning.rds$",
  full.names = TRUE
)

# Combine them into one
phase1_tuning <- purrr::map_dfr(
  phase1_tune_files,
  ~ readr::read_rds(.x) |>               # read one at the time
    as.data.frame()  |>                  # convert matrix to data.frame
    dplyr::mutate(file = basename(.x)),  # record file name
  .id = NULL) |>  
  # add centile to cluster name 
  dplyr::mutate(centile = sub("^c(.*)_tuning.rds$", "\\1", file),
                final_k = as.numeric(
                  names(dplyr::pick(`2`:`10`))[max.col(dplyr::pick(`2`:`10`))])) |> 
  dplyr::select(centile, final_k, `2`:`10`)


# save single file
saveRDS(phase1_tuning,
        file.path(input_dir, 'phase1_tuning.rds'))

# Remove individual files
file.remove(phase1_tune_files)

# Centroids --------------------------------------------------------------------

phase1_centroid_files <- list.files(
    path   = input_dir,
    pattern = "^c.*_centroids.rds$",
    full.names = TRUE
  )

ecdf_list <- phase1_centroid_files |>
  purrr::imap(
    function(path, centile) {
      # Get centile names
      file <- basename(path)
      centile <- sub("^c(.*)_centroids.rds$", "\\1", file)
      
      # Read in centrod list 
      centroids  <- readRDS(path) # list of centroid ecdf
      
      # rename elements: <centile>.<cluster>
      names(centroids) <- paste(centile, seq_along(centroids), sep='.')
      
      return(centroids)}) |>
  # flatten to one list
  purrr::list_c() 


# save single file
saveRDS(ecdf_list,
        file.path(input_dir, 'phase1_centroids.rds'))

# Remove individual files
file.remove(phase1_centroid_files)

# ==============================================================================

