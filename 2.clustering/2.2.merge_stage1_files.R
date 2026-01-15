# ==============================================================================
# EP-means clustering - merge cluster files from stage 1 (prep for stage 2)
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

input_dir <- output_dir <- "~/MPSR/metadata/clusters"

source("2.0.epmeans_helpers.R")

# require(dplyr)
# require(readr)
# require(purrr)

# Identify the files containing the clusters assigned to each CpG in stage 1 ---
stage1_cluster_files <- list.files(
  path = input_dir,
  pattern = "^c.*_k[1-9][0-9]?\\.rds$",
  full.names = TRUE)

# Combine them into one
stage1_clusters <- purrr::map_dfr(
  stage1_cluster_files,
  ~ readr::read_rds(.x) |>               # read one at the time
    as.data.frame()  |>                  # convert matrix to data.frame
    dplyr::mutate(file = basename(.x)),  # record file name
  .id = NULL) |>  
  # add centile to cluster name 
  dplyr::mutate(centile = sub("^c(.*)_k[0-9]+\\.rds$", "\\1", file),
                cluster = paste(centile, V2, sep = '.')) |>
  dplyr::select(cpg = V1, cluster)

# Save single file
saveRDS(stage1_clusters,
        file.path(output_dir, 'stage1_clusters.rds'))

# Remove individual files
file.remove(stage1_cluster_files)

# Identify the tuning files, and combine them ----------------------------------

stage1_tune_files <- list.files(
  path = input_dir,
  pattern = "^c.*_tuning.rds$",
  full.names = TRUE
)

k_eval <- as.character(2:10) # TODO: make this dynamic based on what was used in tuning

# Combine them into one
stage1_tuning <- purrr::map_dfr(
  stage1_tune_files,
  ~ readr::read_rds(.x) |>               # read one at the time
    as.data.frame()  |>                  # convert matrix to data.frame
    dplyr::mutate(file = basename(.x)),  # record file name
  .id = NULL) |>  
  # add centile to cluster name 
  dplyr::mutate(centile = sub("^c(.*)_tuning.rds$", "\\1", file),
                final_k = as.numeric(
                  names(dplyr::pick(dplyr::all_of(k_eval)))[
                    max.col(dplyr::pick(dplyr::all_of(k_eval)))])) |> 
  dplyr::select(centile, final_k, dplyr::all_of(k_eval))


# save single file
saveRDS(stage1_tuning,
        file.path(output_dir, 'stage1_tuning.rds'))

# Remove individual files
file.remove(stage1_tune_files)

# Centroids --------------------------------------------------------------------

stage1_centroid_files <- list.files(
    path   = input_dir,
    pattern = "^c.*_centroids.rds$",
    full.names = TRUE)

ecdf_list <- stage1_centroid_files |>
  purrr::imap(
    function(path, centile) {
      # Get centile names
      file <- basename(path)
      centile <- sub("^c(.*)_centroids.rds$", "\\1", file)
      
      # Read list of centroid ECDF objects 
      centroids  <- readRDS(path)
      
      # Rename elements: <centile>.<cluster>
      names(centroids) <- paste(centile, seq_along(centroids), sep='.')
      
      return(centroids)}) |>
  # Flatten
  purrr::list_c() 

# save single file
saveRDS(ecdf_list,
        file.path(output_dir, 'stage1_centroids.rds'))

# Remove individual files
file.remove(stage1_centroid_files)

# ==============================================================================
# Also merge all the log files

system('cat 2.1.logs/* > 2.1.log')
system('rm -r 2.1.logs')