# ==============================================================================
# EP-means clustering - merge cluster files from stage 1 (prep for stage 2)
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

input_dir <- "~/MPSR/metadata/clusters"

# Meta-data matrix
meta_bins <- readRDS(file.path(dirname(input_dir), 'metadata_mega_ComBat_binned.rds'))
centiles <- sort(unique(meta_bins[,'centile_bins'])) # 400 in this version

# library(dplyr)
require(readr)

# Identify the files containing the clusters assigned to each CpG in phase 1 ---
phase1_cluster_files <- list.files(
  path   = input_dir,
  pattern = "^c.*_k[1-9][0-9]?\\.csv$",
  full.names = TRUE
)

# Combine them into one
phase1_clusters <- readr::read_csv(phase1_cluster_files, id = "file", show_col_types = FALSE)

phase1_clusters  <- phase1_clusters |>
  dplyr::mutate(
    file = basename(file),
    centile = sub("^c(.*)_k[0-9]+\\.csv$", "\\1", file),
    cluster = paste(centile, V2, sep = '.')) |>
  dplyr::select(cpg = V1, cluster)

# Save single file
saveRDS(phase1_clusters,
        file.path(input_dir, 'phase1_clusters.rds'))

# Remove individual files
# file.remove(phase1_cluster_files)

# Identify the tuning files, and combine them ----------------------------------
phase1_tune_files <- list.files(
  path   = input_dir,
  pattern = "^c.*_tuning.csv$",
  full.names = TRUE
)

# phase1_tuning <- readr::read_csv(phase1_tune_files, id = "file", show_col_types = FALSE)

# TMP the first files are done slightly differently (much slower, fix later)
phase1_tuning <- phase1_tune_files |>
  purrr::map(~ readr::read_csv(.x, id = "file", show_col_types = FALSE)) |>
  dplyr::bind_rows()

phase1_tuning <- phase1_tuning |>
  dplyr::mutate(
    file = basename(file),
    file = sub("^c(.*)_tuning.csv$", "\\1", file)) |>
  dplyr::rename(centile = file) |> 
  dplyr::select(-`1`) # also 

# save single file
saveRDS(phase1_tuning,
        file.path(input_dir, 'phase1_tuning.rds'))

# Remove individual files
# file.remove(phase1_tune_files)

# Plot it ------------------------------------
tuning_long <- phase1_tuning |>
  # split "centile" like "01.02" into range_centile = 01, median_centile = 02
  tidyr::separate(
    centile,
    into = c("range_centile", "median_centile"),
    sep  = "\\.",
    convert = TRUE       # make them numeric
  ) |>
  tidyr::pivot_longer(
    cols      = `2`:`10`, # numeric columns
    names_to  = "k",
    values_to = "value"
  ) |>
  dplyr::mutate(k = as.integer(k),
       range_centile  = factor(range_centile,  levels = 1:20),
       median_centile = factor(median_centile, levels = 1:20)) |>
  dplyr::group_by(range_centile, median_centile) |>
  dplyr::mutate(
    optimal_k = value == max(value, na.rm = TRUE)) |>
  dplyr::ungroup()


library(ggplot2)

ggplot(
  tuning_long,
  aes(x = k, y = value,
      group = interaction(range_centile, median_centile),
      color = range_centile,
      alpha = median_centile)
) +
  geom_line() +
  geom_point(data = ~ subset(.x, optimal_k),
             size  = 3.5, 
             shape = 15 ) +
  scale_color_viridis_d(name = "Range centile") +
  scale_alpha_discrete(name = "Median centile", range = c(0.1, 1)) +
  guides(
    color = guide_legend(ncol = 2),   # 2-column color legend
    alpha = guide_legend(ncol = 2)    # 2-column alpha legend
  )

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
# file.remove(phase1_centroid_files)