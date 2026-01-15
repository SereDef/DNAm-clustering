use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

input_dir <- "~/MPSR/metadata"

output_file <- "../metadata/CpG_metadata.rds"

meta_desc <- file.path(input_dir, 'metadata_mega.desc')
meta_desc_ComBat <- file.path(input_dir, 'metadata_mega_ComBat.desc')

meta_filter <- bigmemory::attach.big.matrix(meta_desc, lockfile = TRUE)
meta_ComBat <- bigmemory::attach.big.matrix(meta_desc_ComBat, lockfile = TRUE)

# Which distribution moments to use to determine groups for stage 1 ------------
based_on_sample <- 'mega'
based_on_moment <- c('median', 'range')

based_on <- paste(based_on_moment, based_on_sample, sep = '_')

# Extract only the moments needed in computation
moment_sub <- meta_ComBat[, based_on]

message('N features: ', nrow(moment_sub))
print(summary(moment_sub))

# Clean up range 
# Batch correction results in ranges above 1 for a small number of CpGs. 
# I am not sure why, but I will just remove them:
if ('range' %in% based_on_moment) {
  range_name <- grep('range', based_on, value = TRUE)
  range_var <- moment_sub[, range_name]
  
  cat(sum(range_var > 1.00), ' CpGs have ranges above 1. Removing them...')
  
  moment_sub <- moment_sub[range_var <= 1.00, ]
  moment_sub <- moment_sub[is.finite(moment_sub[, range_name]), ]
}

message('N features: ', nrow(moment_sub))
print(summary(moment_sub))

# Strategy: ====================================================================
# - 5 bins of median centiles
# -- 20 bins of range centiles

moment1 <- grep('median', based_on, value = TRUE)
moment2 <- grep('range', based_on, value = TRUE)

cutoff1 <- 0.20 # 5 bins
cutoff2 <- 0.05 # 20 bins

# Grouping 1
centile_breaks1 <- quantile(moment_sub[, moment1], 
                           probs = seq(0, 1, by = cutoff1), 
                           na.rm = TRUE)

centiles_bins1 <- findInterval(moment_sub[, moment1], vec = centile_breaks1,
                               rightmost.closed = TRUE,  all.inside = TRUE)

centile_bins <- vector(mode='character', nrow(moment_sub))
centile_breaks2 <- vector(mode='list', length(centile_breaks1)-1)

for (c1 in 1:max(centiles_bins1)){

  c1_idx = which(centiles_bins1 == c1)
  
  c1_median <- moment_sub[c1_idx, moment2]
  
  c1_centile_breaks2 <- quantile(c1_median,
                                 probs = seq(0, 1, by = cutoff2),
                                 na.rm = TRUE)
  centile_bins2 <- findInterval(c1_median, vec = c1_centile_breaks2,
                                rightmost.closed = TRUE,  all.inside = TRUE)
  
  centile_bin_names <- paste(sprintf("%02d", c1), 
                             sprintf("%02d", centile_bins2), sep='.')
  
  centile_bins[c1_idx] <- centile_bin_names
  centile_breaks2[[c1]] <- c1_centile_breaks2
}

table(centile_bins)

metadata_binned <- cbind(meta_ComBat[rownames(moment_sub), ], centile_bins)

saveRDS(metadata_binned, file=output_file)

# ==============================================================================
# Plot this 

library(ggplot2)
library(patchwork)

# 1) Tidy data for first-level breaks --------------------------------------
df1 <- tidyr::tibble(bin1 = seq_len(length(centile_breaks1) - 1),
                     xmin = centile_breaks1[-length(centile_breaks1)],
                     xmax = centile_breaks1[-1])

# 2) Tidy data for nested (second-level) breaks ----------------------------
df2 <- purrr::imap_dfr(centile_breaks2, ~ { brks <- .x
                       tidyr::tibble(bin1 = as.integer(.y),
                              bin2 = as.factor(seq_len(length(brks) - 1)),
                              xmin = brks[-length(brks)],
                              xmax = brks[-1]) })

band <- 0.45

# First plot: top-level centiles (one row per bin1) ---------------------
p1 <- ggplot(df1) +
  geom_rect(aes(xmin = xmin, xmax = xmax, 
                ymin = bin1 - band, ymax = bin1 + band),
            fill = "#E9E9E9", color = "#FFFFFF") +
  scale_y_reverse(breaks = df1$bin1, expand = expansion(add = 0.5)) + # labels
  labs(x = "Median", y = "Level-1 bins") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Second plot: nested centiles with color-coded bins --------------------
p2 <- ggplot(df2) +
  geom_rect(aes(xmin = xmin, xmax = xmax, 
                ymin = bin1 - band, ymax = bin1 + band,
                fill = bin2), color = "#FFFFFF") +
  scale_y_reverse(breaks = sort(unique(df2$bin1)), 
                   expand = expansion(add = 0.5)) +
  labs(x = "Range", y = NULL, fill = "Nested bins") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank()) +
  guides(fill = guide_legend(ncol = 2))


plot_title <- sprintf("Binning strategy (%s groups; group size: %s CpGs)", 
                      length(table(centile_bins)),
                      table(centile_bins)[1])

# Combine plots side-by-side --------------------------------------------

(p1 | p2) + plot_layout(widths = c(1, 1.3)) + plot_annotation(title = plot_title) 

