use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

output_dir <- "~/MPSR/metadata"

meta_desc <- file.path(output_dir, 'metadata_mega.desc')
meta_desc_ComBat <- file.path(output_dir, 'metadata_mega_ComBat.desc')

meta_filter <- bigmemory::attach.big.matrix(meta_desc, lockfile = TRUE)
meta_ComBat <- bigmemory::attach.big.matrix(meta_desc_ComBat, lockfile = TRUE)

moment_list <- unique(sapply(strsplit(colnames(meta_filter), '_'), `[`, 1))
subset_list <- unique(sapply(strsplit(colnames(meta_filter), '_'), `[`, 2))

rgba_color <- function(colors, alpha) {
  
  return(
    sapply(colors, function(col) {
      rgb_val <- col2rgb(col) / 255
      rgb(rgb_val[1], rgb_val[2], rgb_val[3], 
          alpha = alpha)})
    )
  
}

hist_moment <- function(moment, meta_x = meta_filter, ylim=c(0, 10), xlim = c(0, 1),
                        colors = c('orange','navyblue','darkred','seagreen'), 
                        show_legend = FALSE,
                        alpha = 0.1) {
  
  # Convert color names/hex to RGB with alpha
  rgba_colors <- rgba_color(colors, alpha)
  
  # Plot the first histogram
  hist(meta_x[, paste0(moment, '_mega')], breaks = 1000,
       col = rgba_colors[1], border = rgba_colors[1],
       freq = FALSE, xlim = xlim, ylim = ylim, xlab = '', 
       main = moment)
  
  # Overlay remaining histograms
  for (subset in 2:length(subset_list)) {
    hist(meta_x[, paste(moment, subset_list[subset], sep='_')], breaks = 1000,
         col = rgba_colors[subset], border = rgba_colors[subset],
         freq = FALSE, add = TRUE)
  }
  # Add legend
  if (show_legend) {
    op <- par(xpd = NA)  # allow drawing outside plot region
    legend("topright", legend = subset_list, 
                          fill = rgba_color(colors, 0.5), ncol = 4)
    par(op)  # restore previous settings
  }
}

par(mfrow = c(2, length(moment_list)))

# Plot the histograms for each moment
hist_moment('mean', ylim = c(0, 9))
hist_moment('median', ylim = c(0, 9))
hist_moment('range', ylim = c(0, 6), xlim = c(-0.001, 2))
hist_moment('variance', ylim = c(0, 20), xlim=c(0, 0.15))
hist_moment('skewness', ylim = c(0, 0.4), xlim=c(-50, 50))
hist_moment('kurtosis', ylim = c(0, 0.3), xlim=c(1, 200), show_legend = TRUE)

hist_moment('mean', meta_x = meta_ComBat, ylim = c(0, 9))
hist_moment('median', meta_x = meta_ComBat, ylim = c(0, 9))
hist_moment('range', meta_x = meta_ComBat, ylim = c(0, 6), xlim = c(-0.001, 2))
hist_moment('variance', meta_x = meta_ComBat, ylim = c(0, 20), xlim=c(0, 0.15))
hist_moment('skewness', meta_x = meta_ComBat, ylim = c(0, 0.55), xlim=c(-50, 50))
hist_moment('kurtosis', meta_x = meta_ComBat, ylim = c(0, 0.45), xlim=c(1, 200))



for (m in colnames(meta_filter)) {
  message(m) 
  print(summary(meta_filter[, m]))
  print(summary(meta_ComBat[, m]))
  cat('\n')
}

colSums(meta_filter[, grep('range', colnames(meta_filter))] > 1)

colSums(meta_ComBat[, grep('range', colnames(meta_filter))] > 1)
