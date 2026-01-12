use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

library(ggplot2)

require(purrr)
require(tidyr)
require(dplyr)
require(Polychrome)

# ------------------------------------------------------------------------------
# Plotting helpers
# ------------------------------------------------------------------------------

colors <- Polychrome::kelly.colors(21)[-1] # avoid white

# colors <- c("red", "blue","gold", "darkgreen", "pink", 
#             "lightblue", "purple", "lightgreen","brown", "orange",
#             "grey","darkblue","yellow", "cyan", "magenta")

tuning_plot <- function(tuning_df) {
  
  # Transform for ggplot to use 
  tuning_long <- tuning_df |>
    # split "centile" like "01.02" into range_centile = 01, median_centile = 02
    tidyr::separate(centile,
                    into = c("range_centile", "median_centile"), sep  = "\\.",
                    # make them numeric
                    convert = TRUE) |>
    tidyr::pivot_longer(cols = `2`:`10`, # numeric columns indicate evaluated k
                        names_to = "k",
                        values_to = "value") |>
    dplyr::mutate(k = as.integer(k),
                  range_centile  = factor(range_centile,  levels = 1:20),
                  median_centile = factor(median_centile, levels = 1:20),
                  k_selected = k == final_k)
  
  range_colors <- setNames(colors, unique(tuning_long$range_centile))
  
  p <- ggplot(tuning_long,
              aes(x = k, y = value, 
                  group = interaction(range_centile, median_centile),
                  color = range_centile)) +
    geom_line(linewidth = 0.3, alpha = 0.5) +
    geom_point(data = ~ subset(.x, k_selected),
               size  = 2.5, aes(shape = median_centile)) +
    scale_color_manual(name = "Range centile", values = range_colors) +
    scale_shape_manual(name = "Median centile", values = 0:20) +
    guides(color = guide_legend(ncol = 2),   # 2-column color legend
           shape = guide_legend(ncol = 2))
  
  return(p)
}

get_empirical_pdf <- function(ecdf_list, x_range = c(0,1), n_samples = 200, type = 'phase1') {
  
  # Check effective range
  obs_x_range <- range(sapply(ecdf_list, function(f) quantile(f, c(0.01, 0.99))))
  message("Observed x range:")
  print(obs_x_range)
  
  # get the empirical PDF back 
  pdf_data <- purrr::map_dfr(seq_along(ecdf_list), ~ {
    
    # Evaluate ECDF at fine grid
    x_seq <- seq(x_range[1], x_range[2], length.out = n_samples)
    cdf_vals <- ecdf_list[[.x]](x_seq)
    
    # Numerical derivative: PDF ≈ ΔCDF/Δx
    pdf_vals <- diff(cdf_vals) / diff(x_seq)  # vector of length n_samples
    x_pdf <- x_seq[-length(x_seq)] + diff(x_seq)[1]/2  # centers
    
    if (type == 'phase1') {
      cluster_comp <- strsplit(names(ecdf_list)[.x], "\\.")[[1]] |> as.numeric()
      range_group <- sprintf("Range: %02d", cluster_comp[1])
      median_group <- sprintf("Median: %02d", cluster_comp[2])
      cluster <- cluster_comp[3]
      
      return(data.frame(x = x_pdf, y = pdf_vals, cluster = cluster
                        range_group = range_group,
                        median_group = median_group))
    } else {
      
      return(data.frame(x = x_pdf, y = pdf_vals, cluster = .x))
    }
  })
  
  return(pdf_data)
}


centroid_plot <- function(pdf_data, output_file = 'centroids_by_range_centile.pdf') {
  
  # Plot on PDF 
  pdf(output_file, width = 10, height = 8)
  
  if ('range_group' %in% names(pdf_data)) {
    
    # median_colors <- setNames(colors, unique(pdf_data$median_group))
    
    purrr::walk(unique(pdf_data$range_group), ~ {
      
      subset_data <- dplyr::filter(pdf_data, range_group == .x)
      
      p <- ggplot(subset_data, aes(x = x, y = y, 
                                   color = factor(median_group), 
                                   linetype = factor(cluster))) +
        geom_line(size = 0.8, alpha = 0.8) +
        scale_color_viridis_d(name = "Median centile") +
        scale_linetype_discrete(name = "Cluster") +
        labs(title = .x, x = "Methylation", y = "Density") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      print(p)
    })
    
  } else {
    cluster_levels <- unique(pdf_data$cluster)
    cluster_colors <- setNames(colors[1:length(cluster_levels)], cluster_levels)
    
    p <- ggplot(pdf_data, aes(x = x, y = y,
                              color = factor(cluster))) +
      geom_line(size = 0.5, alpha = 0.7) +
      scale_color_manual(name = "Cluster", values = cluster_colors) +
      labs(title = 'Final clusters', x = "Methylation", y = "Density") +
      theme_minimal() +
      theme(legend.position = "bottom") + 
      guides(color = guide_legend(ncol = 10))
    
    print(p)
  }
  
  dev.off()
  
  return(invisible(NULL))
  
}

clusters_plot <- function(cpg_data, cluster_data, output_file = "density_by_cluster.pdf", 
                          fixed_x_range = FALSE) {
  
  # TODO: Pre-compute global y-limit for consistent scaling across page
  # y_max <- 200  # max density
  
  pdf(output_file, width = 12, height = 10)
  
  for (range_g in unique(cluster_data$range_centile)) {
    
    par(mfrow = c(4, 5), mar = c(1, 1, 1, 1), oma = c(0, 0, 2, 0), 
        yaxt  = "n")  # suppress y axis ticks & labels
    
    for (median_g in unique(cluster_data$median_centile)) {
      
      cluster_subset <- cluster_data |>
        dplyr::filter(range_centile == range_g, median_centile == median_g)
      
      cpg_idx <- match(cluster_subset$cpg, rownames(cpg_data))
      cpg_idx <- cpg_idx[!is.na(cpg_idx)] # should not be any missing matched but just in case
      
      # Subset valid data only + get colors
      cpg_data_subset <- cpg_data[cpg_idx, ]
      cluster_color <- colors[as.integer(cluster_subset$cluster[cpg_idx])]
      
      # Pre-compute ALL densities first for limits
      all_densities <- apply(cpg_data_subset, 1, density, bw = "SJ", n = 500)
      
      y_range <- c(0, max(sapply(all_densities, function(d) max(d$y)), na.rm = TRUE))
      
      if (fixed_x_range) {
        x_range <- c(0, 1)
      } else {
        x_range <- range(sapply(all_densities, function(d) range(d$x)), na.rm = TRUE)
      }
      
      # Plot frame with axes
      plot(0, type = "n", xlim = x_range, ylim = y_range, main = paste('Median centile:', median_g), 
           cex.main = 0.8, xlab = "Methylation", ylab = "Density")
      
      # Draw densities
      for (i in seq_len(nrow(cpg_data_subset))) {
        lines(all_densities[[i]], 
              col = scales::alpha(cluster_color[i], 0.3), lwd = 0.3)
      }
    }
    
    mtext(paste("Range centile:", range_g), outer = TRUE, cex = 1.2, line = 0)
  }
  
  dev.off()
  
  return(ividible(NULL))
}
