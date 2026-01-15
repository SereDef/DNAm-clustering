use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

library(ggplot2)
library(patchwork)

# require(purrr)
# require(tidyr)
# require(dplyr)

# ------------------------------------------------------------------------------
# Plotting helpers
# ------------------------------------------------------------------------------

# remove Polychrome dependency while mantaining optimal perceptual differences
# colors <- Polychrome::kelly.colors(21)[-1] # avoid white

colors <- setNames(
  c("#3D4C8B",  # deep blue (replace one of the close purples)
    "#F3C300", # vivid yellow
    "#2E8B57", # sea green (swap for one of the very similar yellows)
    "#875692", # strong purple
    "#BE0032", # strong red
    "#A1CAF1", # light blue
    "#C2B280", # sand
    "#F38400", # strong orange
    "#848482", # medium grey
    "#008856", # green
    "#E68FAC", # pink
    "#0067A5", # blue
    "#F6A603", # yellow-orange
    "#B3446C", # purplish red
    "#DCD300", # yellow-green
    "#604E97", # violet
    "#F99379", # light orange
    "#654522", # brown
    "#E25822", # reddish orange
    "#222222"  # black-ish, very dark grey
  ),
  1:20
)

# ==============================================================================
# --- Plot k tuning process ---
# ==============================================================================

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
  
  # range_colors <- setNames(colors, unique(tuning_long$range_centile))
  
  p <- ggplot(tuning_long,
              aes(x = k, y = value, 
                  group = interaction(range_centile, median_centile),
                  color = range_centile)) +
    geom_line(linewidth = 0.3, alpha = 0.5) +
    geom_point(data = ~ subset(.x, k_selected),
               size  = 2.5, aes(shape = median_centile)) +
    scale_color_manual(name = "Range centile", values = colors) +
    scale_shape_manual(name = "Median centile", values = 0:20) +
    guides(color = guide_legend(ncol = 2),   # 2-column color legend
           shape = guide_legend(ncol = 2))
  
  return(p)
}

# ==============================================================================
# --- Plot cluster centroids ---
# ==============================================================================

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
      range_group <- cluster_comp[1]
      median_group <- cluster_comp[2]
      cluster <- cluster_comp[3]
      
      return(data.frame(x = x_pdf, y = pdf_vals, cluster = cluster,
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
    
    purrr::walk(unique(pdf_data$range_group), ~ {
      
      subset_data <- dplyr::filter(pdf_data, range_group == .x)
      
      # pre-compute within-range_group max y to scale offsets
      offset_step <- max(subset_data$y) * 0.3
      
      subset_data$y_offset <- (subset_data$median_group - 1L) * offset_step
      subset_data$new_y <- subset_data$y + subset_data$y_offset
      
      p <- ggplot(subset_data, aes(x = x, y = new_y, 
                                   color = factor(cluster), 
                                   group = interaction(cluster, median_group))) +
        geom_line(linewidth = 0.5, alpha = 0.8) +
        scale_color_manual(name = "Cluster", values = colors) +
        scale_y_continuous(name   = "Density (offset by median group)",
                           breaks = subset_data$y_offset,
                           labels = subset_data$median_group) +
        labs(title = .x, x = "Methylation") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      print(p)
    })
    
  } else {
     
    p <- ggplot(pdf_data, aes(x = x, y = y,
                              color = factor(cluster))) +
      geom_line(size = 0.5, alpha = 0.7) +
      scale_color_manual(name = "Cluster", values = colors) +
      labs(title = 'Final clusters', x = "Methylation", y = "Density") +
      theme_minimal() +
      theme(legend.position = "bottom") + 
      guides(color = guide_legend(ncol = 10))
    
    print(p)
  }
  
  dev.off()
  
  return(invisible(NULL))
  
}

# ==============================================================================
# --- Plot individual CpG densities ---
# ==============================================================================

cpg_subset_densities <- function(cluster_subset, cpg_data, fixed_x_range = TRUE, 
                                 title, color) {
  
  cpg_idx <- match(cluster_subset$cpg, rownames(cpg_data))
  cpg_idx <- cpg_idx[!is.na(cpg_idx)] # should not be any missing matches but just in case
  
  # Subset valid data only + get colors
  cpg_data_subset <- cpg_data[cpg_idx, ]
  
  # Pre-compute ALL densities first for limits (rows = cpgs)
  all_densities <- apply(cpg_data_subset, 1, density, bw = "SJ", n = 500)
  
  y_range <- c(0, max(sapply(all_densities, function(d) max(d$y)), na.rm = TRUE))
  
  if (fixed_x_range) {
    x_range <- c(0, 1)
  } else {
    x_range <- range(sapply(all_densities, function(d) range(d$x)), na.rm = TRUE)
  }
  
  # Plot frame with axes
  plot(0, type = "n", xlim = x_range, ylim = y_range, 
       main = title, 
       cex.main = 0.8, xlab = "Methylation", ylab = "Density")
  
  if (length(color) == 1L) {
    cluster_colors <- rep(colors[color], nrow(cpg_data_subset))
  } else {
    cluster_colors <- colors[as.integer(cluster_subset$cluster[cpg_idx])]
  }
  
  # Draw densities
  for (i in seq_len(nrow(cpg_data_subset))) {
    lines(all_densities[[i]], col = scales::alpha(cluster_colors[i], 0.3), lwd = 0.3)
  }
  
}

clusters_plot <- function(cpg_data, cluster_data, cluster_var = 'p1_cluster', 
                          output_file = "density_by_cluster.pdf", 
                          fixed_x_range = FALSE) {
  
  # TODO: Pre-compute global y-limit for consistent scaling across page
  # y_max <- 200  # max density
  
  pdf(output_file, width = 12, height = 10)
  
  if (cluster_var == 'p2_cluster') {
    n_clusters <- max(as.numeric(cluster_data$p2_cluster))
    for (cluster in 1:n_clusters) {
      
      cluster_subset <- cluster_data |>
        dplyr::filter(p2_cluster == cluster)
      
      cpg_subset_densities(cluster_subset, cpg_data, fixed_x_range = TRUE,
                           title = paste('Cluster:', cluster), 
                           color = cluster)
    }
  } else if (cluster_var == 'p1_cluster') {
    
    for (range_g in unique(cluster_data$range_centile)) {
      
      par(mfrow = c(4, 5), mar = c(1, 1, 1, 1), oma = c(0, 0, 2, 0), 
          yaxt  = "n")  # suppress y axis ticks & labels
      
      for (median_g in unique(cluster_data$median_centile)) {
        
        cluster_subset <- cluster_data |>
          dplyr::filter(range_centile == range_g, median_centile == median_g)
        
        n_clusters <- max(as.integer(cluster_subset$cluster))
        
        cpg_subset_densities(cluster_subset, cpg_data, fixed_x_range = fixed_x_range,
                             title = paste('Median centile:', median_g), 
                             color = 1:n_clusters)
      }
      mtext(paste("Range centile:", range_g), outer = TRUE, cex = 1.2, line = 0)
    }
  }
  
  dev.off()
  
  return(invisible(NULL))
}

# ==============================================================================
# --- Plot cluster representation in EWAS results ---
# ==============================================================================

manhattan <- function(dset, ewas, cluster_var = 'p2_cluster', 
                      thresh_gnmwide = 1e-7,
                      thresh_suggest = 1e-5, 
                      annotate_top = 50) {
  
  # Pre-process cluster information ============================================
  n_clusters <- length(unique(dset[, cluster_var]))
  
  # set CpGs not assigned to cluster as separate category
  dset[, cluster_var][is.na(dset[, cluster_var])] <- 'NA'
  
  cluster_colors <- colors
  cluster_colors[['NA']] <- '#f7f7f7'
  
  # Pre-process p-value information ============================================
  
  dset$log_pvalue <- -log10(dset[, paste('pvalue', ewas, sep="_")])
  
  ewas_name <- gsub('\\.', ' ', ewas)
  
  # Pre-process chromosome positions ===========================================
  
  # order by genomic position
  dset <- dset[order(dset$chr, dset$pos), ]
  
  chr_lengths <- tapply(dset$pos, dset$chr, max)
  chr_offsets <- c(0, cumsum(as.numeric(head(chr_lengths, -1))))
  names(chr_offsets) <- names(chr_lengths)
  
  dset$cum_pos <- dset$pos + chr_offsets[as.character(dset$chr)]
  
  # sig_cpgs <- subset(dset, log_pvalue >= -log10(thresh_gnmwide))
  top_cpgs <- dset[order(-dset$log_pvalue), ][1:annotate_top, ]
  
  ggplot(dset,
         aes(x = cum_pos,
             y = log_pvalue,
             color = factor(.data[[cluster_var]]))) +
    geom_point(alpha = 0.8, size = 0.5) +
    # Color clusters 
    scale_color_manual(values = cluster_colors, name = "Cluster") +
    # Position chr labels
    scale_x_continuous(name = "Chromosome",
                       breaks = tapply(dset$cum_pos, dset$chr, mean),
                       labels = names(chr_lengths)) +
    # Add significance thresholds
    geom_hline(yintercept = -log10(thresh_gnmwide), color = "grey") +
    geom_hline(yintercept = -log10(thresh_suggest), color = "grey", linetype = "dashed") +
    # Annotate clusters
    ggrepel::geom_text_repel(data = top_cpgs, 
                             aes(label = .data[[cluster_var]]),
                             size = 2.5, max.overlaps = Inf, box.padding = 0.3) +
    # Labels, title, legend
    labs(y = expression(-log[10](p)), title = paste(ewas_name, "EWAS")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) + 
    guides(color = guide_legend(ncol = 11))
  
}

cluster_representation <- function(dset, ewas, cluster_var = 'p2_cluster', 
                                   threshold = 1e-7) {
  
  # Pre-process p-value information ============================================
  ewas_name <- gsub('\\.', ' ', ewas)
  dset$signif <- dset[, paste('pvalue', ewas, sep="_")] < threshold
  # Ensure table get all levels always
  dset$signif <- factor(dset$signif, levels = c(FALSE, TRUE))
  
  # Compute log(OR) for comparison =============================================
  # Cross tab for comparison: rows = clusters, columns = FALSE/TRUE (significant)
  crosstab <- as.data.frame.matrix(table(dset[, cluster_var], dset$signif))
  names(crosstab) <- c("n_nonsig", "n_sig")
  crosstab$cluster <- rownames(crosstab)
  
  # total per cluster
  crosstab$n_total <- crosstab$n_sig + crosstab$n_nonsig
  
  # overall rates
  overall_sig_rate <- sum(crosstab$n_sig) / sum(crosstab$n_total)
  overall_nonsig_rate <- 1 - overall_sig_rate
  
  # avoid zeros by adding a small continuity correction (if needed)
  eps <- 0.001
  
  # Compute log(OR) [better properties: symmetric around 0]
  # + ⇒ cluster over‑represented among significant CpGs
  # - ⇒ under‑represented
  
  crosstab$lor <- with(crosstab, {
    # cluster odds of being significant
    odds_cluster  <- (n_sig + eps) / (n_nonsig + eps)
    # overall odds
    total_sig <- sum(n_sig)
    total_nonsig  <- sum(n_nonsig)
    odds_overall  <- (total_sig + eps) / (total_nonsig + eps)
    log(odds_cluster / odds_overall)
  })
  
  # panel 1: log-odds representation
  cluster_lor <- ggplot(crosstab,
                        aes(x = lor, y = reorder(cluster, lor), 
                            fill = lor > 0)) +
    geom_col() +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "#8B0000"),
                      labels = c("Under", "Over"),
                      name = "Represented") +
    labs(x = "Log-odds of significant (vs overall)", y = "Cluster",
         title = paste(ewas, "EWAS")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # panel 2: overall counts per cluster (grey bars)
  cluster_counts <- ggplot(crosstab,
                           aes(x = n_total, y = reorder(cluster, lor))) +
    geom_col(fill = "#EDEDED") +
    labs(x = "Total CpGs per cluster", y = "Cluster") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.y = element_blank())
  
  # patchwork syntax
  (cluster_lor | cluster_counts) + plot_layout(widths = c(2, 1))
 
}

ewas_plot <- function(dset, ewases, cluster_var = 'p2_cluster', 
                      thresh_gnmwide = 1e-7, output_file = 'EWAS_cluster.pdf',
                      ...) {
  
  pdf(output_file, width = 12, height = 16)
  
  for (ewas in ewases) {
    p_manhattan <- manhattan(dset, ewas, cluster_var, 
                             thresh_gnmwide = thresh_gnmwide, ...)
    
    p_represent <- cluster_representation(dset, ewas, cluster_var, 
                                          threshold = thresh_gnmwide)
    
    print((p_manhattan / p_represent) + plot_layout(heights = c(1.5, 2)))
    
  }
  
  dev.off()
  
  return(invisible(NULL))
}