use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

# install.packages("maotai", repos = "http://cran.us.r-project.org")

require(maotai)
require(mcclust)
require(parallel)

# ------------------------------------------------------------------------------
# Helper: run epmeans on a subset of CpGs, then assign ALL CpGs to clusters
# ------------------------------------------------------------------------------

subsample_epmeans <- function(ecdf_list, k, beta, n_cpg) {
  # ecdf_list: list of CpG distribution ECDFs
  # k: number of clusters
  # beta: fraction of CpGs used to fit centroids
  
  stopifnot(beta > 0, beta <= 1)
  
  # 1) sample a subset of CpGs
  n_sub <- max(1L, floor(beta * n_cpg))
  idx_sub <- sample.int(n_cpg, n_sub, replace = TRUE)  # with replacement
  ecdf_sub <- ecdf_list[idx_sub]
  
  # 2) run EP-means on subset (each row is a distribution)
  # - epmeans expects a list of vectors (empirical distributions)
  
  fit_sub <- maotai::epmeans(ecdf_sub, k = k)
  
  # 3) assign ALL CpGs to nearest centroid
  #    Build list of all CpG distributions and re-use the learned centers.
  # NOTE: maotai::epmeans does not expose a "predict", so I find the closest
  #    prediction by minimal Wasserstein/EP distance.
  
  closest_centroid <- function(cpg_ecdf, centers) {
    dists = sapply(centers, function(centroid) {
      # ecdfdist returns a symmetric  matrix of pairwise distances, keep only value
      maotai::ecdfdist(list(cpg_ecdf, centroid), method="Wasserstein")[1,2]
    })
    return(which.min(dists))
  }
  
  centroids <- fit_sub$centers
  
  # Assign a cluster to all cpgs
  fit_all <- sapply(ecdf_list, function(cpg) closest_centroid(cpg, centroids))
  
  return(fit_all)
}

# ------------------------------------------------------------------------------
# Main: stability-based model selection as in Henderson et al. (2015)
# ------------------------------------------------------------------------------

stability_selection <- function(ecdf_list,
                                k_values = 2:10,
                                beta = 0.5,
                                max_iter = 50,
                                seed = 3108,
                                n_cores = 1) {
  # ecdf_list: list of CpG ECDFs
  # k_values: candidate numbers of clusters
  # beta: fraction of CpGs to subsample (0 < beta <= 1)
  # max_iter: number of iterations per k
  # seed: RNG seed for reproducibility
  
  set.seed(seed)
  
  # Precompute:
  # - number of CpGs
  n_cpg <- length(ecdf_list)
  
  # - all pairs of iterations (i.e. max_iter!)
  comb_idx <- as.data.frame(utils::combn(max_iter, 2))
  
  # Mean stability per k value
  S_k <- sapply(k_values, function(k){
    
    message("Evaluating k = ", k)
    
    # 1) Get max_iter clusterings for this k
    clusters_k <- parallel::mclapply(seq_len(max_iter), function(i) {
      subsample_epmeans(ecdf_list, k = k, beta = beta, n_cpg = n_cpg)
    }, mc.cores = n_cores)
    
    # 2) Compute pairwise stability S(Ci, Cj)
    #    VI-based similarity: S = 1 - VI(C1, C2) / (2 * log(k))
    S_vals <- sapply(comb_idx, function(j) {
      vi <- mcclust::vi.dist(clusters_k[[j[1]]], clusters_k[[j[2]]])
      return(1 - vi / (2 * log(k)))
    })
    
    # 3) Average 
    return(mean(S_vals))
  })
  
  names(S_k) <- k_values # For display purposes
  
  return(S_k)
}


fit_epmeas <- function(ecdf_list, k, centile, output_folder) {
  
  # Input check 
  if (!dir.exists(output_folder)) dir.create(output_folder, showWarnings = FALSE,
                                             recursive = TRUE)
  
  message("Fitting ", k, " cluters to the full sample...")
  final_clusters = maotai::epmeans(ecdf_list, k = k) 

  
  clusters <- cbind(names(ecdf_list), final_clusters$cluster)
  
  centroids <- final_clusters$centers
  
  if (centile == 'phase2') {
    names(clusters) <- c('centile','cluster')
    centile_name <- centile 
  } else {
    names(clusters) <- c('cpg','cluster')
    centile_name <- paste0('c', centile) 
  }
  
  # Save output
  saveRDS(clusters, 
          file = file.path(output_folder, paste0(centile_name,'_k',k,'.csv')))

  saveRDS(centroids, 
          file = file.path(output_folder, paste0(centile_name,'_centroids.rds')))
  
  return(final_clusters)
}

# plot_clusters <- function(final_clusters, k, centile, output_folder,
#                           colors = c("red", "pink", "gold", 
#                                      "blue", "lightblue", "darkgreen",
#                                      "purple", "grey","brown",
#                                      "orange","lightgreen","darkblue",
#                                      "yellow", "cyan", "magenta"){
#   message("\nPlotting distributions per cluster ...")
#                             
#   pdf(file.path(output_folder, paste0('c',centile,'_k',k,'.pdf')))
#   
#   # Plotting params
#   plotting_dims <- switch(
#     findInterval(k, c(1, 4, 7, 10, 17), rightmost.closed = TRUE),
#     c(1, k),                     # k < 4
#     c(2, ceiling(k / 2)),        # 4–6
#     c(3, 3),                     # 7–9
#     c(4, 4),                     # 10–16
#     c(4, 5)                      # 17-20
#   )
#   
#   par(mfrow=plotting_dims)
#   
#   message("\nPlotting distributions per cluster ...")
#   for (g in 1:k) {
#     
#     plot.new()
#     plot.window(xlim=c(-0.1, 1.1), ylim = c(0, 70)) #TODO: smart upper limit
#     
#     grp <- df[, which(final_clusters$cluster == g)]
#     
#     group_desc <- paste0("Group ", g, " (n = ", ncol(grp), ")")
#     message(group_desc)
#     
#     # Draw individual CpGs
#     apply(grp, 2, function(x) {
#       rng <- range(x)
#       lines(stats::density(x, bw="SJ", from = rng[1], to = rng[2]), 
#             col=scales::alpha(colors[g], 0.2), lw=0.5)}
#     )
#     # Draw centroids 
#     lines(stats::density(stats::knots(final_clusters$centers[[g]])), 
#           col="black", lw=1.2)
#     
#     # Add ticks
#     axis(1, at = seq(0, 1, 0.1))
#     title(main=group_desc, xlab = "Methylation", ylab="Density")
#   }
#   
#   dev.off()
# }