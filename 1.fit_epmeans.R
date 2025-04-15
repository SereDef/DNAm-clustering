
timepoint <- "birth"
array <- "450K"

output_folder <- "~/MPSR/"

install.packages("maotai", repos = "http://cran.us.r-project.org")

library(maotai)
library(mcclust)

# TODO: remove X, Y chromosomes and control probes (annot files )
# TODO: handle missing values 

# ==============================================================================
# EP-means clustering
# ==============================================================================

message("Loading files...")
clean_metad <- read.csv(file.path(output_folder, paste0("distrib_metad_",array,"_",timepoint,".csv")),
                  row.names = 1)

load(file.path(output_folder, paste0("data_filtered_",array,"_",timepoint,".RData")))

epi_epmeas <- function(df, k_values = 2:10, 
                       title="", 
                       output_folder=getwd(),
                       colors = c("red", "pink", "gold", 
                                 "blue", "lightblue", "darkgreen",
                                 "purple", "grey","brown","orange")
                       ) {
  # TODO: check df is complete (TMP)
  
  message("STEP 1: Estimating k...")
  
  # Partition into random sub-samples w/ replacement (for quicker estimates)
  set.seed(3108)
  n_subsamples <- 100
  size_subsample <-  length(df)/10 
  
  subsamples <- lapply(1:n_subsamples, function(i) sample(df, size_subsample, replace = TRUE))
  
  # Run EP-means on each subsample for a range of k values
  clustering_results <- lapply(k_values, function(k) {
    message("Testing k == ", k)
    lapply(subsamples, function(subsample) maotai::epmeans(subsample, k=k)$cluster)
  })
  
  message("Calculate VI metric...")
  # Calculate VI metric for each pair of clusterings
  vi_values <- sapply(clustering_results, function(clusters_k) {
    utils::combn(length(clusters_k), 2, function(index_pair) {
      mcclust::vi.dist(clusters_k[[index_pair[1]]], clusters_k[[index_pair[2]]])
    })
  })
  
  # Assess Stability Across Different k Values:
  # compute average VI value across all pairs, for each k
  avg_vi_values <- colMeans(vi_values) # sapply(vi_values, mean)
  
  message(avg_vi_values)
  
  #  k with the lowest average VI indicates the most stable clustering
  optimal_k <- k_values[which.max(round(avg_vi_values))]
  
  message("STEP 2: Estimating, ", optimal_k, "cluters...")
  final_clusters = maotai::epmeans(df, k=optimal_k) 
  
  message("\nPlotting distributions per cluster ...")
  pdf(file.path(output_folder, paste0('clust_cent_',title,'_k_',optimal_k,'.pdf')))
  
  # Plotting param (fix later)
  if (optimal_k < 4) { plotting_dims <- c(1, optimal_k)
  } else if (optimal_k < 7) { plotting_dims <- c(2, ceiling(optimal_k/2))
  } else if (optimal_k < 10) { plotting_dims <- c(3, 3)
  } else { plotting_dims <- c(4, 4) }
  
  message("\nPlotting distributions per cluster ...")
  for (g in 1:optimal_k) {
    
    plot.new()
    plot.window(xlim=c(-0.1, 1.1), ylim = c(0, 70)) #TODO: smart upper limit
    
    grp <- df[, which(final_clusters$cluster == g)]
    
    group_desc <- paste0("Group ", g, " (n = ", ncol(grp), ")")
    message(group_desc)
    
    # Draw individual CpGs
    apply(grp, 2, function(x) {
      rng <- range(x)
      lines(stats::density(x, bw="SJ", from = rng[1], to = rng[2]), 
            col=scales::alpha(colors[g], 0.2), lw=0.5)}
    )
    # Draw centroids 
    lines(stats::density(stats::knots(final_clusters$centers[[g]])), 
          col="black", lw=1.2)
    
    # Add ticks
    axis(1, at = seq(0, 1, 0.1))
    title(main=group_desc, xlab = "Methylation", ylab="Density")
  }
  
  dev.off()
  
  # Save output
  write.csv(final_clusters$cluster, 
            file.path(output_folder, paste0('clust_cent_',title,'_k_',optimal_k,'.csv')))
  write.csv(avg_vi_values, 
            file.path(output_folder, paste0('clust_cent_',title,'_k_tuning.csv')))
  
  # final_clusters$avg_vi_values <- avg_vi_values
  return(final_clusters$cluster)

}

# Split into centiles of variable ranges
centiles <- quantile(clean_metad$vrange, seq(0, 1, 0.05)) 

# Run centiles in parallel with slurm tasks 
library(parallel)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# cl <- makeCluster(n_cores) # Create a cluster
# Export needed variables/functions
# clusterExport(cl, c("df_compl", "clean_metad", "centiles", "epi_epmeas"))

results <- parallel::mclapply(
  1:(length(centiles)-1), 
  function(centile) {
  
  message('Subsetting centile...')
  df_compl_ci <- df_compl[, (clean_metad$vrange > centiles[centile]) & 
                            (clean_metad$vrange < centiles[centile+1])]
  message("Data between ", 
          names(centiles)[centile], " and ", names(centiles)[centile+1], 
          " centiles of variable ranges, aka between ", 
          round(centiles[centile],2), " and ", round(centiles[centile+1],2), 
          ". (n=",ncol(df_compl_ci),")")
  
  k_ci <- epi_epmeas(df_compl_ci, k_values = 2:10, 
                     output_folder="~/MPSR/clustering_results",
                     title=paste(substr(
                       names(centiles)[centile], 1, nchar(names(centiles)[centile])-1), sep='_'))
  }, 
  mc.cores = n_cores)

# stopCluster(cl)

