
setwd("~/GENR3/Methylation")

library(maotai)
library(mcclust)

# list.files("./GENR_EPICMETH_Norm")

# list.files("./GENR_450KMETH_Norm_Release3")

message("Loading files...")
load("./GENR_450KMETH_Norm_Release3/GENR_450KMETH_Release3_Betas_ALL_9y_20190813.RData")

# str(x)

# TODO: remove X, Y chromosomes and control probes (annot files )
# TODO: handle missing values 
# TODO: parallelize within 
# TODO: run centiles in parallel with slurm tasks 

# ==============================================================================
# Simple distribution exploration 
# Print distributions 

# pdf('~/DNAm_distrib.pdf')
# par(mfrow=c(2,2))
# 
# lapply(1:10000, function(i){
#   
#   hist(x[i,], breaks=60, xlim=c(0,1), 
#        xlab=row.names(x)[i], 
#        main=row.names(x)[i], border='blue')
#   
#   # fitdistrplus::fitdist(x[1,], "norm")
#   tryCatch(
#     fitdistrplus::descdist(x[i,], discrete = FALSE),
#     error=function(e) NULL
#   )
# })
# 
# 
# dev.off()

# ==============================================================================
# Explore data

# library(moments)
# 
# # Range ------------------------------------------------------------------------
# vrange <- apply(x, 1, function(v) diff(range(v, na.rm=TRUE)))
# summary(vrange)
# hist(vrange, breaks=100, xlim=c(0,1))
# 
# vmeans <- apply(x, 1, mean, na.rm=TRUE)
# hist(vmeans, breaks=100, xlim=c(0,1))
# 
# # Range against mean value 
# plot(vrange, vmeans, cex=1, pch=16, col=scales::alpha("aquamarine4",0.1),
#      xlab = "Variable range", ylab = "Variable mean")
# 
# # Skewness (= asymmetry) -------------------------------------------------------
# # Negative skew = tail is on the left side (towards negative values).
# # Positive skew = tail is on the right side (towards positive values).
# 
# skew <- function(v) {
#   n <- length(v)
#   vmean <- mean(v, na.rm=TRUE)
#   vsd <- sd(v, na.rm=TRUE)
#   
#   sk <- (n * sum((v - vmean)^3)) / ((n - 1) * (n - 2) * vsd^3)
#   return(sk)
# }
# 
# # vskew <- apply(sub, 1, skew)
# vskew <- apply(x, 1, moments::skewness, na.rm=TRUE)
# summary(vskew)
# 
# # what percentage is not symmetrical 
# # 1 - sum(vskew>-1 & vskew <1)/length(vskew)
# 
# hist(vskew, breaks=100, xlim=c(-15,15))
# abline(v = c(-1, -0.5, 0.5, 1), col = c('pink','blue','blue','pink'), 
#        lwd = 1.7)
# 
# # Skewness against mean value 
# plot(vskew, vrange, cex=1, pch=16, col=scales::alpha("aquamarine4",0.1),
#      xlab = "Variable skewness", ylab = "Variable range")
# 
# # Kurtosis ( = heavy- or light-tailed relative to a normal distribution) -------
# # Relative to a normal distribution (kurtosis = 3)
# # < 3 = playkurtic (fewer and less extreme outliers)
# # > 3 = leptokurtic (more outliers)
# 
# kurt <- function(v) {
# 
#   vmean <- mean(v, na.rm=TRUE)
#   vsd <- sd(v, na.rm=TRUE)
#   
#   kt <- sum((v - vmean)^4) / (99 * vsd^4) - 3
#   return(kt)
# }
# 
# # vkurt <- apply(sub, 1, kurt)
# vkurt <- apply(x, 1, moments::kurtosis, na.rm = TRUE)
# summary(vkurt)
# hist(vkurt, breaks=100, xlim=c(0,280))
# abline(v = c(3), col = 'blue', lwd = 1.5)
# 
# # what percentage is leptokurtic 
# # sum(vkurt>3.5)/length(vkurt)
# # sum(vkurt<2.5)/length(vkurt)
# 
# metad <- data.frame(vrange, vskew, vkurt, vmeans)
# row.names(metad) <- row.names(x)
# 
# write.csv(metad, "distrib_metad_450K_9y.csv")
# # clean
# rm(vrange, vskew, vkurt, vmeans)

# plot mean value against the range 

# ==============================================================================
# EP-means clustering
# ==============================================================================

# transform data to list of vectors
# does not handle missing data (TMP)
message("Transforming...")
df_compl <- as.data.frame(t(x[rowSums(is.na(x))==0, ]))
rm(x)

# metadata
metad <- read.csv("~/MPSR/distrib_metad_450K_9y.csv", row.names = 1)
clean_metad <- metad[names(df_compl), ]

# subset large range (above median)
# df_wider <- df[, metad$vrange > 0.19 ]

# does not handle missing data (TMP)
# df_compl <- df[, colSums(is.na(df))== 0]

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
  subsamples <- lapply(1:n_subsamples, function(i) sample(df, length(df), replace = TRUE))
  
  # Run EP-means on each subsample for a range of k values
  clustering_results <- lapply(k_values, function(k) {
    lapply(subsamples, function(subsample) maotai::epmeans(subsample, k=k)$cluster)
  })
  
  # Calculate VI metric for each pair of clusterings
  vi_values <- sapply(clustering_results, function(clusters_k) {
    utils::combn(length(clusters_k), 2, function(index_pair) {
      mcclust::vi.dist(clusters_k[[index_pair[1]]], clusters_k[[index_pair[2]]])
    })
  })
  
  # Assess Stability Across Different k Values:
  # compute average VI value across all pairs, for each k
  avg_vi_values <- sapply(vi_values, mean)
  
  #  k with the lowest average VI indicates the most stable clustering
  optimal_k <- k_values[which.min(avg_vi_values)]
  
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
  return(final_clusters$cluster)

}

# Split into centiles of variable ranges
centiles <- quantile(clean_metad$vrange, seq(0, 1, 0.05)) 

# TODO: run centiles in parallel with slurm tasks 

lapply(c(1:(length(centiles)-1)), function(c) {
  
  df_compl_ci <- df_compl[, (clean_metad$vrange > centiles[c]) & 
                            (clean_metad$vrange < centiles[c+1])]
  message("Data between ", 
          names(centiles)[c], " and ", names(centiles)[c+1], 
          "centiles of variable ranges, aka between ", 
          round(centiles[c],2), " and ", round(centiles[c+1],2), 
          ". (n=",ncol(df_compl_ci),")")
  
  k_ci <- epi_epmeas(df_compl_ci, k_values = 2:10, 
                     output_folder="~/MPSR/clustering",
                     title=paste(substr(
                       names(centiles)[c], 1, nchar(names(centiles)[c])-1), sep='_'))
  
})

# ==============================================================================

# median split - not enough ----------------------------------------------------
# df_compl_narrow <- df_compl[, clean_metad$vrange < 0.17] # median
# k_narrow <- epi_epmeas(df_compl_narrow, title="_narrow")
# 
# message("Wide range...")
# df_compl_wide <- df_compl[, clean_metad$vrange > 0.17]
# k_wide <- epi_epmeas(df_compl_wide, title="_wide")
# 
# message("Everything...")
# k_all <- epi_epmeas(df_compl)

# ==============================================================================
# END
