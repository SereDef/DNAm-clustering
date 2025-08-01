# Install if needed
# install.packages(c("bigmemory"))

# =================================== CONFIG ===================================
dataset <- "GENR"
timepoint <- "birth"
normalization <- "funcnorm"

batch_corrected <- TRUE

use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

input_dir <- '~/MPSR/metadata'
output_dir <- "~/MPSR/metadata"

# n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
# chunk_size <- 1055 # ceiling(n_probes / n_cores) / 2

message("\n=====================================================================")
message("Dataset: ", dataset, 
        "\nArray: 450K & EPICv1",
        "\nTime point: ", timepoint, 
        "\nNormalization: ", normalization)
message("=====================================================================")

if (batch_corrected) {
  input_file_name <-  "metadata"
  plots_color <- 'lightblue'
} else {
  input_file_name <- "metadata_nobatch"
  plots_color <- 'aquamarine4'
}
get_meta_matrix <- function(array, batch_corrected=TRUE) {

  metadt_file <- paste(dataset, array, timepoint, normalization, input_file_name, sep = "_")
  
  desc <- file.path(input_dir, paste0(metadt_file, ".desc"))
  
  message('\nLoading metadata...')
  mat <- as.matrix(bigmemory::attach.big.matrix(desc))
  colnames(mat) <- paste(colnames(mat), array, sep='_')
  
  return(mat)
}

meta_450K <- get_meta_matrix("450K", batch_corrected = batch_corrected)
meta_EPIC <- get_meta_matrix("EPIC", batch_corrected = batch_corrected)

meta <- cbind(meta_450K, meta_EPIC)
rm(meta_450K, meta_EPIC)


# Explore data------------------------------------------------------------------

moment_diff <-function(moment, data = meta, xlim=c(-1,1)) {
  diff <- data[, paste0(moment, '_450K')] - data[, paste0(moment, '_EPIC')]
  print(summary(diff))
  hist(diff, breaks=100, main=moment, xlim=xlim,
       xlab = "Higher in EPIC   <--     Difference     -->   Higher in 450K", 
       col=plots_color, border = plots_color)
  return(diff)
}
plot_moment_diff <- function(){
  par(mfrow = c(2, 2))
  mean_diff <- moment_diff('mean', xlim=c(-1,1))
  range_diff <- moment_diff('range', xlim=c(-1,1))
  skew_diff <- moment_diff('skewness', xlim=c(-40,40))
  kurt_diff <- moment_diff('kurtosis', xlim=c(-1000,1000))
  
}
plot_moment_diff()

# Range against mean value
plot_moments <- function(moment_list = c('mean','range','skewness','kurtosis'),
                         data = meta){
  # Set layout: 2 rows, 4 columns
  par(mfrow = c(2, length(moment_list)))
  
  # First row
  for (array in c('450K', 'EPIC')){
    for (moment in moment_list) {
      stat <- data[,paste(moment, array, sep='_')]
      hist(stat, col = plots_color, border = plots_color, breaks=80, 
           main = paste(array, moment),
           xlab = moment)
      if (moment == 'skewness') {
        abline(v = c(-1, -0.5, 0.5, 1), 
               col = c('pink2','blue','blue','pink2'),
               lwd = 1)
        
        # what percentage is not symmetrical
        percent_bad <- (1 - sum(stat>-1 & stat <1) / length(stat))*100
        
        usr <- par("usr") # this subplot
        x_pos <- usr[2]   # right edge
        y_pos <- usr[4]   # top edge
        
        text(x = x_pos-5, y = y_pos, 
             labels = paste0(round(percent_bad), '%\nasymmetrical'),
             adj = c(0.5, 1), col = "pink4", cex = 1) 
      }
      if (moment == 'kurtosis') {
        abline(v = c(3), col = 'blue', lwd = 1)
        
        # what percentage is leptokurtic
        bad_percent <- sum(stat>3.5) / length(stat) *100
        
        # what percentage is pleiokurtic
        # sum(vkurt<2.5)/length(vkurt)
        
        usr <- par("usr") # this subplot
        x_pos <- usr[2]   # right edge
        y_pos <- usr[4]   # top edge
        
        text(x = x_pos-100, y = y_pos, 
             labels = paste0(round(bad_percent), '%\nleptokurtic'),
             adj = c(0.5, 1), col = "blue4", cex = 1)
        
      }
    }
  }
}
plot_moments()

plot_scatter <- function(x, y, data = meta){
  # Set layout
  par(mfrow = c(1, 2))
  
  for (array in c('450K', 'EPIC')){
    plot(data[,paste(x, array, sep='_')], 
         data[,paste(y, array, sep='_')],
         cex=1, pch=16, col=scales::alpha(plots_color, 0.1),
         xlab = x, ylab =y, main = array)
    
    if (x == 'skewness') {
      abline(v = c(0), col = c('orange'), lwd = 1)
    }
  }
}
# Range against the mean value
plot_scatter('range', 'mean')
# Skewness against mean value
plot_scatter('skewness', 'mean')
plot_scatter('skewness', 'range')

# ==============================================================================

tot_range <- c(meta[,'range_450K'], meta[,'range_EPIC'])

centile_breaks <- quantile(tot_range, probs = seq(0, 1, by = 0.1), na.rm = TRUE)

centile_bins_450K <- findInterval(meta[,'range_450K'], vec = centile_breaks, 
                                  rightmost.closed = TRUE,  all.inside = TRUE)
centile_bins_EPIC <- findInterval(meta[,'range_EPIC'], vec = centile_breaks, 
                                  rightmost.closed = TRUE,  all.inside = TRUE)

ovlp <- table(centile_bins_450K, centile_bins_EPIC)
ovlp

diff <- centile_bins_450K - centile_bins_EPIC
table(diff)

tolerate = 0
message(round(length(diff[ diff %in% c(-tolerate:tolerate)]) /length(diff)*100),'% are ok.')


# calc_range_centiles <- function(vrange, cent=0.01) {
#   
#   centile_breaks <- quantile(vrange, 
#                              probs = seq(0, 1, by = cent), 
#                              na.rm = TRUE)
#   centile_bins <- findInterval(vrange, 
#                                vec = centile_breaks, 
#                                rightmost.closed = TRUE, 
#                                all.inside = TRUE)
#   return(centile_bins)
#   
# }
# 
# range_bins_450K <- calc_range_centiles(K[,'range'], cent=0.05)
# range_bins_EPIC <- calc_range_centiles(E[,'range'], cent=0.05)

ovlp <- table(range_bins_450K, range_bins_EPIC)
ovlp

diff <- range_bins_450K - range_bins_EPIC
table(diff)

message(round(length(diff[ diff %in% c(0)]) /length(diff)*100),'% are ok.')

