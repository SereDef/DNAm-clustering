use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

output_dir <- "~/MPSR/metadata"

meta_desc <- file.path(output_dir, 'metadata_mega.desc')
meta_desc_ComBat <- file.path(output_dir, 'metadata_mega_ComBat.desc')

meta_filter <- bigmemory::attach.big.matrix(meta_desc, lockfile = TRUE)
meta_ComBat <- bigmemory::attach.big.matrix(meta_desc_ComBat, lockfile = TRUE)

based_on = 'range'

moment_sub <- meta_ComBat[, c(paste0(based_on,'_mega'), 'median_mega')]

message('N features: ', nrow(moment_sub))
print(summary(moment_sub))

# clean up range 
if (based_on == 'range') {
  moment_sub <- moment_sub[moment_sub[,'range_mega'] <= 1, ]
  moment_sub <- moment_sub[is.finite(moment_sub[,'range_mega']), ]
}

message('N features: ', nrow(moment_sub))
print(summary(moment_sub))

# Range centiles 
centile_breaks1 <- quantile(moment_sub[,'range_mega'], 
                           probs = seq(0, 1, by = 0.05), 
                           na.rm = TRUE)

centiles_bins1 <- findInterval(moment_sub[,'range_mega'], vec = centile_breaks1,
                             rightmost.closed = TRUE,  all.inside = TRUE)

# moment_sub <- cbind(moment_sub, centiles_bins1)

centile_bins <- vector(mode='character', nrow(moment_sub))

for (c1 in 1:max(centiles_bins1)){

  c1_idx = which(centiles_bins1 == c1)
  
  # Median centiles 
  c1_median <- moment_sub[c1_idx, 'median_mega']
  
  centile_breaks2 <- quantile(c1_median, 
                              probs = seq(0, 1, by = 0.05), 
                              na.rm = TRUE)
  centile_bins2 <- findInterval(c1_median, vec = centile_breaks2,
                                rightmost.closed = TRUE,  all.inside = TRUE)
  
  centile_bin_names <- paste(sprintf("%02d", c1), 
                             sprintf("%02d", centile_bins2), sep='.')
  
  centile_bins[c1_idx] <- centile_bin_names
  
}

table(centile_bins)

moment_sub <- cbind(moment_sub, centile_bins)

saveRDS(moment_sub, file=file.path(output_dir, 'metadata_mega_ComBat_binned.rds'))

# Median centiles 
# centile_breaks2 <- quantile(moment_sub[,'median_mega'], 
#                             probs = seq(0, 1, by = 0.05), 
#                             na.rm = TRUE)
# centile_bins2 <- findInterval(moment_sub[,'median_mega'], vec = centile_breaks2,
#                              rightmost.closed = TRUE,  all.inside = TRUE)

