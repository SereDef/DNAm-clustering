
clust_results_path <- './DNAm-clustering/birth_450K_5cent_15maxk'

clust_results_files <- list.files(path = clust_results_path, pattern = '[0-9].csv')


clust_reslts <- data.frame()

for (f in seq_along(clust_results_files)) {
  
  cf = read.csv(file.path(clust_results_path, clust_results_files[f]))
  
  cent <- as.integer(gsub('range-cent', '', 
               stringr::str_split_i(clust_results_files[f], '_', 1)))
  cf$X <- cent
  
  names(cf) <- c('centile','CpG','cluster')
  
  cf$cent_clust <- paste0(sprintf("%02d", cf$centile), '_', 
                          sprintf("%02d", cf$cluster))
  
  cf <- cf[,c('CpG','centile','cluster','cent_clust')]
  
  clust_reslts <- rbind(clust_reslts, cf)
}

clust_reslts <- clust_reslts[order(clust_reslts$cent_clust),]

# table(clust_reslts$cent_clust)

write.csv(clust_reslts, file.path(clust_results_path, 'all_cluster_results.csv'), 
          row.names=FALSE)



