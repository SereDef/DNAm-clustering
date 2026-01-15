# ===================== Functional analysis (1) mQTLs ==========================

use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

# I process each centile as a task
array_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID', unset = 1))
# I use CPUs to run iterations in parallel
n_cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', unset = 1))

# Set-up
input_dir <- "~/MPSR/metadata/clusters"
output_dir <- "~/MPSR/funct_analysis"

comparison_types <- c('within_range_centile', 'across_centiles')

comparison <- comparison_types[array_id]

# Read in clusters 
clusters <- readRDS(file.path(input_dir, 'phase1_clusters.rds'))

## mQTL ANALYSES ---------------------------------------------------------------

# TODO: specify the source
mqtls <- data.table::fread(file.path(output_dir, 'mqtls.cvs.gz'))

mqtls <- mqtls[!duplicated(mqtls$cpg),] # clean-up duplicates

mqtls$mQTL <- as.factor(ifelse((mqtls$clumped == TRUE) & (mqtls$pval < 9*10^-8), 
                               'Yes', 'No'))

cat(as.integer(table(mqtls$mQTL)['Yes']), 'mQTLs identified.')

# Make sure 0 = No
mqtls$mQTL <- relevel(mqtls$mQTL, ref = 'No')

table(mqtls$mQTL, useNA = "ifany") # TODO: Leave NA in there or set at No? 

# waldo_overlap <- function(x, y, max_diffs = 5) {
#   cat("Overlap summary\n")
#   cat("---------------\n")
#   cat("length(x):", length(x), "\n")
#   cat("length(y):", length(y), "\n")
#   cat("intersect:", length(intersect(x, y)), "\n")
#   cat("only in x:", length(setdiff(x, y)), "\n")
#   cat("only in y:", length(setdiff(y, x)), "\n\n")
#   
#   print(waldo::compare(x, y, max_diffs = max_diffs))
# }

# example
# waldo_overlap(clusters_annotated$cpg, mqtls$cpg)

# Merge
clust_data <- merge(clusters, mqtls, by = 'cpg', all.x = TRUE)
rm(clusters, mqtls)

table(clust_data$mQTL, useNA = "ifany") # TODO: Leave NA in there or set at No? 

# Analysis =====================================================================

clusters <- sort(unique(clust_data$cluster))

# Set-up parallel processing
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

result_list <- parallel::mclapply(
  clusters, 
  function(clust_i) {
    
    # Define type of comparison (i.e. within or across all centiles)
    if (comparison == 'within_range_centile') {
  
      # Subset centile data
      range_centile <- stringr::str_extract(clust_i, "^[0-9]+")
      
      cent_i_df <- clust_data[stringr::str_starts(clust_data$cluster, 
                                                  pattern = range_centile), ]

      reference <- paste0(range_centile, '_other')
      
    } else {
      cent_i_df <- clust_data
      reference <- 'all_other'
    }
    
    # Run Fisher test for each cluster
    ct <- table(relevel(as.factor(ifelse(cent_i_df$cluster == clust_i, clust_i, reference)),
                        ref = reference), # making sure reference is constant
                  cent_i_df$mQTL)
      
    fisher_test <- stats::fisher.test(ct)
      
    group_result <- data.frame(cluster = clust_i,
                               comparison = reference,
                               estimate = fisher_test$estimate,
                               CI_lower = fisher_test$conf.int[1],
                               CI_upper = fisher_test$conf.int[2],
                               p_value = fisher_test$p.value,
                               n_yes_group = ct[clust_i, 'Yes'],
                               n_no_group = ct[clust_i, 'No'],
                               n_yes_other = ct[reference, 'Yes'],
                               n_no_other = ct[reference, 'No'])
      
      return(group_result)

}, mc.cores = n_cores)

# Save results 
results <- do.call(rbind, result_list)
rownames(results) <- NULL

# Add FDR correction
results$p_fdr <- p.adjust(p = results$p_value, method = "fdr")

write.csv(results, file.path(output_dir, 
                             paste0('mqtl_analysis_', comparison, '.csv')), 
          row.names=FALSE) 
