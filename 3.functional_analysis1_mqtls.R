# ===================== Functional analysis (1) mQTLs ==========================

args <- commandArgs(trailingOnly = TRUE)
comparison <- args[1]

comparison <- match.arg(comparison, c('within_centile', 'across_centiles'))

clust_results_path <- './birth_450K_5cent_15maxk'
functional_utils_path <- './../funct_analysis'

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Load annotation data 
# data(list = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(Locations)
data(Other)
annotation <- cbind(as.data.frame(Locations), as.data.frame(Other))
annotation$CpG <- rownames(annotation)
rm(Locations, Other)

# Load clusters 
clust_results <- read.csv(file.path(clust_results_path, 'all_cluster_results.csv'))

# Merge 
clust_annotated <- merge(clust_results, annotation, by = 'CpG', all.x = TRUE)

rm(clust_results, annotation)

## mQTL ANALYSES ---------------------------------------------------------------

# TODO: specify the source
mqtls <- data.table::fread(file.path(functional_utils_path, 'mqtls.cvs.gz'))

mqtls <- mqtls[!duplicated(mqtls$cpg),] # clean-up duplicates

mqtls$mQTL <- as.factor(ifelse((mqtls$clumped == TRUE) & (mqtls$pval < 9*10^-8), 
                               'Yes', 'No'))

cat(as.integer(table(mqtls$mQTL)['Yes']), 'mQTLs identified.')

# Make sure 0 = No
mqtls$mQTL <- relevel(mqtls$mQTL, ref = 'No')

# Merge
clust_data <- merge(clust_annotated, mqtls, 
                    by.x = "CpG", by.y = 'cpg', all.x = TRUE)
rm(clust_annotated, mqtls)

table(clust_data$mQTL, useNA = "ifany") # TODO: Leave NA in there or set at No? 


# Analysis =====================================================================

# Always chunk by centile 
centiles <- sort(unique(clust_data$centile))

# Set-up parallel processing
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

result_list <- parallel::mclapply(
  centiles, 
  function(cent_i) {
    
    # Define type of comparison (i.e. within or across all centiles)
    if (comparison == 'within_centile') {
  
      # Subset centile data
      cent_i_df <- clust_data[clust_data$centile == cent_i, ]
      
      groups <- unique(cent_i_df$cent_clust)
      reference <- paste0(sprintf("%02d", cent_i), '_other')
      
    } else {
      
      cent_i_df <- clust_data
      
      groups <- unique(clust_data[clust_data$centile == cent_i, ]$cent_clust)
      reference <- 'all_other'
    }
    
    # Run Fisher test for each k group in the centile
    cent_result_list <- lapply(groups, function(group_k) {
      
      # Confusion matrix
      ct <- table(relevel(as.factor(ifelse(cent_i_df$cent_clust == group_k, 
                                           group_k, reference)), 
                          ref = reference), # making sure reference is constant
                  cent_i_df$mQTL)
      
      fisher_test <- stats::fisher.test(ct)
      
      group_result <- data.frame(group = group_k,
                                 comparison = reference,
                                 estimate = fisher_test$estimate,
                                 CI_lower = fisher_test$conf.int[1],
                                 CI_upper = fisher_test$conf.int[2],
                                 p_value = fisher_test$p.value,
                                 n_yes_group = ct[group_k, 'Yes'],
                                 n_no_group = ct[group_k, 'No'],
                                 n_yes_other = ct[reference, 'Yes'],
                                 n_no_other = ct[reference, 'No'])
      
      return(group_result)
    })

    # Combine into a data frame
    cent_i_results <- do.call(rbind, cent_result_list)
    
}, mc.cores = n_cores)

# Save results 
results <- do.call(rbind, result_list)

# Add FDR correction
results$p_fdr <- p.adjust(p = results$p_value, method = "fdr")


write.csv(results, file.path(dirname(clust_results_path), 
                             paste0('mqtl_analysis_', comparison, '.csv')), 
          row.names=FALSE) 


# Summarize results
# df_intracentile_mqtl_summary <- df_intracentile_mqtl_results %>%
#   summarize(sign_tot = sum(p_value < 0.05),
#             sign_pos = sum(p_value < 0.05 & estimate_odds_ratio > 1),
#             sign_neg = sum(p_value < 0.05 & estimate_odds_ratio < 1),
#             insign_tot = sum(p_value >= 0.05),
#             insign_pos = sum(p_value >= 0.05 & estimate_odds_ratio > 1),
#             insign_neg = sum(p_value >= 0.05 & estimate_odds_ratio < 1),
#             fdr_sign_tot = sum(pv_fdr < 0.05),
#             fdr_sign_pos = sum(pv_fdr < 0.05 & estimate_odds_ratio > 1),
#             fdr_sign_neg = sum(pv_fdr < 0.05 & estimate_odds_ratio < 1),
#             fdr_insign_tot = sum(pv_fdr >= 0.05),
#             fdr_insign_pos = sum(pv_fdr >= 0.05 & estimate_odds_ratio > 1),
#             fdr_insign_neg = sum(pv_fdr >= 0.05 & estimate_odds_ratio < 1)
#   )

