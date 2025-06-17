# ================= Functional analysis (2) EWAS catalog =======================

args <- commandArgs(trailingOnly = TRUE)
comparison <- args[1]

comparison <- match.arg(comparison, c('within_centile', 'across_centiles'))

min_association_count <- 50 # Do not look at traits that are associated with < than n CpGs

clust_results_path <- './birth_450K_5cent_15maxk'
functional_utils_path <- './../funct_analysis'

# Load clusters 
clust_results <- read.csv(file.path(clust_results_path, 'all_cluster_results.csv'))

## Keep only non-cross-reactive probes #TODO: remove them before running clusters?
crossreactive_probes <- unique(maxprobes::xreactive_probes(array_type = "450K"))

clust_results_filtered <- clust_results[!clust_results$CpG %in% crossreactive_probes, ]

rm(clust_results, crossreactive_probes)

## EWAS CATALOG ANALYSES -------------------------------------------------------

# TODO: specify the source
# Load the EWAS Catalog data
ewas_studies <- data.table::fread(file.path(functional_utils_path, 'ewascatalog-studies.txt'), 
                                  select = c('StudyID','Trait'))
ewas_results <- data.table::fread(file.path(functional_utils_path, 'ewascatalog-results.txt'), 
                                  select = c('StudyID','CpG'))

ewas_catalog <- merge(ewas_results, ewas_studies, by='StudyID', all.x=TRUE)
rm(ewas_results, ewas_studies)

# Merge
clust_data <- merge(clust_results_filtered, ewas_catalog, 
                    by = "CpG", all.x = TRUE, all.y=FALSE)
rm(clust_results_filtered, ewas_catalog)

# Clean-up ---------------------------------------------------------------------

# NOTE: this is just the most popular traits MORE CLEANING NEEDED
clust_data$Trait[clust_data$Trait %in% c('sex',
                                         'Infant sex')] <- 'Sex'
clust_data$Trait[clust_data$Trait == 'age'] <- 'Age'
clust_data$Trait[clust_data$Trait == 'gestational age'] <- 'Gestational age'
clust_data$Trait[clust_data$Trait %in% c('birthweight', 
                                         'birth weight',
                                         'Birth weight')] <- 'Birthweight'
clust_data$Trait[clust_data$Trait %in% c('Incident Type 2 Diabetes',
                                         'Prevalent Type 2 Diabetes (Self-report)',
                                         'Type II diabetes',
                                         'Type II diabetes mellitus',
                                         'Type 2 diabetes')] <- 'Type 2 Diabetes'
clust_data$Trait[clust_data$Trait %in% c('Gestational diabetes mellitus',
                                         'gestational diabetes milletus',
                                         'Intrauterine exposure to diabetes')] <- 'Gestational Diabetes Mellitus'
clust_data$Trait[clust_data$Trait %in% c('Maternal overweight/obesity',
                                         'maternal pre-pregnancy body mass index',
                                         'Pre-pregnancy body mass index',
                                         'Maternal body mass index')] <- 'Maternal BMI'
clust_data$Trait[clust_data$Trait %in% c('Maternal smoking in pregnancy',
                                         'maternal smoking', 'Maternal smoking',
                                         'Maternal smoking during pregnancy',
                                         'Smoking during pregnancy',
                                         'Prenatal smoke exposure')] <- 'Maternal smoking (pregnancy)'
clust_data$Trait[clust_data$Trait == 'schizophrenia'] <- 'Schizophrenia' 
clust_data$Trait[clust_data$Trait == 'alzheimers disease braak stage'] <- "Alzheimer's disease"
clust_data$Trait[clust_data$Trait == 'body mass index'] <- 'BMI' 
clust_data$Trait[clust_data$Trait %in% c('Tobacco smoking',
                                         'smoking',
                                         'Smoking pack-years')] <- 'Smoking' 
clust_data$Trait[clust_data$Trait %in% c('Early spontaneous preterm birth',
                                         'very early preterm birth (<28 weeks)',
                                         'early preterm birth (<28 weeks)',
                                         'idiopathic preterm birth',
                                         'preterm birth')] <- 'Preterm birth' 


traits <- sort(table(clust_data$Trait), decreasing = TRUE)
message(length(traits), 'unique CpG - trait associations.')

# Cleaning helpers
# traits[1:50] # check most popular associations
# traits[50:100]
# traits[100:150]
# 
# names(traits)[grep('smok', names(traits), ignore.case=TRUE)] # find re-wordings

# Duplicated CpGs - trait associations -----------------------------------------
cpg_trait <- paste0(clust_data$CpG, clust_data$Trait)
message('Removing ', sum(duplicated(cpg_trait)), ' duplicate associations')

clust_data <- clust_data[!duplicated(cpg_trait),]

# Clean out traits associates with < 10 CpGs -----------------------------------
# This mostly for dimensionality reduction
traits <- sort(table(clust_data$Trait), decreasing = TRUE)

rare_associations <- traits[traits < min_association_count]

message('Removing ', length(rare_associations), ' single CpG - trait associations.')

clust_data <- clust_data[!(clust_data$Trait %in% names(rare_associations)),]

traits <- sort(table(clust_data$Trait), decreasing = TRUE)
message(length(traits), ' traits to examine.')

# Analysis =====================================================================

centiles <- sort(unique(clust_data$centile))

# Set-up parallel processing
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

result_list <- parallel::mclapply(
  names(traits),
  function(trait) {
    
    message(trait)
    
    trait_result_list <- lapply(centiles, function(cent_i) {
      
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
      cent_result_list <- lapply(sort(groups), function(group_k) {
        
        # Confusion matrix -- making sure reference is constant & matrix structure is 2x2
        ct <- table(factor(ifelse(cent_i_df$cent_clust == group_k, group_k, reference),
                           levels = c(reference, group_k)), 
                    factor(ifelse(cent_i_df$Trait == trait,'Yes', 'No'),
                           levels = c('No','Yes')))
        
        tryCatch({
          fisher_test <- stats::fisher.test(ct)
          
          group_result <- data.frame(group = group_k,
                                   comparison = reference,
                                   trait = trait,
                                   estimate = fisher_test$estimate,
                                   CI_lower = fisher_test$conf.int[1],
                                   CI_upper = fisher_test$conf.int[2],
                                   p_value = fisher_test$p.value,
                                   n_yes_group = ct[group_k, 'Yes'],
                                   n_no_group = ct[group_k, 'No'],
                                   n_yes_other = ct[reference, 'Yes'],
                                   n_no_other = ct[reference, 'No'])
        }, error = function(e) {
          message(e$message)
          group_result <- data.frame(group = group_k,
                                     comparison = reference,
                                     trait = trait,
                                     estimate = NA,
                                     CI_lower = NA,
                                     CI_upper = NA,
                                     p_value = NA,
                                     n_yes_group = ct[group_k, 'Yes'],
                                     n_no_group = ct[group_k, 'No'],
                                     n_yes_other = ct[reference, 'Yes'],
                                     n_no_other = ct[reference, 'No'])
        })
        # cat('.')
        return(group_result)
      })
  
    # Combine into a data frame
    cent_i_results <- do.call(rbind, cent_result_list)
    
    return(cent_i_results)
    })
    
    # Combine into a data frame
    trait_t_results <- do.call(rbind, trait_result_list)
    
    return(trait_t_results)
    
}, mc.cores = n_cores)

# Save results 
results <- do.call(rbind, result_list)

# Add FDR correction
results$p_fdr <- p.adjust(p = results$p_value, method = "fdr")


write.csv(results, 
          file = file.path(dirname(clust_results_path),
                           paste0('ewatcatalog_analysis_', comparison, '.csv')), 
          row.names=FALSE)

