# ================= Functional analysis (2) EWAS catalog =======================

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

min_association_count <- 10 # Do not look at traits that are associated with < than n CpGs

# Read in clusters 
clusters <- readRDS(file.path(input_dir, 'phase1_clusters.rds'))

## EWAS CATALOG ANALYSES -------------------------------------------------------

# TODO: specify the source
# Load the EWAS Catalog data
ewas_studies <- data.table::fread(file.path(output_dir, 'ewascatalog-studies.txt'), 
                                  select = c('StudyID','Trait'))
ewas_results <- data.table::fread(file.path(output_dir, 'ewascatalog-results.txt'), 
                                  select = c('StudyID','CpG'))

ewas_catalog <- merge(ewas_results, ewas_studies, by='StudyID', all.x=TRUE)
rm(ewas_results, ewas_studies)

# Merge
clust_data <- merge(clusters, ewas_catalog, 
                    by.x = "cpg", by.y = "CpG", all.x = TRUE, all.y=FALSE)
rm(clusters, ewas_catalog)

# Clean-up ---------------------------------------------------------------------

traits <- sort(table(clust_data$Trait), decreasing = TRUE)
message(length(traits), ' unique CpG - trait associations.')

# Cleaning helpers
# traits[1:50] # check most popular associations
# traits[50:100]
# traits[100:150]
# 
# names(traits)[grep('smok', names(traits), ignore.case=TRUE)] # find re-wordings

clean_trait <- function(new_name, old_names){
  # Super-assigned: modifies clust_data in parent environment
  clust_data$Trait[clust_data$Trait %in% old_names] <<- new_name
}

clean_trait('Sex', c('sex', 'Infant sex'))
clean_trait('Age', c('age',
                     'Age 4 vs age 0'))
clean_trait('Gestational age', c('gestational age'))
clean_trait('Birthweight', c('birthweight', 
                             'birth weight',
                             'Birth weight'))
clean_trait('Type 2 Diabetes', c('Incident Type 2 Diabetes',
                                 'Prevalent Type 2 Diabetes (Self-report)',
                                 'Type II diabetes',
                                 'Type II diabetes mellitus',
                                 'Type 2 diabetes'))
clean_trait('Gestational Diabetes', c('Gestational diabetes mellitus',
                                      'gestational diabetes milletus',
                                      'Intrauterine exposure to diabetes'))
clean_trait('Maternal BMI', c('Maternal overweight/obesity',
                              'maternal pre-pregnancy body mass index',
                              'Pre-pregnancy body mass index',
                              'Maternal body mass index',
                              'Maternal underweight'
                              ))
clean_trait('Maternal smoking (pregnancy)', c('Maternal smoking in pregnancy',
                                              'maternal smoking', 'Maternal smoking',
                                              'Maternal smoking during pregnancy',
                                              'Smoking during pregnancy',
                                              'Prenatal smoke exposure'))
clean_trait('Schizophrenia', c('schizophrenia'))
clean_trait('Alzheimers disease', c('alzheimers disease braak stage'))
clean_trait('BMI', c('body mass index',
                     'Body mass index',
                     'Body mass index change'))
clean_trait('Smoking', c('Tobacco smoking',
                         'smoking',
                         'Smoking pack-years'))
clean_trait('Preterm birth', c('Early spontaneous preterm birth',
                               'very early preterm birth (<28 weeks)',
                               'early preterm birth (<28 weeks)',
                               'idiopathic preterm birth',
                               'preterm birth'))
clean_trait('Waist circumference', c('waist circumference',
                                     'Waist circumfrence'))
clean_trait('Gene expression', grep('Gene expression', unique(clust_data$Trait), value=TRUE))
clean_trait('ADHD', c('Attention deficit hyperactivity disorder'))


traits <- sort(table(clust_data$Trait), decreasing = TRUE)
message(length(traits), ' unique CpG - trait associations.')

# cat(names(traits), sep = '\n')
# Cleaning helpers
traits[1:50] # check most popular associations
# traits[50:100]
# traits[100:150]
# 
# names(traits)[grep('smok', names(traits), ignore.case=TRUE)] # find re-wordings

# Duplicated CpGs - trait associations -----------------------------------------
cpg_trait <- paste0(clust_data$cpg, clust_data$Trait)
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

# Plotting 
trait = 'Birthweight' # 'Gestational age'

data_desc <- '~/MPSR/data/batch_corrected/mega_ComBat.desc'
data <- bigmemory::attach.big.matrix(data_desc, lockfile = TRUE)

cpgs_associated <- clust_data$cpg[clust_data$Trait == trait] |>
  unique() |>           # remove duplicates
  (\(x) x[!is.na(x)])() # Ensure no NA

# subset rows of big.matrix by rownames
cpg_idx <- which(rownames(data) %in% cpgs_associated)
data_associated  <- data[cpg_idx, , drop = FALSE] 

cent_associated <- clust_data$cluster

ggplot(
  df_long,
  aes(x = value,
      group = row,          # one density per row
      color = group_var)    # color by other variable
) +
  geom_density() +
  labs(color = "Group") +
  theme_minimal() +
  theme(legend.position = "right")


df_long <- as.data.frame(data_associated) |>
  dplyr::mutate(row = cpgs_associated) |>
  tidyr::pivot_longer(
    cols      = -row,
    names_to  = "id",
    values_to = "value"
  )

ggplot(df_long, aes(x = value, color = row)) +
  geom_density() +
  labs(color = "Row") +
  theme_minimal() + 
  theme(legend.position = "none")


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

