use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

# Install if needed
# install.packages(c("bigmemory", "moments","foreach", "future", "doFuture", "progressr"))

# Parallel-safe reproducibility
library(doRNG)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# =================================== CONFIG ===================================

args <- commandArgs(trailingOnly = TRUE)
which_matrix <- args[1] # before or after batch correction

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
chunk_size <- 1055 # ceiling(n_probes / n_cores) / 2

input_dir <- file.path('~/MPSR/data', which_matrix)
output_dir <- "~/MPSR/metadata"

if (which_matrix == "clean") {
  input_file_name <- "mega_filtered"
  output_file_name <- "metadata_mega"
} else {
  input_file_name <- "mega_ComBat"
  output_file_name <- "metadata_mega_ComBat"
}
# methyl_file <- paste(dataset, array, timepoint, normalization, input_file_name, sep = "_")
# metadt_file <- paste(dataset, array, timepoint, normalization, output_file_name, sep = "_")

x_desc <- file.path(input_dir, paste0(input_file_name, ".desc"))
meta_desc <- file.path(output_dir, paste0(output_file_name, ".desc"))

message("\n=====================================================================")
# message("Dataset: ", dataset, 
#         "\nArray: ", array, 
#         "\nTime point: ", timepoint, 
#         "\nNormalization: ", normalization)
message(which_matrix)
message("=====================================================================")

message('\nLoading methylation data...')
x <- bigmemory::attach.big.matrix(x_desc)
message(' * Data dimentions:')
dim(x)

# Load subset info
load('~/MPSR/data/clean/mega_info.RData') # args, plate_info

# ==============================================================================
# Explore data

compute_moments <- function(data) {
  
  # * Mean
  vmeans <- apply(data, 1, mean, na.rm = TRUE)
  # * Median
  vmedians <- apply(data, 1, stats::median, na.rm = TRUE)
  # * Range 
  vrange <- apply(data, 1, function(v) diff(range(v, na.rm = TRUE)))
  # * Variance
  vvariances <- apply(data, 1, stats::var, na.rm = TRUE)
  # * Skewness (= asymmetry)
  #   Negative skew = tail is on the left side (towards negative values)
  #   Positive skew = tail is on the right side (towards positive values)
  vskew <- apply(data, 1, moments::skewness, na.rm = TRUE)
  # * Kurtosis (= heavy- or light-tailed)
  #   Relative to a normal distribution (kurtosis = 3)
  #   < 3 = playkurtic (fewer and less extreme outliers)
  #   > 3 = leptokurtic (more outliers)
  vkurt <- apply(data, 1, moments::kurtosis, na.rm = TRUE)
  
  return(cbind(vmeans, vmedians, vrange, vvariances, vskew, vkurt))
}

moment_list <- c('mean', 'median', 'range', 'variance', 'skewness', 'kurtosis')

subset_list <- c('mega', paste0(args$dataset, args$array))

all_moments <- do.call(paste, c(expand.grid(moment_list, subset_list), sep='_'))

cpg_list <- rownames(x)

# Clean up previous matrices if necessary
if (!file.exists(meta_desc)) {
  # Create output big.matrix
  n_probes <- length(cpg_list)
  n_moments <- length(all_moments)
  
  message('Generate output metadata matrix...')
  corrected_x <- bigmemory::filebacked.big.matrix(
    nrow = n_probes,
    ncol = n_moments,
    type = "double",
    backingpath = output_dir,
    descriptorfile = basename(meta_desc),
    backingfile = gsub('.desc', '.bin', basename(meta_desc)),
    dimnames = list(cpg_list, all_moments)
  )
}

# Chunks of CpGs to process in parallel
chunk_seq <- split(cpg_list,
                   ceiling(seq_along(cpg_list) / chunk_size))

# Parallel-safe reproducibility
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Progress updates do not really work on this SLURM set-up...
message('Computing metadata...')
# Set up parallel backend
future::plan(future::multisession, workers = n_cores)
doFuture::registerDoFuture()

foreach(chunk = chunk_seq,
        .packages = c("bigmemory", "moments"),
        .export = c("x_desc", "meta_desc"),
        .combine = 'c') %dorng% {
          
          # Note: with future framework I need to re-attach matrices
          # this adds a little overhead but ensures the memory pointers are not
          # corrupted -> recommended pipeline on SLURM clusters 
          x <- bigmemory::attach.big.matrix(x_desc, readonly = TRUE, lockfile = TRUE)
          meta_x <- bigmemory::attach.big.matrix(meta_desc, lockfile = TRUE)
        
          data_chunk <- x[chunk, ]
          
          # Extract moments for the full and the subsets 
          # TODO: generalize to any number of subsets
          chunk_moments <- cbind(
            compute_moments(data_chunk), # everything
            compute_moments(data_chunk[, args[1,'start']:args[1,'end']]), #GENR 450K
            compute_moments(data_chunk[, args[2,'start']:args[2,'end']]), #ALSP 450K
            compute_moments(data_chunk[, args[3,'start']:args[3,'end']]) #GENR EPIC
          )
          
          # Write corrected values
          meta_x[chunk, ] <- chunk_moments
          
          # Suppress result collection
          invisible(NULL)
        }

warnings()

# Close parallel process
future::plan(future::sequential)

message("Done!\nSummary:")

meta_x <- bigmemory::attach.big.matrix(meta_desc, lockfile = TRUE)
for (m in all_moments) {
  message(m) 
  print(summary(meta_x[, m]))
  cat('\n')
}

# ==============================================================================

# Simple distribution exploration 
# Print distributions 

# pdf(file.path(output_folder, 'DNAm_distrib_birth.pdf'))
# 
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
# dev.off()


