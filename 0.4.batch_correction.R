use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

# ==============================================================================
# SET-UP
# ==============================================================================

options(startup.messages = FALSE) # avoid Loading package messages on each worker

# Install if needed
# install.packages(c("bigmemory", "foreach", "future", "doFuture", "doRNG", "progressr"))
# BiocManager::install("sva")

# Parallel-safe reproducibility
library(doRNG)
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
# Note: I will use '%dorng%' from the 'doRNG' package instead of '%dopar%'. 
# This ensures that proper, parallel-safe random numbers are produced.

input_dir <- "~/MPSR/data/clean"

output_dir <- "~/MPSR/data/batch_corrected"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
chunk_size <- 1055 # ceiling(n_probes / n_cores) / 2

# ==============================================================================
# BATCH INFO
# ==============================================================================

# Load plate info
load(file.path(input_dir, 'mega_info.RData')) # args, plate_info

# batch <- plate_info[, c("Dataset", "Array", "Plate")]
# This version of ComBat only allows one batch variable
# However, plates are unique per cohort / array
batch <- plate_info$Plate

# ==============================================================================
# BUILD OUPUT MATRIX
# ==============================================================================

# Load filtered mega_matrix
x_desc <- file.path(input_dir, 'mega_filtered.desc')
x <- bigmemory::attach.big.matrix(x_desc)

# Initialize mega matrix
xbc_desc <- file.path(output_dir, "mega_ComBat.desc")

if (file.exists(xbc_desc)) {
  warning(' * "', xbc_desc, '" already exists, removing it.')
  file.remove(xbc_desc)
  file.remove(gsub('.desc', '.bin', xbc_desc))
}

message('Generate (mega) output methylation matrix...')
xbc <- bigmemory::filebacked.big.matrix(
  nrow = nrow(x),
  ncol = ncol(x),
  type = "double",
  backingpath = output_dir,
  descriptorfile = basename(xbc_desc),
  backingfile = gsub('.desc', '.bin', basename(xbc_desc)),
  dimnames = dimnames(x)
)

# ==============================================================================
# RUN BATCH CORRECTION
# ==============================================================================
message('Running analysis...')

# Chunks of CpGs to process in parallel
chunk_seq <- split(rownames(x),
                   ceiling(seq_along(rownames(x)) / chunk_size))

# Set up parallel backend
future::plan(future::multisession, workers = n_cores)
doFuture::registerDoFuture()

# Run ComBat in parallel
foreach::foreach(chunk = chunk_seq,
                 .packages = c("bigmemory", "sva"),
                 .export = c("batch", "x_desc", "xbc_desc"),
                 .combine = 'c') %dorng% {
                   
                   # Note: with future framework I need to re-attach matrices
                   # this adds a little overhead but ensures the memory pointers are not
                   # corrupted -> recommended pipeline on SLURM clusters 
                   x <- bigmemory::attach.big.matrix(x_desc, readonly = TRUE, lockfile = TRUE)
                   xbc <- bigmemory::attach.big.matrix(xbc_desc, lockfile = TRUE)
                   
                   adjusted_chunk <- NULL
                   
                   tryCatch({
                     data_chunk <- x[chunk, ]
                     
                     # Run ComBat (non-parametric to avoid normality assumption)
                     adjusted_chunk <- suppressWarnings(suppressMessages(
                       sva::ComBat(dat = data_chunk, batch = batch, par.prior = FALSE)))
                     
                     
                   }, error = function(e) {
                     message("Chunk failed: ", conditionMessage(e))
                   })
                   
                   if (!is.null(adjusted_chunk)) {
                     # Write corrected values
                     xbc[chunk, ] <- adjusted_chunk
                   }
                   
                   # Suppress result collection
                   invisible(NULL)
                 }

future::plan(future::sequential)

warnings()

# checking ---------------------------------------------------------------------

message('Checking batch correction...')


# TODO 20K most variable (variance )

set.seed(123)
random_cpgs <- sample(rownames(xbc), 1000)
# Check these are the same 
print(random_cpgs[1:3])

x_subset <- x[random_cpgs, ]
xbc_subset <- xbc[random_cpgs, ]

# Compute correlations row-wise
cors <- mapply(function(orig, corr) cor(orig, corr, method='spearman'),
               split(x_subset, row(x_subset)),
               split(xbc_subset, row(xbc_subset)))

# Name and print results nicely
# names(cors) <- random_cpgs
# print(cors_450K)
message('Pre-post correction correlations (sample n = 1000):')
print(summary(cors))
