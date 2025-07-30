# ==============================================================================
# DNAm matrices tend to be very large and come in various formats, so I will
# first format them and save them as a file-backed big.matrix object
# While at it, I also fetch and format batch (i.e. plate) information

library(dplyr)
library(bigmemory)

convert_matrix <- function(mat, 
                           dataset = c('GENR', 'ALSPAC'), 
                           array = c('EPIC','450K'),
                           timepoint = 'birth',
                           normalization = 'funcnorm',
                           outp_path = '~/MPSR/data') {
  
  dataset <- match.arg(dataset)
  array <- match.arg(array)
  
  data_name <- paste(dataset, array, timepoint, normalization, sep = "_")
  
  x <- as.big.matrix(mat, 
                     type = type(mat),
                     backingpath = outp_path,
                     backingfile = paste0(data_name, ".bin"),
                     descriptorfile = paste0(data_name, ".desc"))
  invisible(NULL)
}

# ==============================================================================
# GENR EPIC @birth ~ functionally normalized together with ALSPAC

og_methyldata_file <- '/home/projects/EPIC_FC/GENR_EPIC/EPIC_release1_birth_funcnorm_randslide_BETA_240124.rds'
mat <- readRDS(og_methyldata_file)

convert_matrix(mat,
               dataset = "GENR",
               array = "EPIC",
               timepoint ="birth",
               normalization ="funcnorm")

# Plate info
selection_file <- '~/selection_files/Selection_GENR_MethylEPIC_release1_birth_20230717.sav'

plate_info <- haven::read_spss(selection_file, col_select = c('SampleID','Sample_Plate')) %>%
  transmute(Sample_ID = as.character(SampleID), 
            Sample_Plate = as.factor(Sample_Plate))

if (length(unique(plate_info$Sample_ID)) != nrow(plate_info)) stop("Duplicated samples")

table(plate_info$Sample_Plate)

saveRDS(plate_info, file = '~/MPSR/data/GENR_EPIC_birth_plates.rds')

# ==============================================================================
# GENR 450K @birth ~ functionally normalized together with ALSPAC

og_methyldata_file <- '/home/projects/EPIC_FC/GENR_450K/ch3_genr_180321.R'
load(og_methyldata_file) # 1 object: ch3_genr

# selection_file <- '~/selection_files/Selection_GENR_450kmeth_release3_birth_20230608.sav'
# plate_info <- haven::read_spss(selection_file, col_select =c('Sample_ID','Sample_Plate'))

selection_file <- '/home/projects/EPIC_FC/GENR_450K/jointQC_samplesfile_180529.RData'
load(selection_file) # 1 object: samples

birth_samples <- samples[samples$Period=="Birth",]
rm(samples)

# Keep only birth samples (and turn to matrix)
mat <- as.matrix(ch3_genr[, as.character(birth_samples$Sample_ID)])
rm(ch3_genr)

convert_matrix(mat = mat,
               dataset = "GENR",
               array = "450K",
               timepoint ="birth",
               normalization ="funcnorm")

# Plate info
plate_info <- birth_samples %>%
  transmute(Sample_ID = as.character(Sample_ID),
            Sample_Plate = as.factor(as.character(Sample_Plate)))

if (length(unique(plate_info$Sample_ID)) != nrow(plate_info)) stop("Duplicated samples")

table(plate_info$Sample_Plate)

saveRDS(plate_info, file = '~/MPSR/data/GENR_450K_birth_plates.rds')

# ==============================================================================
# ALSPAC 450K @birth ~ functionally normalized together with GENR

og_methyldata_file <- '/home/projects/EPIC_FC/ALS_450K/meth.rds'
mat <- readRDS(og_methyldata_file)

selection_file <- '/home/projects/EPIC_FC/ALS_450K/samplesheet.csv'
plate_info <- read.csv(selection_file)

birth_samples <- plate_info[plate_info$time_point == 'cord', 'Sample_Name']

mat <- mat[, birth_samples]

convert_matrix(mat = mat,
               dataset = "ALSP",
               array = "450K",
               timepoint ="birth",
               normalization ="funcnorm")

# How to use plate info??
