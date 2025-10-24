# ==============================================================================
# DNAm matrices tend to be very large and come in various formats, so I will
# first format them and save them as a file-backed big.matrix object
# While at it, I also fetch and format batch (i.e. plate) information

library(dplyr)
library(bigmemory)

outp_path <- '~/MPSR/data/raw'
outp_path_plate_info <- file.path(outp_path,'plate_info')

# dir.create(outp_path, showWarnings = FALSE, recursive = TRUE)
dir.create(outp_path_plate_info, showWarnings = FALSE, recursive = TRUE)

convert_matrix <- function(mat, 
                           dataset = c('GENR', 'ALSP'), 
                           array = c('EPIC','450K'),
                           timepoint = 'birth',
                           normalization = 'funcnorm',
                           output_path = outp_path) {
  
  dataset <- match.arg(dataset)
  array <- match.arg(array)
  
  data_name <- paste(dataset, array, timepoint, normalization, sep = "_")
  
  x <- as.big.matrix(mat, 
                     type = typeof(mat),
                     backingpath = outp_path,
                     backingfile = paste0(data_name, ".bin"),
                     descriptorfile = paste0(data_name, ".desc"))
                     # dimnames = list(rownames(mat), colnames(mat))
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
if (!all(plate_info$Sample_ID == colnames(mat))) stop("Samples need reordering.")

plate_info <- plate_info[match(colnames(mat), plate_info$Sample_ID), ]
if (!all(plate_info$Sample_ID == colnames(mat))) stop("Samples need reordering.")

table(plate_info$Sample_Plate)

saveRDS(plate_info, file = file.path(outp_path_plate_info, 'GENR_EPIC_birth_plates.rds'))

rm(mat)
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
if (!all(plate_info$Sample_ID == colnames(mat))) stop("Samples need reordering.")

table(plate_info$Sample_Plate)

saveRDS(plate_info, file = file.path(outp_path_plate_info, 'GENR_450K_birth_plates.rds'))

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

plate_info_clean <- plate_info %>%
  filter(Sample_Name %in% birth_samples) %>%
  transmute(Sample_ID = as.character(Sample_Name),
            Sample_Plate = as.factor(as.integer(gsub('\\D', '', BCD_plate))))

if (length(unique(plate_info_clean$Sample_ID)) != nrow(plate_info_clean)) stop("Duplicated samples")
if (!all(plate_info_clean$Sample_ID == colnames(mat))) stop("Samples need reordering.")

table(plate_info_clean$Sample_Plate)

saveRDS(plate_info_clean, file = file.path(outp_path_plate_info, 'ALSP_450K_birth_plates.rds'))
