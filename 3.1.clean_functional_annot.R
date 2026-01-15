# ==============================================================================
# Post-processing: clean functional files used for annotation
# ==============================================================================
use_library = '/home/s.defina/R/x86_64-pc-linux-gnu-library/4.4'
.libPaths(use_library)

sumstats_dir <- "~/MPSR/funct_analysis/PrentalRiskFactors_sumstats/"

# What i am working with you guys:
# 
# ├── (E)PACE-HypertensiveDisorders
# │   ├── GestHTEandPE_CpGs.xlsx
# │   ├── hyp_hype201912634d_supp5.xlsx
# │   └── PE_CpGs.xlsx
# ├── PACE-Birthweight
# │   ├── BirthweightEWAS_450kresults_exclCrossReactiveProbes.csv
# │   └── ReadMe.txt
# ├── PACE-C-section
# │   ├── MetaC_Fixed.out
# │   └── MetaC_Fixed_txt.txt
# ├── PACE-Gestationalage
# │   ├── GestationalageEWAS_450kmeta-analysisresult.xlsx
# │   └── ReadMe.txt
# ├── PACE-Gestationaldiabetes
# │   ├── DC190524SupplementaryData2.xlsx
# │   ├── PACE_GDM_Model2_NomSig.csv
# │   └── Thumbs.db
# ├── PACE-Glycemicdysregulation
# │   ├── fastingglucose_EWAS_MRC.txt.gz
# │   ├── fastinginsulin_EWAS_MRC.txt.gz
# │   ├── glucose_EWAS_MRC.txt.gz
# │   ├── insulin_EWAS_MRC.txt.gz
# │   ├── OGTT_EWAS_MRC.txt.gz
# │   └── ReadMe.txt
# ├── PACE-MaternalBMI
# │   ├── BMI_cells_1_random1.Rdata
# │   └── OVERorOBESE_cells_1_all1.Rdata
# ├── PACE-Maternaleducation
# │   ├── mateduewas_summarystatistics.zip
# │   └── ReadMe.txt
# ├── PACE-Maternal&paternalage
# │   ├── acel14194-sup-0005-tabless1-s15_SEE_TABLE_S11.xlsx
# │   ├── DSF12_paternal_newborn_main_EPIC.csv
# │   ├── DSF13_paternal_newborn_main_sexadj.csv
# │   ├── DSF2_maternal_newborn_main.csv
# │   ├── DSF5_maternal_newborn_main_EPIC.csv
# │   ├── DSF6_maternal_newborn_main_sexadj.csv
# │   ├── DSF9_paternal_newborn_main.csv
# │   └── ReadMe.txt
# └── PGS-ADHD,ASDandSCZ(Isabelpaper)
#     ├── adhd_all_elena_2024-07-02.RData
#     ├── adhd_exclGenREpic_elena_2024-07-02.RData
#     ├── asd_all_elena_2024-07-02.RData
#     ├── asd_exclGenREpic_elena_2024-07-02.RData
#     ├── scz_all_elena_2024-07-02.RData
#     └── scz_exclGenREpic_elena_2024-07-02.RData

read_ss_file <- function(file, dir = sumstats_dir) {
  
  if (grepl('.xlsx$', file)) {
    ss_raw <- readxl::read_excel(file.path(dir, file))
  } else if (grepl('.csv$|.txt$|.gz$', file)) {
    ss_raw <- data.table::fread(file.path(dir, file))
  } else if (grepl('.Rdata$|.RData$', file)) {
    tmp_env <- new.env(parent = emptyenv())
    obj_name <- load(file.path(dir, file), envir = tmp_env)
    if (length(obj_name) != 1L) stop("RData file expected to contain one object")
    ss_raw <- tmp_env[[obj_name]]
  }
  
  return(ss_raw)
}

# Hypertensive disorders in pregnancy only have the top 1000-ish CpGs, leave them out
# ss_raw <- read_ss_file('(E)PACE-HypertensiveDisorders/GestHTEandPE_CpGs.xlsx')
# ss_raw <- read_ss_file('(E)PACE-HypertensiveDisorders/PE_CpGs.xlsx')

# ------------------------------------------------------------------------------
ss_birthweight <- read_ss_file(
  'PACE-Birthweight/BirthweightEWAS_450kresults_exclCrossReactiveProbes.csv') |>
  dplyr::select(cpg = MarkerName, pvalue_Birthweight = P.value)

ss_c.section <- read_ss_file('PACE-C-section/MetaC_Fixed_txt.txt') |>
  dplyr::select(cpg = rs_number, `pvalue_C-Section` = `p-value`)

ss_gest.age <- read_ss_file(
  'PACE-Gestationalage/GestationalageEWAS_450kmeta-analysisresult.xlsx') |>
  dplyr::select(cpg = CpGID, pvalue_Gestational.age = PVALUE_FE)

ss_gest.diabetes <- read_ss_file(
  'PACE-Gestationaldiabetes/PACE_GDM_Model2_NomSig.csv') |>
  dplyr::select(cpg = MarkerName, pvalue_Gestational.diabetes = P_value)

# Glycemic dysregulation

ss_glucose.fasting <- read_ss_file(
  'PACE-Glycemicdysregulation/fastingglucose_EWAS_MRC.txt.gz') |>
  dplyr::select(cpg = CpG, pvalue_Glucose.fasting = P_FDR)

ss_insulin.fasting <- read_ss_file(
  'PACE-Glycemicdysregulation/fastinginsulin_EWAS_MRC.txt.gz') |>
  dplyr::select(cpg = CpG, pvalue_Insulin.fasting = P_FDR)

ss_glucose <- read_ss_file(
  'PACE-Glycemicdysregulation/glucose_EWAS_MRC.txt.gz') |>
  dplyr::select(cpg = CpG, pvalue_Glucose = P_FDR)

ss_insulin <- read_ss_file(
  'PACE-Glycemicdysregulation/insulin_EWAS_MRC.txt.gz') |>
  dplyr::select(cpg = CpG, pvalue_Insulin = P_FDR)

ss_OGTT <- read_ss_file(
  'PACE-Glycemicdysregulation/OGTT_EWAS_MRC.txt.gz') |>
  dplyr::select(cpg = CpG, pvalue_OGTT = P_FDR)

# Maternal BMI (... what about the cohort specific effects?)

ss_maternalBMI <- read_ss_file('PACE-MaternalBMI/BMI_cells_1_random1.Rdata') |>
  dplyr::select(cpg = MarkerName, pvalue_Maternal.BMI = `P-value`)
 
ss_maternalBMI2 <- read_ss_file('PACE-MaternalBMI/OVERorOBESE_cells_1_all1.Rdata') |>
  dplyr::select(cpg = MarkerName, pvalue_Maternal.Overweight = `P-value`)

# Maternal education

# There are too many models here... picking one? Skip for now  
# f <- 'PACE-Maternaleducation/mateduewas_summarystatistics.zip'
# files <- unzip((file.path(sumstats_dir, f)), list = TRUE)
# for file in files$Name read_ss_file()?

# Parental age 

# Only has a subset?? EPIC and sexadj models? ...leave them for now
# ss_maternalage <- read_ss_file(
#   'PACE-Maternal&paternalage/DSF2_maternal_newborn_main.csv') |>
#   dplyr::select(cpg = CpG, pvalue_Maternal.age = pval_FDR)

# ss_paternalage <- read_ss_file(
#   'PACE-Maternal&paternalage/DSF9_paternal_newborn_main.csv') |>
#   dplyr::select(cpg = CpG, pvalue_Paternal.age = pval_FDR)

# Mental heath

# exclGenREpic ??? This is also only top hits??

ss_ADHD <- read_ss_file(
  'PGS-ADHD,ASDandSCZ(Isabelpaper)/adhd_exclGenREpic_all_elena_2024-07-02.RData') |>
  dplyr::select(cpg = MarkerName, pvalue_ADHD.PGS = P.value)

ss_ASD <- read_ss_file(
  'PGS-ADHD,ASDandSCZ(Isabelpaper)/asd_exclGenREpic_all_elena_2024-07-02.RData') |>
  dplyr::select(cpg = MarkerName, pvalue_ASD.PGS = P.value)

ss_SCZ <- read_ss_file(
  'PGS-ADHD,ASDandSCZ(Isabelpaper)/scz_exclGenREpic_all_elena_2024-07-02.RData') |>
  dplyr::select(cpg = MarkerName, pvalue_SCZ.PGS = P.value)

# ------------------------------------------------------------------------------
merge_and_rm <- function(df_names, ..., env = parent.frame()) {
  # get the data frames 
  dfs <- mget(df_names, envir = env)
  
  # merge them
  merged <- Reduce(function(x, y) merge(x, y, ...), dfs)
  
  # remove originals from caller's environment
  rm(list = df_names, envir = env)
  
  return(merged)
}

df_list <- grep('^ss_', ls(), value = TRUE)

all_ss <- merge_and_rm(df_list, by = 'cpg', all = TRUE)

summary(all_ss)

# Find sustats files with only significant values  
top_hits_only <- names(all_ss)[sapply(all_ss, function(s) {
  is.numeric(s) && max(s, na.rm = TRUE) <= 0.05
  })]

top_hits_only

# Remove them (for now)
all_ss[top_hits_only] <- NULL

# Add annotation 
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

annot <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) |>
  as.data.frame() |>
  dplyr::select(cpg = Name, chr, pos)

# Why are there NAs?? how to handle them??
sumstats <- merge(annot, all_ss, by = 'cpg', all.x = TRUE)

saveRDS(sumstats, 'summstats_prenatalrisk.rds')

