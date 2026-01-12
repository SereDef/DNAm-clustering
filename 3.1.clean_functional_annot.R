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

# I start this this because it has position (I think...?)
f <- 'PACE-Gestationalage/GestationalageEWAS_450kmeta-analysisresult.xlsx'
ss_raw <- readxl::read_excel(file.path(sumstats_dir, f))

ss <- ss_raw |> 
  dplyr::select(cpg = CpGID, chr = CHR, pos = MAPINFO_Hg19, 
                pvalue_GestationalAge = PVALUE_FE)

# ------------------------------------------------------------------------------
f <- 'PACE-Birthweight/BirthweightEWAS_450kresults_exclCrossReactiveProbes.csv'
ss_raw <- data.table::fread(file.path(sumstats_dir, f))

ss <- ss_raw |>
  dplyr::select(cpg = MarkerName, pvalue_Birthweight = P.value) |>
  dplyr::full_join(ss, by='cpg')

# ------------------------------------------------------------------------------
f <- 'PACE-Maternaleducation/mateduewas_summarystatistics.zip'
ss_raw <- unzip((file.path(sumstats_dir, f)), list = TRUE)

# ff1 <- "PACE-Birthweight/BirthweightEWAS_450kresults_exclCrossReactiveProbes.csv"

# ss_raw <- read.csv(file.path(sumstats_dir, ff1))

