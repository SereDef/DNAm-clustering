# ==============================================================================
# Find overlapping CpGs between 450k and EPIC

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("minfi", 
                       "IlluminaHumanMethylation450kanno.ilmn12.hg19", 
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Load 450k and EPIC annotation
anno450k <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annoEPIC <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Get CpG IDs
cpg450k <- rownames(anno450k)
cpgEPIC <- rownames(annoEPIC)

# 1. Overlapping CpGs between arrays
overlapping_cpgs <- intersect(cpg450k, cpgEPIC)
length(overlapping_cpgs)  # ~450,000

# 2. Filter out probes on X and Y chromosomes
autosomal_cpgs <- overlapping_cpgs[
  !(anno450k[overlapping_cpgs, "chr"] %in% c("chrX", "chrY"))
]
length(autosomal_cpgs)

# 3. Remove cross-reactive probes ----------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("perishky/maxprobes")

## Keep only non-cross-reactive probes
xreactive_probes_450k <- unlist(unique(maxprobes::xreactive_probes(
  array_type = '450K')))
xreactive_probes_EPIC <- unlist(unique(maxprobes::xreactive_probes(
  array_type = 'EPIC')))

xreactive_probes <- union(xreactive_probes_450k, xreactive_probes_EPIC)

final_cpgs <- setdiff(autosomal_cpgs, xreactive_probes)
length(final_cpgs)

# detection P value..?
# any(grepl('^rs', final_cpgs)) # FALSE
  
writeLines(final_cpgs, "~/MPSR/data/clean_cpg_list.txt")

# NOTE: 425 cpgs in the list were not present in the EPIC file...
# removing them for now, till I figure out the root of the problem 
# x <- attach.big.matrix("GENR_EPIC_birth_funcnorm.desc")
# cpgs_to_keep <- setdiff(cpgs_to_keep, setdiff(cpgs_to_keep, rownames(x)))
# writeLines(cpgs_to_keep, "~/MPSR/data/clean_cpg_list.txt")

