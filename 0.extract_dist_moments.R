timepoint <- "birth"
array <- "450K"

setwd("~/GENR3/Methylation")
message('Loading methylation data...')
load(paste0("./GENR_",array,"METH_Norm_Release3/GENR_",array,"METH_Release3_Betas_ALL_",timepoint,"_20190813.RData"))

output_folder <- "~/MPSR/"

# ==============================================================================
# Remove X, Y chromosomes and control probes (from annot files )

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # Chromosome annotation

# load annotation data
data(list = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
data(Locations)
data(Other)
annotation <- cbind(as.data.frame(Locations), as.data.frame(Other))

XY_chr_probes = rownames(annotation)[annotation$chr %in% c('chrX', 'chrY')]

other_chr_probes = rownames(annotation)[!annotation$chr %in% c('chrX', 'chrY')]
rm(annotation, Locations, Other)

message("Removing X and Y chromosomes...")
dim(x)
x <- x[!rownames(x) %in% XY_chr_probes, ]
dim(x)
x <- x[!rownames(x) %in% other_chr_probes, ]

## load cross-reactive probes from the max probes dataset for both 450k and epic
# library(maxprobes)
# xloci_EPIC <- unique(unlist(maxprobes::xreactive_probes(array_type = "EPIC")))
# xloci_450K <- unique(maxprobes::xreactive_probes(array_type = "450K"))
# xloci <- unique(c(xloci_EPIC, xloci_450K))
# 
# ## remove from keep the cross-reactive probes
# keep <- keep[!(keep %in% xloci)]


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

# ==============================================================================
# Explore data

library(moments)

message('Computing metadata...')
# Range ------------------------------------------------------------------------
vrange <- apply(x, 1, function(v) diff(range(v, na.rm=TRUE)))
summary(vrange)
# hist(vrange, breaks=100, xlim=c(0,1))

vmeans <- apply(x, 1, mean, na.rm=TRUE)
# hist(vmeans, breaks=100, xlim=c(0,1))

# Range against mean value
# plot(vrange, vmeans, cex=1, pch=16, col=scales::alpha("aquamarine4",0.1),
#      xlab = "Variable range", ylab = "Variable mean")

# Skewness (= asymmetry) -------------------------------------------------------
# Negative skew = tail is on the left side (towards negative values).
# Positive skew = tail is on the right side (towards positive values).

# skew <- function(v) {
#   n <- length(v)
#   vmean <- mean(v, na.rm=TRUE)
#   vsd <- sd(v, na.rm=TRUE)
# 
#   sk <- (n * sum((v - vmean)^3)) / ((n - 1) * (n - 2) * vsd^3)
#   return(sk)
# }

# vskew <- apply(sub, 1, skew)
vskew <- apply(x, 1, moments::skewness, na.rm=TRUE)
summary(vskew)

# what percentage is not symmetrical
# 1 - sum(vskew>-1 & vskew <1)/length(vskew)

# hist(vskew, breaks=100, xlim=c(-15,15))
# abline(v = c(-1, -0.5, 0.5, 1), col = c('pink','blue','blue','pink'),
#        lwd = 1.7)

# Skewness against mean value
# plot(vskew, vrange, cex=1, pch=16, col=scales::alpha("aquamarine4",0.1),
#      xlab = "Variable skewness", ylab = "Variable range")

# Kurtosis ( = heavy- or light-tailed relative to a normal distribution) -------
# Relative to a normal distribution (kurtosis = 3)
# < 3 = playkurtic (fewer and less extreme outliers)
# > 3 = leptokurtic (more outliers)

# kurt <- function(v) {
# 
#   vmean <- mean(v, na.rm=TRUE)
#   vsd <- sd(v, na.rm=TRUE)
# 
#   kt <- sum((v - vmean)^4) / (99 * vsd^4) - 3
#   return(kt)
# }

# vkurt <- apply(sub, 1, kurt)
vkurt <- apply(x, 1, moments::kurtosis, na.rm = TRUE)
summary(vkurt)
# hist(vkurt, breaks=100, xlim=c(0,280))
# abline(v = c(3), col = 'blue', lwd = 1.5)

# what percentage is leptokurtic
# sum(vkurt>3.5)/length(vkurt)
# sum(vkurt<2.5)/length(vkurt)

metad <- data.frame(vrange, vskew, vkurt, vmeans)
row.names(metad) <- row.names(x)

# clean
rm(vrange, vskew, vkurt, vmeans)

message("Transforming dataframe...")

cat(dim(x))
# transform data to list of vectors
# ---> does not handle missing data (TMP: removing it)
df_compl <- as.data.frame(t(x[rowSums(is.na(x))==0,]))

cat(dim(df_compl))

rm(x)

save(df_compl, file=file.path(output_folder, paste0("data_filtered_",array,"_",timepoint,".RData")))

# metadata
clean_metad <- metad[names(df_compl), ]
rm(metad)

write.csv(clean_metad, file.path(output_folder, paste0("distrib_metad_",array,"_",timepoint,".csv")))
