################################################################################
#
# Starting from a file with mutation frequency for every sample, this program 
# does two things:
#
# 1) compute the combined values per series/variant (ALPHA_A1 + ALPHA_A_2 + ...)
# 2) Add the lethal/non-lethal status information

BD_DATA <- ""
DB_DIR <- ""

if( .Platform$OS.type == "windows" ){
  setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
  source("C:/lab/TR_ERRORS/R_files/seqan-test/seqan-functions-v3.R")
  source("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/R_files/getLethals.R")
  source("R_files/covid-functions-v2.R")
  
  BD_DATA = "C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/CAND_OBS_CALLS_BY_VARIANT"
  DB_DIR = "C:/lab/TR_ERRORS/DATA/Coronavirus/DB/ensembl_ASM985889v3/"
  
  fname_lethal_stops <- "C:/lab/TR_ERRORS/DATA/Coronavirus/DB/ensembl_ASM985889v3/res_lethal_stop.bin"
  fname_vcf <- "C:/lab/TR_ERRORS/DATA/Coronavirus/DB/ensembl_ASM985889v3/variants/vcf_usher_ensembl.tab"
  
} else {
  options(width=Sys.getenv("COLUMNS"))
  setwd("/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2")
  source("/mnt/scratch/goutlab/RollingCircle/R_files/seqan-functions-v3.R")
  source("R_files/covid-functions-v2.R")
  
  BD_DATA <- "/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/CAND_OBS_CALLS_BY_VARIANT"
  DB_DIR <- "/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/DB/ensembl_ASM985889v3"
  
  fname_vcf <- "/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/DB/ensembl_ASM985889v3/variants/vcf_usher_ensembl.tab"
}


vBases = getBases()
vMuts = getMutsFromBases(vBases)

# This file contains the unfiltered mutation frequencies for every possible mutation in every sample
fname <- paste(BD_DATA, "all-frequencies.tab.gz", sep="/")
tmuts <- read.table(fname, h=T)

# Extracting the list of sample names from the columns names
vcn <- colnames(tmuts)
wmut <- grepl("_muts", vcn, fixed = TRUE)
wcov <- grepl("_cov", vcn, fixed = TRUE)
vcr <- vcn[which(wmut==T)]
sampleNames <- substr(vcr, 1, nchar(vcr)-5)

# Computing the total across all samples:
tmuts$total_muts <- apply(tmuts[ , wmut], 1, sum)
tmuts$total_cov <- apply(tmuts[ , wcov], 1, sum)

# Computing the totals per variant:
v_variants <- c("USA", "ALPHA", "DELTA")
for(variant in v_variants) {
  pattern <- paste(variant, "_*_muts", sep = "")
  wmut <- grepl(glob2rx(pattern), vcn, fixed = F)
  pattern <- paste(variant, "_*_cov", sep = "")
  wcov <- grepl(glob2rx(pattern), vcn, fixed = F)
  cn <- paste(variant, "muts", sep="_")
  tmuts[ , cn] <- apply(tmuts[ , wmut], 1, sum)
  cn <- paste(variant, "cov", sep="_")
  tmuts[ , cn] <- apply(tmuts[ , wcov], 1, sum)
}


# Computing the totals per series
v_series <- c("USA_A", "USA_B", "ALPHA_A", "ALPHA_B", "DELTA_A", "DELTA_B")
for(series in v_series) {
  pattern <- paste(series, "_*_muts", sep = "")
  wmut <- grepl(glob2rx(pattern), vcn, fixed = F)
  pattern <- paste(series, "_*_cov", sep = "")
  wcov <- grepl(glob2rx(pattern), vcn, fixed = F)
  cn <- paste(series, "muts", sep="_")
  tmuts[ , cn] <- apply(tmuts[ , wmut], 1, sum)
  cn <- paste(series, "cov", sep="_")
  tmuts[ , cn] <- apply(tmuts[ , wcov], 1, sum)
}


# Get a list of lethal mutations
lethals <- sars_get_lethal(
  fname_vcf = fname_vcf,
  fname_stop_pos = fname_stop_pos
)


# Add the lethal/non-lethal status for each mutation
tmuts$lethal = F
tmuts$lethal[which( is.element(tmuts$uid, lethals) == T )] = T

fname <- paste(BD_DATA, "all-frequencies-combined-and-lethal.tab.gz", sep="/")
gzf <- gzfile(fname, "w")
write.table(tmuts, gzf, col.names=T, row.names=F, quote=F, sep="\t")
close(gzf)
