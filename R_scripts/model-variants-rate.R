################################################################################
#
# This script attempts to model how variation in mutation rate can be explained 
# by variation in the frequency of mutations inside the ORF1ab.
# The goal is to find variants of ORF1ab that are likely to cause increased or 
# decreased mutation rate.

if( .Platform$OS.type == "windows" ){
  setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
  source("C:/lab/TR_ERRORS/R_files/seqan-test/seqan-functions-v3.R")
  source("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/R_files/getLethals.R")
  source("R_files/covid-functions-v2.R")
  
  BD_DATA = "C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/"
  DB_DIR = "C:/lab/TR_ERRORS/DATA/Coronavirus/DB/ensembl_ASM985889v3/"

} else {
  options(width=Sys.getenv("COLUMNS"))
  setwd("/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2")
  source("/mnt/scratch/goutlab/RollingCircle/R_files/seqan-functions-v3.R")
  source("R_files/covid-functions-v2.R")
  
  BD_DATA <- "/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/"
  DB_DIR <- "/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/DB/ensembl_ASM985889v3"
}


library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(genomation)
library(Biostrings)
library(tidyr)
library(readxl)

vBases = getBases()
vMuts = getMutsFromBases(vBases)

################################################################################
# This file contains the frequencies of every mutation in every sample
fAllFrequencies <- paste(BD_DATA, "all-frequencies.csv", sep = "/")
tAll <- as_tibble( read.csv(fAllFrequencies) )
cNames <- colnames(tAll)
cnCov <- cNames[grep("_cov$", colnames(tAll) )]
sampleNames <- substr(cnCov, 1, nchar(cnCov)-4)
for(sampleName in sampleNames){
  cn <- paste(sampleName, "_freq", sep = "")
  cnCov <- paste(sampleName, "cov", sep = "_")
  cnMut <- paste(sampleName, "muts", sep = "_")
  tAll[ , cn] <- tAll[ , cnMut] / tAll[ , cnCov]
}

# Keeping only variants inside ORF1ab
tAll <- tAll[which(tAll$position < 21560) , ]
# Keeping only nonsynonymous (and nonsense) variants:
load("C:/lab/TR_ERRORS/DATA/Coronavirus/DB/protein_predict.bin")
res$uid <- paste( start(res), as.character(res$REF), as.character(res$varAllele), sep = "_" )
uid_ns <- res$uid[which(res$CONSEQUENCE!="synonymous")]
tAll <- tAll[which( is.element(tAll$uid, uid_ns)==T) , ]

################################################################################
# Reading the mutation spectrum for all the samples
ts <- read.csv("spectrum-lethal_COV-1000_RATE-0.001_HS-10.csv")

Variants <- c("ALPHA", "USA", "DELTA")
Replicates <- c("A", "B")


################################################################################
# Keeping only positions that are covered by at least 10 reads, have a relative 
# change in frequency of at least 20% and reach an absolute frequency of at 
# least 5% in one of the passages.
MIN_COV <- 10
minDelta <- 0.2
minAbsolute <- 0.05

tres <- data.frame(
  position = numeric(),
  uid = character(),
  cor = numeric(),
  pval = numeric(),
  best = logical(),
  seriesName = character(),
  varMut = numeric()
)

# Now, computing the correlations
for(Variant in Variants){
  for(Replicate in Replicates){
    seriesName <- paste(Variant, Replicate, sep = "_")

    tpos <- getChangingAlleleFrequencies(tAll = tAll, seriesName = seriesName, 
                                         MIN_COV = MIN_COV, minDelta = minDelta, minAbsolute = minAbsolute)
    
    
    sPattern <- paste(seriesName, ".*_freq", sep = "")
    tpr <- tpos[ , c(1:4, grep(sPattern, colnames(tpos)))]
    
    tsr <- ts[which(ts$variant==Variant & ts$replicate==Replicate & is.na(ts$passage)==F & ts$mutation=="N->N") , c("passage", "rate")]
    
    
    cnFreq <- grep(sPattern, colnames(tpos))
    tpr <- tpos[ , c(1, 2, cnFreq)]
    #colnames(tpr)[2:ncol(tpr)] <- seq(from=1, to=(ncol(tpr)-1), by = 1)
    
    tpr$cor <- NA
    tpr$pval <- NA
    for(i in (1:nrow(tpr))){
      tpr$cor[i] <- cor( as.numeric(tpr[i,grep("_freq", colnames(tpr))]) , tsr$passage)
      f <- summary(lm( as.numeric(tpr[i,grep("_freq", colnames(tpr))])~tsr$passage))$fstatistic
      pval <- pf(f[1],f[2],f[3],lower.tail=F)
      tpr$pval[i] <- pval
    }
    
    tmp <- tpr[ , c("position", "uid", "cor", "pval")]
    tmp$best <- F
    tmp$best[which.max(abs(tmp$cor))] <- T
    tmp$seriesName <- seriesName
    tmp$varMut <- sd(tsr$rate)
      
    tres <- rbind(tres, tmp)
  }
}

################################################################################
# Variants that are found in multiple series:

tres$OK <- F
tres$multiple <- F
tres$sign <- "N"
tres$sign[which(tres$cor>0)] <- "P"
uids <- unique(tres$uid)
for(uid in uids){
  wu <- which(tres$uid == uid)
  tmp <- tres[wu , ]
  if( min(tmp$pval)<0.05 && length(unique(tmp$sign))==1 ){
    tres$OK[wu] <- T
    if( length(which(tmp$pval<0.05)) > 1 ){
      tres$multiple[wu] <- T
    }
  }
}

to <- as.data.frame(tres[which(tres$OK==T) , ])
to <- to[order(to$position) , ]
to


# Selecting the mutation with the best correlation (for each series)
tb <- tres[which(tres$best == T) , ]
length(which(tb$position < 21560))
# 5 out of 6 series have their most-correlated mutation inside ORF1ab


################################################################################
# Using simulations to test how likely it is to obtain 5 out of 6 in the previous step
NB_RANDOM <- 1000
vSuccess <- rep(NA, NB_RANDOM)

for(STEP in (1:NB_RANDOM)){
  for(Variant in Variants){
    for(Replicate in Replicates){
      seriesName <- paste(Variant, Replicate, sep = "_")
      w <- which(tres$seriesName == seriesName)
      tres$best[w] <- sample(tres$best[w])
    }
  }
  tb <- tres[which(tres$best == T) , ]
  nbSuccess <- length(which(tb$position < 21560))
  vSuccess[STEP] <- nbSuccess
}

pvalue <- length(which(vSuccess>=5)) / length(vSuccess)
pvalue
# --> 0.03
