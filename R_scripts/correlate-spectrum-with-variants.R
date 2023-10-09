################################################################################
#
# This script searches for variants whose frequency changes correlate with 
# changes in the mutation rate across the passages of each series.


####################################################################################################
# First part: location of files:
setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
source("../R_files/seqan-functions.R")
source("../R_files/covid-functions.R")

FILE_PROTEIN_PREDICT <- "../DB/protein_predict.rds" # R binary file with the results from predicting mutation consequence (synonymous, non-synonymous, ...)
FILE_ALL_FREQUENCIES <- "all-frequencies.rds"
FILE_SPECTRUM <- "spectrum-lethal-COV-1000_RATE-0.001_HS-1000.csv"
FILE_NSPS <- "misc/nsp.csv"

ORF1ab_END <- 21560
EXCLUDE_C_TO_T <- T

library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(genomation)
library(Biostrings)
library(tidyr)
library(readxl)

vBases = getBases()
vMuts = getMutsFromBases(vBases)
if( EXCLUDE_C_TO_T == T ){ 
  vMuts <- vMuts[which(vMuts != "C->T")]
  FILE_SPECTRUM <- "spectrum-lethal-no_C_to_T-COV-1000_RATE-0.001_HS-1000.csv"
}

################################################################################
# This file contains the frequencies of every mutation in every sample
tAll <- readRDS(FILE_ALL_FREQUENCIES)
cNames <- colnames(tAll)
cnCov <- cNames[grep("_cov$", colnames(tAll) )]
sampleNames <- substr(cnCov, 1, nchar(cnCov)-4)
for(sampleName in sampleNames){
  cn <- paste(sampleName, "_freq", sep = "")
  cnCov <- paste(sampleName, "cov", sep = "_")
  cnMut <- paste(sampleName, "muts", sep = "_")
  tAll[ , cn] <- tAll[ , cnMut] / tAll[ , cnCov]
}

# Keeping only nonsynonymous (and nonsense) variants:
res <- readRDS(FILE_PROTEIN_PREDICT)
res$uid <- paste( start(res), as.character(res$REF), as.character(res$varAllele), sep = "_" )
uid_ns <- res$uid[which(res$CONSEQUENCE!="synonymous")]
tAll <- tAll[which( is.element(tAll$uid, uid_ns)==T) , ]

################################################################################
# Reading the mutation spectrum for all the samples
ts <- read.csv(FILE_SPECTRUM)

Variants <- c("ALPHA", "USA", "DELTA")
Replicates <- c("A", "B")


################################################################################
# Keeping only positions that are covered by at least 10 reads, have a relative 
# change in frequency of at least 0.1 and reach an absolute frequency of at 
# least 5% in one of the passages.
MIN_COV <- 100
minDelta <- 0.01
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
  #for(Replicate in Replicates){
    #seriesName <- paste(Variant, Replicate, sep = "_")
    seriesName <- Variant
    
    tpos <- getChangingAlleleFrequencies(tAll = tAll, seriesName = seriesName, 
                                         MIN_COV = MIN_COV, minDelta = minDelta, minAbsolute = minAbsolute)
    
    
    sPattern <- paste(seriesName, ".*_freq", sep = "")
    tpr <- tpos[ , c(1:4, grep(sPattern, colnames(tpos)))]
    
    #tsr <- ts[which(ts$variant==Variant & ts$replicate==Replicate & is.na(ts$passage)==F & ts$mutation=="N->N") , c("passage", "rate")]
    tsr <- ts[which(ts$variant==Variant & is.na(ts$passage)==F & ts$mutation=="N->N") , c("passage", "rate")]
    
    cnFreq <- grep(sPattern, colnames(tpos))
    tpr <- tpos[ , c(1, 2, cnFreq)]
    #colnames(tpr)[2:ncol(tpr)] <- seq(from=1, to=(ncol(tpr)-1), by = 1)
    
    tpr$cor <- NA
    tpr$pval <- NA
    for(i in (1:nrow(tpr))){
      #tpr$cor[i] <- cor( as.numeric(tpr[i,grep("_freq", colnames(tpr))]) , tsr$passage)
      tpr$cor[i] <- cor( as.numeric(tpr[i,grep("_freq", colnames(tpr))]) , tsr$rate)
      #f <- summary(lm( as.numeric(tpr[i,grep("_freq", colnames(tpr))])~tsr$passage))$fstatistic
      f <- summary(lm( as.numeric(tpr[i,grep("_freq", colnames(tpr))])~tsr$rate))$fstatistic
      pval <- pf(f[1],f[2],f[3],lower.tail=F)
      tpr$pval[i] <- pval
    }
    
    tmp <- tpr[ , c("position", "uid", "cor", "pval")]
    tmp$best <- F
    tmp$best[which.max(abs(tmp$cor))] <- T
    tmp$seriesName <- seriesName
    tmp$varMut <- sd(tsr$rate)
      
    tres <- rbind(tres, tmp)
#  }
}

tres$adj <- p.adjust(tres$pval, method = "BH")

# Selecting the mutation with the best correlation (for each series)
#tb <- tres[which(tres$best == T) , ]
#length(which(tb$position < ORF1ab_END))


# Generating a table with all mutations inside ORF1ab that have a significant correlation with mutation rate
t <- tres[which(tres$position < ORF1ab_END & tres$pval < 0.05) , ]
t <- tres[which(tres$pval < 0.05) , ]
#t <- tres[which(tres$adj < 0.05) , ]
dim(t)

vnb <- by(t$pval, t$uid, min)
tnb <- as_tibble( data.frame(uid = as.character(names(vnb)), 
                  pval = as.numeric(vnb)
    ))
stmp <- unlist(strsplit(tnb$uid, "_"))
tnb$position <- as.numeric(stmp[seq(from=1, to=length(stmp)-2, by = 3)])
tnb$bFrom <- stmp[seq(from=2, to=length(stmp)-1, by = 3)]
tnb$bTo <- stmp[seq(from=3, to=length(stmp), by = 3)]
tnb$cor <- NA
tnb$series <- NA
tnb$geneName <- NA
tnb$proteinLoc <- NA
tnb$refAA <- NA
tnb$varAA <- NA
for(i in (1:nrow(tnb))){
  uid <- tnb$uid[i]
  w <- which(t$uid == uid)
  tnb$series[i] <- paste(t$seriesName[w], round(t$cor[w],2), collapse = ",")
  
  vCor <- t$cor[w]
  tnb$cor[i] <- vCor[which.max(abs(vCor))]
  
  pos <- tnb$position[i]
  bFrom <- tnb$bFrom[i]
  bTo <- tnb$bTo[i]
  w <- which(start(res)==pos & res$REF==bFrom & res$varAllele==bTo)
  if( length(w) > 0 ){
    tnb$refAA[i] <- paste(res$REFAA[w], collapse = ",")
    tnb$varAA[i] <- paste(res$VARAA[w], collapse = ",")
    tnb$geneName[i] <- unique(res$Name[w])
    tnb$proteinLoc[i] <- unique(as.numeric(res$PROTEINLOC[w]))
  }
}

tnsp <- read.csv(FILE_NSPS, h = T)
#tnb$nsp <- NA
#tnb$protLoc <- NA
for(i in (1:nrow(tnb))){
  pos <- tnb$position[i]
  w <- which(tnsp$start<= pos & tnsp$end>=pos)
  if( length(w) == 1 ){
    #tnb$nsp[i] <- tnsp$nsp[w]
    tnb$geneName[i] <- tnsp$nsp[w]
    tnb$proteinLoc[i] <- ceiling(( (pos - tnsp$start[w]) + 1 ) / 3)
  } else {
    cat("NOT INSIDE NSP.\n")
  }
}


tFinal <- tnb[ , c("position", "bFrom", "bTo", "cor", "pval", "geneName", "proteinLoc", "refAA", "varAA", "series")]
tFinal$mutation <- paste(tFinal$geneName, ":", tFinal$refAA, tFinal$proteinLoc, tFinal$varAA, sep = "")
tfo <- tFinal[order(tFinal$position, decreasing = F) , ]
write.csv(tfo, file = "SuppTable_variants_fidelity-real.csv", row.names = F)


tmp <- tfo[which(tfo$geneName != "S") , ]
prop.test( length(which(tmp$cor<0)), nrow(tmp), p = 0.5 )
length(which(tmp$cor<0))
nrow(tmp)
