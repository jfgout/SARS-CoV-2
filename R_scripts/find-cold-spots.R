################################################################################
# Looking for locations with abnomrlally low mutation rate

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


library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(genomation)
library(Biostrings)

library(readxl)

vBases = getBases()
vMuts = getMutsFromBases(vBases)

# This file contains the unfiltered mutation frequencies for every possible mutation in every sample
fname <- paste(BD_DATA, "all-frequencies-combined-and-lethal.tab.gz", sep="/")
tmuts <- read.table(fname, h=T)

# Extracting the list of sample names from the columns names
vcn <- colnames(tmuts)
wmut <- grepl("_muts", vcn, fixed = TRUE)
vcr <- vcn[which(wmut==T)]
sampleNames <- substr(vcr, 1, nchar(vcr)-5)
sampleNames <- sampleNames[which(sampleNames != "total")]

for(sampleName in sampleNames){
  cn_rate <- paste(sampleName, "rate", sep="_")
  cn_err <- paste(sampleName, "muts", sep = "_")
  cn_obs <- paste(sampleName, "cov", sep = "_")
  tmuts[ , cn_rate] <- tmuts[ , cn_err] / tmuts[ , cn_obs]
}


tSpec <- read.csv("spectrum-lethal_COV-1000_RATE-0.001_HS-1e+10.csv")
MIN_COV <- 1000

ll_cold_spots <- list()
series <- c("ALPHA_A", "ALPHA_B", "USA_A", "USA_B", "DELTA_A", "DELTA_B")
for(series_name in series){
  cat("Processing ", series_name, "\n")
  series_max <- 7
  series_cov <- paste(series_name, "_", seq(from=1,to=series_max,by=1), "_cov", sep = "")
  series_mut <- paste(series_name, "_", seq(from=1,to=series_max,by=1), "_muts", sep = "")
  tmuts$mut <- apply(tmuts[ , series_mut], 1, sum)
  tmuts$cov <- apply(tmuts[ , series_cov], 1, sum)
  tmuts$rate <- tmuts$mut / tmuts$cov

  tmr <- tmuts[which(tmuts$cov>MIN_COV & tmuts$rate<1/MIN_COV) , c("position", "uid", "ref", "bTo", "lethal", "mut", "cov", "rate")]
  tsr <- tSpec[which(tSpec$sampleName == series_name) , ]
  tmr <- mutation_rate_deviation(tsr = tsr, tmr = tmr)
  tmr$padj <- p.adjust(tmr$pval, method = "fdr")
  ll_cold_spots[[series_name]] <- tmr
}

save(ll_cold_spots, file = "ll_cold_spots.bin")

MAX_FC <- 0.01
MAX_PVAL <- 0.0001
v_cold_spots <- c()

for(series_name in series){
  tmrs <- ll_cold_spots[[series_name]]
  tms <- tmrs[which(tmrs$FC < MAX_FC & tmrs$padj < MAX_PVAL) , ]
  v_cold_spots <- c(v_cold_spots, tms$uid)
  cat(series_name, ": ", nrow(tms), "\n")
}

length(v_cold_spots)
vcsu <- unique(v_cold_spots)
length(unique(vcsu))

tcs <- data.frame(
  uid = character(),
  position = numeric(),
  bFrom = character(),
  bTo = character(),
  nbObserved = numeric(),
  nbExpected = numeric(),
  padj = numeric(),
  minDecrease = numeric(),
  seriesName = character()
)

for(series_name in series){
  tmrs <- ll_cold_spots[[series_name]]
  tms <- tmrs[which(is.element(tmrs$uid, vcsu) == T) , 
              c("uid", "position", "ref", "bTo", "mut", "nbExpected", "padj") ]
  
  tms$minDecrease <- tms$nbExpected / tms$mut
  wzero <- which(tms$mut == 0)
  tms$minDecrease[wzero] <- tms$nbExpected[wzero]
  tms$seriesName <- series_name
  colnames(tms) <- colnames(tcs)
  tcs <- rbind(tcs, tms)
  
  ttc <- tms
  ttc$FC <- ttc$nbObserved / ttc$nbExpected
  frac1 <- length(which(tmrs$nbExpected > 1 & tmrs$FC < 1)) / length(which(tmrs$nbExpected > 1))
  frac2 <- length(which(ttc$nbExpected > 1 & ttc$FC < 1)) / length(which(ttc$nbExpected > 1))
  cat(series_name, ": ", frac1, " vs ", frac2, "\n")
}

write.csv(tcs, file = "cold-spots-all-series.csv")

MAX_FC <- 0.1
MAX_ADJ_PVAL <- 0.001
tmr <- ll_cold_spots[["DELTA_B"]]
tmr$mutation <- paste(tmr$ref, "->", tmr$bTo, sep = "")
tms <- tmr[which(tmr$FC<MAX_FC & tmr$padj<MAX_ADJ_PVAL) , ]
dim(tms)
max(tms$FC)

# Loading the information about synonymous vs non-synonymous mutations
load("../DB/protein_predict.bin")
res$uid <- paste( start(res), as.character(res$REF), as.character(res$varAllele), sep = "_" )
tr <- data.frame(uid = res$uid,
                 consequence = res$CONSEQUENCE
)
tmr <- merge(tmr, tr, by = "uid")
tms <- tmr[which(tmr$FC<MAX_FC & tmr$padj<MAX_ADJ_PVAL) , ]


tsr$nbTot <- 0 # Total number of positions with this type of mutation
tsr$nbCold <- 0 # Total number of cold spots with this type of mutation
tsr$nbLethal_all <- 0 # Number of lethal mutations of this type
tsr$nbLethal_cold <- 0 # Number of cold lethal mutations among cold spots of this type
tsr$nbSynAll <- 0 # Number of synonymous mutations of this type
tsr$nbSynCold <- 0 # Number of synonymous cold spots
tsr$nbNonsenseAll <- 0
tsr$nbNonsenseCold <- 0

tsr$lethalEnrichement <- 0 # Excess of cold spots among lethal mutations
tsr$lethalPval <- 0
tsr$nonsenseEnrichement <- 0
tsr$nonsensePval <- 0
tsr$synDepletion <- 0
tsr$synPval <- 0

for(mutation in unique(tms$mutation)){
  wm <- which(tsr$mutation == mutation)
  tsr$nbTot[wm] <- length(which(tmr$mutation == mutation))
  tsr$nbCold[wm] <- length(which(tms$mutation == mutation))
  tsr$nbLethal_all[wm] <- length(which(tmr$mutation == mutation & tmr$lethal == T))
  tsr$nbLethal_cold[wm] <- length(which(tms$mutation == mutation & tms$lethal == T))
  
  tsr$nbSynAll[wm] <- length(which(tmr$mutation == mutation & tmr$consequence == "synonymous"))
  tsr$nbSynCold[wm] <- length(which(tms$mutation == mutation & tms$consequence == "synonymous"))
  tsr$nbNonsenseAll[wm] <- length(which(tmr$mutation == mutation & tmr$consequence == "nonsense"))
  tsr$nbNonsenseCold[wm] <- length(which(tms$mutation == mutation & tms$consequence == "nonsense"))
  
  pp <- prop.test( c(tsr$nbLethal_all[wm], tsr$nbLethal_cold[wm]), 
                          c(tsr$nbTot[wm], tsr$nbCold[wm])
                          )
  tsr$lethalEnrichement[wm] <- as.numeric(pp$estimate[2]/pp$estimate[1])
  tsr$lethalPval[wm] <- as.numeric(pp$p.value)
  
  pp <- prop.test( c(tsr$nbNonsenseAll[wm], tsr$nbNonsenseCold[wm]),
                  c(tsr$nbTot[wm], tsr$nbCold[wm])
                  )
  tsr$nonsenseEnrichement[wm] <- as.numeric(pp$estimate[2]/pp$estimate[1])
  tsr$nonsensePval[wm] <- as.numeric(pp$p.value)

  pp <- prop.test( c(tsr$nbSynAll[wm], tsr$nbSynCold[wm]),
                   c(tsr$nbTot[wm], tsr$nbCold[wm])
  )
  tsr$synDepletion[wm] <- as.numeric(pp$estimate[1]/pp$estimate[2])
  tsr$synPval[wm] <- as.numeric(pp$p.value)
  
}

write.csv(tsr, file = "cold-spots-enrichements.csv", row.names = F)

tsr[ , c("mutation", "lethalEnrichement", "nonsenseEnrichement", "synDepletion")]

