################################################################################
# Looking at mutations for large variation in frequency

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


library("tidyr")
library("directlabels")

MIN_COV <- 100
MIN_RATE <- 0.05
MAX_RATE <- 0.2
vp <- seq(from=1, to=7, by=1)
series <- c("ALPHA_A", "ALPHA_B", "USA_A", "USA_B", "DELTA_A", "DELTA_B")
pdf(file = "high_frequency_variants-v5.pdf", width = 12, height = 8)
for(seriesName in series){
  cn_covs <- paste(seriesName, "_", seq(from=1, to=7), "_cov", sep="")
  tm <- tmuts
  tm$minCov <- apply(tm[ , cn_covs], 1, min)
  tmr <- tm[which(tm$minCov > MIN_COV) , ]
  cn_rates <- paste(seriesName, "_", seq(from=1, to=7), "_rate", sep="")
  tmr$minRate <- apply(tmr[ , cn_rates], 1, min)
  tmr$maxRate <- apply(tmr[ , cn_rates], 1, max)
  tmr$dRate <- tmr$maxRate - tmr$minRate
  
  #tmh <- tmr[which(tmr$minRate > MIN_RATE & tmr$minRate < MAX_RATE) , ]
  #tmh <- tmr[which(tmr$maxRate > MIN_RATE & tmr$minRate < MAX_RATE & tmr$maxRate<0.2) , ]
  #tmh <- tmr[which(tmr$maxRate > 0.02 & tmr$minRate < 0.001 & tmr$maxRate<0.5) , ]
  #tmh <- tmr[which( (tmr$minRate<0.02 & tmr$maxRate>0.05) | (tmr$maxRate>0.98 & tmr$minRate<0.95) ) , ]
  tmh <- tmr[which(tmr$dRate > 0.05) , ]
  
  tmp <- tmh[ , c("uid", cn_rates)]
  colnames(tmp)[2:8] <- seq(from=1, to=7, by=1)
  td <- gather(tmp, key = "passage", value = "frequency", -uid)
  gp <- ggplot(data = td, aes(x = as.numeric(passage), y = frequency, group = uid)) +
      geom_line( aes(color = uid) ) + 
    geom_point( aes(color = uid) ) +
    geom_dl(aes(label = uid, color = uid), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
    scale_x_continuous(breaks=vp, labels=vp, expand=c(0,0.15), limits=c(1,8)) +
    theme(legend.position="none") + 
    ggtitle(seriesName) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Passage")

  plot(gp)
}
dev.off(dev.cur())

plot(gp)
library(plotly)


gp <- ggplot(data = td, aes(x = passage, y = frequency, group = uid)) +
  geom_line( aes(color = uid) ) + 
  geom_point( aes(color = uid) ) +
  theme(legend.position="none") + 
  ggtitle(seriesName) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Passage")

ggplotly(gp)

MIN_COV <- 100
seriesName = "DELTA_B"
cn_covs <- paste(seriesName, "_", seq(from=1, to=7), "_cov", sep="")
tm <- tmuts
tm$minCov <- apply(tm[ , cn_covs], 1, min)
tmr <- tm[which(tm$minCov > MIN_COV) , ]
cn_rates <- paste(seriesName, "_", seq(from=1, to=7), "_rate", sep="")
tmr$minRate <- apply(tmr[ , cn_rates], 1, min)
tmr$maxRate <- apply(tmr[ , cn_rates], 1, max)
tmr$sdRate <- apply(tmr[ , cn_rates], 1, sd)
tmr$maxFC <- tmr$maxRate/tmr$minRate

#tt <- tmr[which(tmr$sdRate>0.000001 & tmr$sdRate<0.01) , ]
tt <- tmr[which(tmr$maxFC > 2 & tmr$minRate > 0.001) , ]

POS_START_SPIKE <- 21563
POS_END_SPIKE <- 25384

plot(tt$maxFC~tt$position, log = "y")
abline(v=POS_START_SPIKE, col="red")
abline(v=POS_END_SPIKE, col="red")


tSpike <- tt[which(tt$position >= POS_START_SPIKE & tt$position <= POS_END_SPIKE) , ]
tRest <- tt[which(tt$position < POS_START_SPIKE | tt$position > POS_END_SPIKE) , ]
tt$inSpike <- F
tt$inSpike[which(tt$position >= POS_START_SPIKE & tt$position <= POS_END_SPIKE)] <- T
boxplot(log(tt$maxFC)~tt$inSpike)
t.test(log(tSpike$maxFC), log(tRest$maxFC))


tmr$ap <- round(tmr$position/1000, 0)
vm <- by(tmr$sdRate, tmr$ap, mean)
tt <- data.frame(
  position = as.numeric(names(vm)),
  sd = as.numeric(vm)
)
plot(tt$sd~tt$position)
