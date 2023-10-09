################################################################################
# Mutation spectrum as a function of the secondary structure.

MIN_COV <- 1000
MAX_RATE <- 1/MIN_COV
MAX_BEFORE_HS = 1e10

DATA_SET = c("LAN", "CAO")[1]

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


# Some data about the genome structure:
pos_single <- c()

if( DATA_SET == "LAN" ) {
  # This one is a list of single-stranded regions (From: https://www.nature.com/articles/s41467-021-22785-x - SupData 5)
  ts <- read.table("../Structure/2022_NatCom_Lan/structure.ct", sep = "\t", h = F)
  colnames(ts) <- c("position", "base", "prev_index", "next_index", "nbp", "nm")
  pos_single <- ts$position[which(ts$nbp==0)]
  length(pos_single)
}

if( DATA_SET == "CAO" ) {
  # Cao dataset:
  pos_single <- c()
  ts <- read.csv("../Structure/2021_Cao_NatCom/single_stranded_regions.csv", h = T)
  for(i in (1:nrow(ts))) {
    start <- ts$start[i]
    end <- ts$end[i]
    vtmp <- seq(from=start, to=end, by=1)
    pos_single <- c(pos_single, vtmp)
    length(pos_single)
  }
}

# This file contains the unfiltered mutation frequencies for every possible mutation in every sample
fname <- paste(BD_DATA, "all-frequencies-combined-and-lethal.tab.gz", sep="/")
tmuts <- read.table(fname, h=T)

# Extracting the list of sample names from the columns names
vcn <- colnames(tmuts)
wmut <- grepl("_muts", vcn, fixed = TRUE)
wcov <- grepl("_cov", vcn, fixed = TRUE)
vcr <- vcn[which(wmut==T)]
sampleNames <- substr(vcr, 1, nchar(vcr)-5)

# Subsetting only lethal mutations
tl <- tmuts[which(tmuts$lethal == T) , ]

#tl <- tmuts # !!!!! NO LONGER USING THE LETHAL MUTATIONS ONLY

# Adding the single/double strand information:
tl$single_strand = F
tl$single_strand[which(is.element(tl$position, pos_single) == T)] = T


# This is the data.frame that will contain the final result (= spectrum for each sample)
tres <- data.frame(
  sampleName = character(),
  mutation = character(),
  single_strand = logical(),
  nbMut = numeric(),
  nbCalls = numeric(),
  bFrom = character(),
  bTo = character()
)

# Computing the spectrum one sample at a time:
for( sampleName in sampleNames) {
  cat("Working on ", sampleName, " ...\n")
  tmp <- data.frame(
    sampleName = character(),
    mutation = character(),
    single_strand = logical(),
    nbMut = numeric(),
    nbCalls = numeric(),
    bFrom = character(),
    bTo = character()
  )
  
  cn_muts <- paste(sampleName, "muts", sep="_")
  cn_cov <- paste(sampleName, "cov", sep="_")
  tlr <- tl[ , c("position", "uid", "ref", "bTo", "single_strand",
                 cn_muts,
                 cn_cov
                 )]
  
  tlr$freq <- tlr[ , cn_muts] / tlr[ , cn_cov]
  tlr <- tlr[which(tlr[ , cn_cov] >= MIN_COV & tlr$freq <= MAX_RATE & tlr[,cn_muts]<=MAX_BEFORE_HS) , ]
  
  for(bFrom in vBases) {
    for(bTo in vBases[which(vBases != bFrom)]) {
      mutation <- paste(bFrom, "->", bTo, sep = "")
      tbm <- tlr[which(tlr$ref == bFrom & tlr$bTo == bTo) , c("position", cn_muts, cn_cov, "single_strand")]
      
      for(single in c(T,F)){
        
        tbmc <- tbm[which(tbm$single_strand == single) , ]
        
        nbMut <- 0
        cov <- 0
        if( nrow(tbmc) > 0 ){
          cov <- as.numeric(sum(as.numeric(tbmc[ , cn_cov])))
          nbMut <- sum(tbmc[,cn_muts])
        }

        tmpc <- data.frame(
            sampleName = sampleName,
            mutation = mutation,
            single_strand = single,
            nbMut = nbMut,
            nbCalls = cov,
            bFrom = bFrom,
            bTo = bTo
          )
          
          tres <- rbind(tres, tmpc)
      } # Done with single T/F
    } # Done with this bTo
  } # Done with this bFrom
} # Done with this sample

tres$rate <- tres$nbMut / tres$nbCalls
tres$rate[which(tres$nbCalls == 0)] <- 0

tres$single_strand <- factor(tres$single_strand, levels = c(T,F))

tres$ci_min <- 0
tres$ci_max <- 0
for(i in (1:nrow(tres))) {
  pp <- prop.test(tres$nbMut[i], tres$nbCalls[i])
  tres$ci_min[i] <- as.numeric(pp$conf.int[1])
  tres$ci_max[i] <- as.numeric(pp$conf.int[2])
}

fName <- paste("spectra_single_vs_double_strand-", DATA_SET, ".csv", sep = "")
write.csv(tres, file = fName)

fname <- paste("spectra-single-vs-double-strand-", DATA_SET, ".pdf", sep = "")
pdf(file = fname, width = 12, height = 8)

for(sampleName in sampleNames) {
  tr <- tres[which(tres$sampleName == sampleName) , ]

  graphTitle <- paste("Mutation spectrum vs single-strand status - ", sampleName, sep = "")
  
  gp <- ggplot(tr, aes(x = mutation, y = rate, fill = single_strand)) +
    geom_bar(stat="identity", position=position_dodge() ) + 
    geom_errorbar(aes(ymin=ci_min, ymax=ci_max), width=.2, position=position_dodge(.9), color = "grey20") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(graphTitle)
  
  plot(gp)
}
dev.off(dev.cur())

png(file = "C_U_strand.png", width = 800, height = 400)
tr <- tres[which(tres$sampleName == "DELTA_B") , ]
graphTitle <- paste("Mutation spectrum vs single-strand status - ", sampleName, sep = "")
gp <- ggplot(tr, aes(x = mutation, y = rate, fill = single_strand)) +
  geom_bar(stat="identity", position=position_dodge() ) + 
  geom_errorbar(aes(ymin=ci_min, ymax=ci_max), width=.2, position=position_dodge(.9), color = "grey20") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(graphTitle)
plot(gp)
dev.off(dev.cur())

