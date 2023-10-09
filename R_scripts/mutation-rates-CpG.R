################################################################################
# Mutation rate at CpG vs non CpG sites.

MIN_COV <- 1000
MAX_RATE <- 1/MIN_COV
MAX_BEFORE_HS = 1e10

BD_DATA <- ""
DB_DIR <- ""

if( .Platform$OS.type == "windows" ){
  setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
  source("C:/lab/TR_ERRORS/R_files/seqan-test/seqan-functions-v3.R")
  source("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/R_files/getLethals.R")
  source("../R_files/covid-functions.R")
  
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

genome <- readDNAStringSet("../DB/ensembl_ASM985889v3/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
positions_CpG <- start(matchPattern("CG", genome[[1]]))

tg <- unique(tmuts[ , c("position", "ref")])
rownames(tg) <- as.numeric(tg$position)
tg$CpG <- F
w <- which(tg$ref == "C")
ww <- which(tg$ref[w+1]=="G")
tg$CpG[w[ww]] <- T
#tt <- tg[which(tg$CpG == T) , ]
positions_CpG <- tg$position[which(tg$CpG==T)]


# Extracting the list of sample names from the columns names
vcn <- colnames(tmuts)
wmut <- grepl("_muts", vcn, fixed = TRUE)
wcov <- grepl("_cov", vcn, fixed = TRUE)
vcr <- vcn[which(wmut==T)]
sampleNames <- substr(vcr, 1, nchar(vcr)-5)

# Subsetting only lethal mutations
tl <- tmuts[which(tmuts$lethal == T) , ]
tl <- tmuts
tl$CpG <- F
tl$CpG[which(is.element(tl$position, positions_CpG)==T)] = T

# This is the data.frame that will contain the final result (= spectrum for each sample)
tres <- data.frame(
  sampleName = character(),
  mutation = character(),
  CpG = logical(),
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
    singl_strand = logical(),
    nbMut = numeric(),
    nbCalls = numeric(),
    bFrom = character(),
    bTo = character()
  )
  
  cn_muts <- paste(sampleName, "muts", sep="_")
  cn_cov <- paste(sampleName, "cov", sep="_")
  tlr <- tl[ , c("position", "uid", "ref", "bTo", "CpG",
                 cn_muts,
                 cn_cov
                 )]
  
  tlr$freq <- tlr[ , cn_muts] / tlr[ , cn_cov]
  tlr <- tlr[which(tlr[ , cn_cov] >= MIN_COV & tlr$freq <= MAX_RATE & tlr[,cn_muts]<=MAX_BEFORE_HS) , ]
  
  for(bFrom in vBases) {
    for(bTo in vBases[which(vBases != bFrom)]) {
      mutation <- paste(bFrom, "->", bTo, sep = "")
      tbm <- tlr[which(tlr$ref == bFrom & tlr$bTo == bTo) , c("position", cn_muts, cn_cov, "CpG")]
      
      for(CpG in c(T,F)){
        
        tbmc <- tbm[which(tbm$CpG == CpG) , ]
        
        nbMut <- 0
        cov <- 0
        if( nrow(tbmc) > 0 ){
          cov <- as.numeric(sum(as.numeric(tbmc[ , cn_cov])))
          nbMut <- sum(tbmc[,cn_muts])
        }

        tmpc <- data.frame(
            sampleName = sampleName,
            mutation = mutation,
            CpG = CpG,
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

tres <- tres[which(tres$bFrom == "C") , ]
tres$mutation <- factor(tres$mutation, levels = c("C->T", "C->A", "C->G"))
tres$rate <- tres$nbMut / tres$nbCalls
tres$rate[which(tres$nbCalls == 0)] <- 0

tres$CpG <- factor(tres$CpG, levels = c(F, T))

tres$ci_min <- 0
tres$ci_max <- 0
for(i in (1:nrow(tres))) {
  pp <- prop.test(tres$nbMut[i], tres$nbCalls[i])
  tres$ci_min[i] <- as.numeric(pp$conf.int[1])
  tres$ci_max[i] <- as.numeric(pp$conf.int[2])
}

#pdf(file = "spectra-single-vs-double-strand-Lan.pdf", width = 12, height = 8)
pdf(file = "spectrum-CpG.pdf", width = 12, height = 8)

for(sampleName in sampleNames) {
  tr <- tres[which(tres$sampleName == sampleName) , ]

  graphTitle <- paste("C->N mutation rates in CpG vs outside CpG - ", sampleName, sep = "")
  
  gp <- ggplot(tr, aes(x = mutation, y = rate, fill = CpG)) +
    geom_bar(stat="identity", position=position_dodge() ) + 
    geom_errorbar(aes(ymin=ci_min, ymax=ci_max), width=.2, position=position_dodge(.9), color = "grey20") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(graphTitle)
  
  plot(gp)
}
dev.off(dev.cur())



png(file = "CpG.png", width = 800, height = 400)
tr <- tres[which(tres$sampleName == "DELTA_B") , ]
graphTitle <- paste("C->N mutation rates in CpG vs outside CpG - ", sampleName, sep = "")
gp <- ggplot(tr, aes(x = mutation, y = rate, fill = CpG)) +
  geom_bar(stat="identity", position=position_dodge() ) + 
  geom_errorbar(aes(ymin=ci_min, ymax=ci_max), width=.2, position=position_dodge(.9), color = "grey20") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(graphTitle)
plot(gp)

dev.off(dev.cur())
