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


# Build a list of contexts (a list of positions for each possibl triplet)
genomeFile <- paste(DB_DIR, "Sars_cov_2.ASM985889v3.dna.toplevel.fa", sep = "/")
genome <- readDNAStringSet(genomeFile)
ll_triplets <- list()
for(focalBase in vBases) {
  for(b5p in vBases) {
    for(b3p in vBases) {
      triplet <- paste(b5p, focalBase, b3p, sep = "")
      vm <- matchPattern(triplet, genome[[1]])
      ll_triplets[[triplet]] <- start(vm) + 1
    }
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

# Get a list of lethal mutations
lethals <- sars_get_lethal(
  fname_vcf = fname_vcf,
  fname_stop_pos = fname_stop_pos
)


# Add the lethal/non-lethal status for each mutation
tmuts$lethal = F
tmuts$lethal[which( is.element(tmuts$uid, lethals) == T )] = T

# Subsetting only lethal mutations or not
USE_ONLY_LETHALS <- F
tl <- tmuts
if( USE_ONLY_LETHALS == T ){
  tl <- tmuts[which(tmuts$lethal == T) , ]
}


# This is the data.frame that will contain the final result (= spectrum for each sample)
tres <- data.frame(
  sampleName = character(),
  mutation = character(),
  triplet = character(),
  b5p = character(),
  b3p = character(),
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
    b5p = character(),
    b3p = character(),
    nbMut = numeric(),
    nbCalls = numeric(),
    bFrom = character(),
    bTo = character()
  )
  
  cn_muts <- paste(sampleName, "muts", sep="_")
  cn_cov <- paste(sampleName, "cov", sep="_")
  tlr <- tl[ , c("position", "uid", "ref", "bTo", 
                 cn_muts,
                 cn_cov
                 )]
  
  tlr$freq <- tlr[ , cn_muts] / tlr[ , cn_cov]
  tlr <- tlr[which(tlr[ , cn_cov] >= MIN_COV & tlr$freq <= MAX_RATE & tlr[,cn_muts]<=MAX_BEFORE_HS) , ]
  
  for(bFrom in vBases) {
    for(bTo in vBases[which(vBases != bFrom)]) {
      mutation <- paste(bFrom, "->", bTo, sep = "")
      tbm <- tlr[which(tlr$ref == bFrom & tlr$bTo == bTo) , c("position", cn_muts, cn_cov)]
      
      for(b5p in vBases) {
        for(b3p in vBases) {
          triplet <- paste(b5p, bFrom, b3p, sep = "")
          vpos <- ll_triplets[[triplet]]
          tbmc <- tbm[which( is.element(tbm$position, vpos) == T ) , ]
          
          cov <- 0
          nbMut <- 0
          
          if( nrow(tbmc) > 0 ){
            cov <- as.numeric(sum(as.numeric(tbmc[ , cn_cov])))
            nbMut <- sum(tbmc[,cn_muts])
          }
            
          tmpc <- data.frame(
            sampleName = sampleName,
            mutation = mutation,
            triplet = triplet,
            b5p = b5p,
            b3p = b3p,
            nbMut = nbMut,
            nbCalls = cov,
            bFrom = bFrom,
            bTo = bTo
          )
          
          tres <- rbind(tres, tmpc)
        } # Done with this triplet
      } # Done with this b5p
    } # Done with this bTo
  } # Done with this bFrom
} # Done with this sample

tres$rate <- tres$nbMut / tres$nbCalls
tres$rate[which(tres$nbCalls == 0)] <- 0

pdfFileName <- paste("context-heatmaps-all-cov", MIN_COV, ".pdf", sep = "")
if( USE_ONLY_LETHALS == T ){
  pdfFileName <- paste("context-heatmaps-lethal-cov", MIN_COV, ".pdf", sep = "")
}
pdf(file =  pdfFileName, width=12, height=8)
for(bFrom in vBases){
  
  tmp <- tres[which(tres$sampleName == "total" & tres$bFrom == bFrom) , ]
  tmp$context <- paste(tmp$b5p, tmp$b3p, sep = "")
  
  graphTitle = paste("Context dependent rates for ", bFrom, " -> N", sep = "")
  gp <- ggplot(data = tmp, aes(bTo, context, fill= rate)) + 
    geom_tile() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(graphTitle)
  
  plot(gp)
  
}
dev.off(dev.cur())


pdfFileName <- paste("triplets-all-cov", MIN_COV, ".pdf", sep = "")
if( USE_ONLY_LETHALS == T ){
  pdfFileName <- paste("triplets-lethal-cov", MIN_COV, ".pdf", sep = "")
}
pdf(file =  pdfFileName, width=12, height=8)
for(bFrom in vBases){
  
  tmp <- tres[which(tres$sampleName == "total" & tres$bFrom == bFrom) , ]
  tmp$context <- paste(tmp$b5p, tmp$b3p, sep = "")
  
  graphTitle = paste("Context dependent rates for ", bFrom, " -> N", sep = "")
  myDodge = 0.975
  gp <- ggplot(tmp, aes(x = mutation, y = rate, fill=triplet)) +
    geom_col(colour="black",width=0.85, size=0.02, position=position_dodge(myDodge)) +
    geom_text(aes(label = triplet, y = 0), angle = 90, hjust = 1.2, size = 2, position = position_dodge2(myDodge, preserve = "single")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(graphTitle)
  plot(gp)
  
}
dev.off(dev.cur())
