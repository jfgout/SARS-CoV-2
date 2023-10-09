################################################################################
# Estimating the nucleotic content under a purely mutational process:

MIN_COV <- 1000
MAX_RATE <- 1/MIN_COV
MAX_BEFORE_HS = 1e10

USE_ONLY_LETHAL <- T

setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")

source("../R_files/covid-functions.R")
source("../R_files/seqan-functions.R")

library(Biostrings)
library(ggplot2)
library(tidyr)


vBases = getBases()
vMuts = getMutsFromBases(vBases)


# A quick way to get the list of nucleotides in the genome
genome <- readDNAStringSet("../DB/ensembl_ASM985889v3/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
seq <- as.character(genome[[1]])

vGenome <- strsplit(seq, split = "")

fname <- paste("spectrum-COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
tres <- read.csv(fname)
# Using data from the DELTA_B series because that's the one we sequences ultra deep
spectrum <- tres[which(tres$series == "DELTA_B" & is.na(tres$passage)==T & tres$lethal==T & tres$mutation!="N->N") , ]

spectrum$rate <- spectrum$rate * 10

NB_GENERATIONS <- 4000

composition <- data.frame(
  A = rep(NA, NB_GENERATIONS),
  T = rep(NA, NB_GENERATIONS),
  C = rep(NA, NB_GENERATIONS),
  G = rep(NA, NB_GENERATIONS)
  )

getATCG <- function(vGenome, vBases){
  nbBases <- length(vBases)
  genomeSize <- length(vGenome)
  vRes <- rep(0, nbBases)
  for(i in (1:nbBases)){
    vRes[i] <- length(which(vGenome == vBases[i])) / genomeSize
  }
  vRes
}
composition[1,] <- getATCG(vGenome = vGenome, vBases = vBases)

mutateGenome <- function(vGenome, spectrum, vBases){
  for(bFrom in vBases){
    wb <- which(vGenome == bFrom)
    for(bTo in vBases[which(vBases != bFrom)]){
      rate <- spectrum$rate[which(spectrum$bFrom==bFrom & spectrum$bTo==bTo)]
      rb <- rbinom(n=length(wb), size = 1, prob = rate)
      wmutated <- which(rb==1)
      if( length(wmutated) >= 1 ){
        vGenome[wb[wmutated]] <- bTo
      }
    }
  }

  vGenome
}

for(i in (2:NB_GENERATIONS)){
  vGenome <- mutateGenome(vGenome = vGenome, spectrum = spectrum, vBases = vBases)
  composition[i,] <- getATCG(vGenome = vGenome, vBases = vBases)
  if( i%%1000 == 0){
    print(i)
  }
}

composition$generation <- seq(from = 1, to = nrow(composition), by = 1)

tt <- as.data.frame(pivot_longer(data = composition, names_to = "nucleotide", values_to = "frequency", cols = vBases))
gp <- ggplot(tt, aes(x = generation, y=frequency, color = nucleotide)) +
    geom_point()
plot(gp)

