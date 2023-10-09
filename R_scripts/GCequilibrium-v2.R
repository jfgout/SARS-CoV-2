################################################################################
# Estimating the nucleotide content under a purely mutational process:

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


fname <- paste("spectrum-COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
tres <- read.csv(fname)

# Using data from the DELTA_B series because that's the one we sequences ultra deep
spectrum <- tres[which(tres$series == "DELTA_B" & is.na(tres$passage)==T & tres$lethal==T & tres$mutation!="N->N") , ]

# Using the combined rate from all the series:
spectrum <- tres[which(tres$series == "total" & tres$lethal==T & tres$mutation!="N->N") , ]

spectrum$rate <- spectrum$rate * 20

NB_GENERATIONS <- 5000

composition <- data.frame(
  A = rep(NA, NB_GENERATIONS),
  C = rep(NA, NB_GENERATIONS),
  G = rep(NA, NB_GENERATIONS),
  T = rep(NA, NB_GENERATIONS)
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

af <- alphabetFrequency(genome[[1]])[1:4]
composition[1,] <- af

mutateGenome <- function(af, spectrum, vBases){
  for(bFrom in vBases){
    for(bTo in vBases[which(vBases != bFrom)]){
      rate <- spectrum$rate[which(spectrum$bFrom==bFrom & spectrum$bTo==bTo)]
      rb <- rbinom(n=af[bFrom], size = 1, prob = rate)
      nb <- length(which(rb>0))
      if( nb >= 1 ){
        af[bFrom] <- af[bFrom] - nb
        af[bTo] <- af[bTo] + nb
      }
    }
  }
  af
}

for(i in (2:NB_GENERATIONS)){
  af <- mutateGenome(af = af, spectrum = spectrum, vBases = vBases)
  composition[i,] <- af
  if( i%%1000 == 0){
    print(i)
  }
}

composition$generation <- seq(from = 1, to = nrow(composition), by = 1)
composition$total <- composition$A + composition$C + composition$G + composition$T
#composition$GC <- (composition$G+composition$C) / composition$total
for(base in vBases){
  composition[ , base] <- composition[ , base] / composition$total
}

cpr <- composition[seq(from=1, to=NB_GENERATIONS, by=100) , ]

tt <- as.data.frame(pivot_longer(data = cpr, names_to = "nucleotide", values_to = "frequency", cols = vBases))
gp <- ggplot(tt, aes(x = generation, y=frequency, color = nucleotide)) +
    geom_point()
plot(gp)

library(tidyr)
FILE_PROTEIN_PREDICT <- "../DB/protein_predict.rds"
res <- readRDS(FILE_PROTEIN_PREDICT)

res$numc <- 1
res$numc[which(res$CONSEQUENCE=="synonymous")] <- 0
res$pos <- start(res)
vc <- by(res$numc, res$pos, sum)

fourFoldSites <- as.numeric(names(vc[which(vc==0)]))
fourFoldSequences <- genome[[1]][fourFoldSites]
af4 <- alphabetFrequency(fourFoldSequences)[1:4]
af <- alphabetFrequency(genome[[1]])[1:4]

af4/sum(af4)
af/sum(af)


