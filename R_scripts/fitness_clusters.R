# Testing the hypothesis that deleterious mutations tend to 
# be clustered in the genome.

setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
library(ggplot2)
library(tidyverse)

MAX_DISTANCE <- 4
MAX_FITNESS <- 0

tf <- read.csv("fitness-mapped-to-protein.csv")
tf$w[which(tf$w<0)] <- 0
tf$w[which(tf$w>1.5)] <- 1.5

tf$puid <- paste(tf$Name, tf$proteinLoc, sep = "_")
tfr <- tf[ , c("puid", "w")]
tfr$nb <- 1
vnb <- by(tfr$nb, tfr$puid, sum)

vmin <- by(tfr$w, tfr$puid, min)
tmin <- data.frame(puid = names(vmin), w = as.numeric(vmin))
tmin$geneName <- unlist(strsplit(tmin$puid, "_", fixed = T))[seq(from=1, to=2*nrow(tmin)-1, by=2)]
tmin$proteinLoc <- as.numeric(unlist(strsplit(tmin$puid, "_", fixed = T))[seq(from=2, to=2*nrow(tmin), by=2)])


################################################################################
# getNbInClusters
# Finds the number of mutations that have a fitness value lower than "maxW"
# and are no more than "maxDistance" away from each other.
getNbInClusters <- function(tf, maxDistance = 10, maxW = 0){
  nbInClusters <- 0
  tfr <- tf[which(tf$w <= maxW) , ]
  for(geneName in unique(tfr$geneName)){
    tfg <- tfr[which(tfr$geneName == geneName) , ]
    nb <- nrow(tfg)
    if( nb > 1 ){
      tfg <- tfg[order(tfg$proteinLoc, decreasing = F) , ]
      tfg$dist <- NA
      tfg$dist[1:(nb-1)] <- tfg$proteinLoc[2:nb] - tfg$proteinLoc[1:(nb-1)]
      nbInClusters <- nbInClusters + length(which(tfg$dist<=maxDistance))
    }
  }
  nbInClusters
}

nbInClustersReal <- getNbInClusters(tf = tmin, maxDistance = MAX_DISTANCE, maxW = MAX_FITNESS)

# Randomizing the fitness values to test of clustering decreases
trand <- tmin
NB_RAND <- 1000
vRand <- rep(NA, NB_RAND)
for(iRand in (1:NB_RAND)){
  trand$w <- sample(trand$w)
  vRand[iRand] <- getNbInClusters(tf = trand, maxDistance = MAX_DISTANCE, maxW = MAX_FITNESS)
}

length(which(vRand > nbInClustersReal))
hist(vRand, xlim = c(min(vRand), max(max(vRand, nbInClustersReal))))
rug(x = nbInClustersReal, col = "red")

# --> Clear pattern of clusterisation


################################################################################
# Looking at it the other way: distribution of difference in fitness values 
# between consecutive mutations that are near each other vs far apart.

MAX_DISTANCE <- 4
MIN_DISTANCE <- 10

vw_near <- c()
vw_far <- c()
tfr <- tf[which(tf$consequence == "nonsynonymous") , ]
colnames(tfr)[which(colnames(tfr)=="Name")] <- "geneName"
for(geneName in unique(tfr$geneName)){
  tfg <- tfr[which(tfr$geneName == geneName) , ]
  nb <- nrow(tfg)
  if( nb > 1 ){
    tfg <- tfg[order(tfg$proteinLoc, decreasing = F) , ]
    tfg$dist <- NA
    tfg$dist[1:(nb-1)] <- tfg$proteinLoc[2:nb] - tfg$proteinLoc[1:(nb-1)]
    tfg$deltaw <- NA
    tfg$deltaw[1:(nb-1)] <- abs(tfg$w[2:nb] - tfg$w[1:(nb-1)])
    vw_near <- c(vw_near, tfg$deltaw[which(tfg$dist <= MAX_DISTANCE)])
    vw_far <- c(vw_far, tfg$deltaw[which(tfg$dist > MIN_DISTANCE)])
  }
}

length(vw_near)
length(vw_far)
t.test(vw_near, vw_far)
