#######################################################################################
# Estimating fitness values for each mutation + amount of variation between replicates

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

MIN_COV <- (-1)
MIN_RATE <- (-0.1)
MAX_RATE <- 1.1
vp <- seq(from=1, to=7, by=1)
series <- c("ALPHA_A", "ALPHA_B", "USA_A", "USA_B", "DELTA_A", "DELTA_B")

seriesName <- "DELTA_B"
NB_PASSAGES <- 7
cnMuts <- paste(seriesName, seq(from = 1, to = NB_PASSAGES, by = 1), "muts", sep = "_")
cnObs <- paste(seriesName, seq(from = 1, to = NB_PASSAGES, by = 1), "cov", sep = "_")

#MIN_RATE <- 0.0001
tm <- tmuts
seriesSampleNames <- paste(seriesName, "_", seq(from = 1, to = NB_PASSAGES, by = 1), sep = "")
rates <- paste(seriesSampleNames, "_rate", sep = "")
tm$maxRate <- apply(tm[ , rates], 1, max)
tmr <- tm[which(tm$maxRate > MIN_RATE) , ]
dim(tmr)


get_w <- function(vv, u, NB_PASSAGES){
  vMuts <- vv[1:NB_PASSAGES]
  vObs <- vv[(NB_PASSAGES+1):(NB_PASSAGES*2)]
  nb <- length(vObs)
  vFreqs <- vMuts / vObs
  
  vw <- (vFreqs[2:nb]-u) / vFreqs[1:(nb-1)]
  vw[which(vw==Inf)] <- 2
  vw[which(vw==(-Inf))] <- 0
  cc <- cor(as.numeric(vw), (1:6))
  vWeights <- rep(0, nb-1)
  for(i in (2:nb)){
      vWeights[i-1] <- min(vObs[(i-1):i]) * u
  }
  w <- sum(vw*vWeights) / sum(vWeights)
  
  mw <- max(vWeights)
  
  #wmax <- which.max(vWeights)
  #w <- as.numeric(vw[wmax])
  
  
  res <- c(w = w, mw = mw, cc = cc)
  res
}

get_w_var <- function(vRates, nbPassages){
  vDif <- vRates[2:length(nbPassages)] / vRates[1:(length(nbPassages)-1)]
  wneg <- which(vDif < 1.0)
  if( length(wneg) > 0 ){
    vDif[wneg] <- 1/vDif[wneg]
  }
  w <- sum(vDif)
}

tmr$w <- NA
#tmr$w <- apply(tmr[ , rates], 1, get_w, NB_PASSAGES)
#tmr$w <- apply(tmr[ , c(cn_err, cn_obs)], 1, get_w, NB_PASSAGES)


tSpec <- read.csv("spectrum-lethal_COV-1000_RATE-0.001_HS-1e+10.csv")
tmr$mut <- paste(tmr$ref, "->", tmr$bTo, sep = "")
allMutations <- getMutsFromBases(vBases)


################################################################################
# Computing the fitness for each of the 12 possible type of base-substitution 
# independently. That's because u is a parameter to pass to the function that 
# computes the fitness. This avoids having to find the value of u for each new 
# entry.
tmr$w <- NA
tmr$mw <- NA
tmr$cc <- NA
for(mut in allMutations){
  u <- tSpec$rate[which(tSpec$sampleName == seriesName & tSpec$mutation == mut)]
  wm <- which(tmr$mut == mut)
  #tmr[wm,w] <- 
   ares <- apply(tmr[wm , c(cnMuts, cnObs)], 1, get_w, u, NB_PASSAGES)
   tmr$w[wm] <- ares["w" , ]
   tmr$mw[wm] <- ares["mw" , ]
   tmr$cc[wm] <- ares["cc" , ]
}

MIN_WEIGHT <- 3
#MIN_COR <- 0.5
tmrn <- tmr[which(is.nan(tmr$w) == F & tmr$mw > MIN_WEIGHT ) , ]
#tmrn <- tmr[which(is.nan(tmr$w) == F & tmr$mw > MIN_WEIGHT & tmr$cc > MIN_COR) , ]
dim(tmrn)
boxplot(tmrn$w~tmrn$lethal)
t.test(tmrn$w[which(tmrn$lethal==T)], tmrn$w[which(tmrn$lethal==F)])

MAX_W <- 1.5
tmp <- tmrn
length(which(tmp$lethal==T))
tmp$w[which(tmp$w<0)] <- 0
tmp$w[which(tmp$w>MAX_W)] <- MAX_W
gp <- ggplot(tmp, aes(x = w, fill = lethal, colour = lethal)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 20)
plot(gp)

library(tidyverse)
res <- readRDS("C:/lab/TR_ERRORS/DATA/Coronavirus/DB/protein_predict.rds")
#res <- readRDS("C:/lab/TR_ERRORS/DATA/Coronavirus/DB/protein_predict-ORF1ab.rds")

tc <- data.frame(
  uid = as.character(res$uid), 
  geneID = as.character(res$GENEID), 
  Name = res$Name,
  proteinLoc = as.numeric(res$PROTEINLOC), 
  refAA = as.character(res$REFAA),
  varAA = as.character(res$VARAA),
  consequence = as.character(res$CONSEQUENCE)
  )


tmrc <- merge(tmrn, tc, by = "uid")
boxplot(tmrc$w~tmrc$consequence)

plot(density(tmrc$w[which(tmrc$consequence=="synonymous")]))
lines(density(tmrc$w[which(tmrc$consequence=="nonsynonymous")]), col = "red")
t.test(tmrc$w[which(tmrc$consequence=="nonsynonymous")] , tmrc$w[which(tmrc$consequence=="synonymous")])

for(mut in allMutations){
  tt <- tmrc[which(tmrc$mut==mut) , ]
  MIN_OBS <- 10
  if( length(which(tt$consequence=="nonsynonymous"))>MIN_OBS && length(which(tt$consequence=="synonymous"))>MIN_OBS ){
    cat(mut, ":\n")
    tt <- t.test(tt$w[which(tt$consequence=="nonsynonymous")] , tt$w[which(tt$consequence=="synonymous")])
    print(tt)
    cat("\n---------\n")
  }
}


#MAX_W <- 1.5
tmrc$w[which(tmrc$w<0)] <- 0
tmrc$w[which(tmrc$w>MAX_W)] <- MAX_W
tmrc$SYN <- F
tmrc$SYN[which(tmrc$consequence == "synonymous")] <- T
tmrc$SYN <- factor(tmrc$SYN)

myColors <- c("red", "steelblue")

tmp <- tmrc
# https://r-charts.com/distribution/histogram-group-ggplot2/
gp <- ggplot(data = tmp, aes(x = w, fill = SYN, colour = SYN)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15)
plot(gp)
t.test(tmrc$w[which(tmrc$consequence=="synonymous")] , tmrc$w[which(tmrc$consequence!="synonymous")])

tmp <- tmrc[which(tmrc$consequence != "nonsense") , ]
gp <- ggplot(data = tmp, aes(x = w, fill = consequence, colour = consequence)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15) + 
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.background = element_rect(colour="white", fill="white"),
    legend.position=c(0.8,0.8),
    legend.box.background=element_rect(),legend.box.margin=margin(5,5,5,5)
    )
plot(gp)
write.csv(tmp[ , c("uid", "consequence", "w")], file = "fitness-vs-consequence-no-nonsense.csv", row.names = F)


tmpct <- tmrc[which(tmrc$mut == "C->T" & tmrc$consequence != "nonsense"), ]
gp <- ggplot(data = tmpct, aes(x = w, fill = consequence, colour = consequence)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15)
plot(gp)
t.test(tmpct$w[which(tmpct$consequence=="synonymous")] , tmpct$w[which(tmpct$consequence!="synonymous")])


tmp <- tmrc[which(tmrc$mut == "C->T" & tmrc$consequence == "nonsynonymous"), ]
tmp$type <- "structural"
tmp$type[which(tmp$position<21560)] <- "nonstructural"
tmp$type <- factor(tmp$type)
gp <- ggplot(data = tmp, aes(x = w, fill = type, colour = type)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15)
plot(gp)



tmpct <- tmrc[which(tmrc$consequence != "nonsense"), ]
tmpct$mutation_type <- "Other"
tmpct$mutation_type[which(tmpct$mut=="C->T")] <- "C-to-U"
tmpct$mutation_type <- factor(tmpct$mutation_type)
gp <- ggplot(data = tmpct, aes(x = w, fill = mutation_type, colour = mutation_type)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15)
plot(gp)


t1 <- tmrc[which(tmrc$consequence=="synonymous" & tmrc$mut=="C->T") , ]
t2 <- tmrc[which(tmrc$consequence=="synonymous" & tmrc$mut!="C->T") , ]
v1 <- c( length(which(t1$w == 0)) , length(which(t2$w == 0)) )
v2 <- c( nrow(t1), nrow(t2) )
prop.test(v1, v2)
length(which(t1$w==0))/nrow(t1)
length(which(t2$w==0))/nrow(t2)


tmp <- tmrc
t.test(tmp$w[which(tmp$consequence=="synonymous" & tmp$mut=="C->T")] , 
       tmp$w[which(tmp$consequence=="synonymous" & tmp$mut!="C->T")])

v1 <- c( length(which(tmp$w==0 & tmp$mut=="C->T")) , length(which(tmp$w==0 & tmp$mut!="C->T")) )
v2 <- c( length(which(tmp$mut=="C->T")) , length(which(tmp$mut!="C->T")) )
prop.test(v1, v2)

tmpct <- tmrc[tmrc$mut == "C->T" , ]
tmp <- tmrc[which(tmrc$mut != "C->T") , ]

vs1 <- c( length(which(tmp$w==0 & tmp$SYN==T)) , length(which(tmpct$w==0 & tmpct$SYN==T)) )
vs2 <- c( length(which(tmp$SYN==T)) , length(which(tmpct$SYN==T)))
prop.test(vs1, vs2)

vn1 <- c( length(which(tmp$w==0 & tmp$SYN==F)) , length(which(tmpct$w==0 & tmpct$SYN==F)) )
vn2 <- c( length(which(tmp$SYN==F)) , length(which(tmpct$SYN==F)))
prop.test(vn1, vn2)



MAX_LETHAL <- 0.1
tsyn <- tmrc[which(tmrc$consequence=="synonymous") , ]
tnon <- tmrc[which(tmrc$consequence!="synonymous") , ]
vLethal <- c( length(which(tsyn$w < MAX_LETHAL)) , length(which(tnon$w < MAX_LETHAL)) )
vDist <- c( nrow(tsyn), nrow(tnon) )
prop.test(vLethal, vDist)

NEUTRAL_LOWER_RANGE <- 0.6
NEUTRAL_UPPER_RANGE <- 1.4
vNeutral <- c( length(which( tsyn$w > NEUTRAL_LOWER_RANGE & tsyn$w < NEUTRAL_UPPER_RANGE)) , 
               length(which( tnon$w > NEUTRAL_LOWER_RANGE & tnon$w < NEUTRAL_UPPER_RANGE)) )

prop.test( vNeutral, vDist)


tExp <- tmrc[ , c("uid", "position", "ref", "bTo", "mut", "w", "geneID", "Name", "refAA", "proteinLoc", "varAA", "consequence")]
write.csv(tExp, file = "fitness-mapped-to-protein.csv", row.names = F)
#write.csv(tExp, file = "fitness-mapped-to-protein-ORF1ab-w20.csv", row.names = F)


################################################################################
# Highest fitness per codon:
tExp <- read.csv("fitness-mapped-to-protein.csv")
tExp$puid <- paste(tExp$geneID, tExp$position, tExp$proteinLoc, sep = "_")
tns <- tExp[which(tExp$consequence == "nonsynonymous") , ]

vns <- by(tns$w, tns$puid, max)
tns <- data.frame(
  puid = names(vns),
  w = as.numeric(vns)
)
vs <- strsplit(tns$puid, split = "_")
vsu <- unlist(vs)
tns$geneID <- vsu[seq(from=1, to=length(vsu)-2, by = 3)]
tns$gPos <- as.numeric(vsu[seq(from=2, to=length(vsu)-1, by = 3)])
tns$protLoc <- as.numeric(vsu[seq(from=3, to=length(vsu), by = 3)])

tnsr <- tns[order(tns$geneID, tns$protLoc, decreasing = F) , c("geneID", "gPos", "protLoc", "w")]
write.csv(tnsr, file="fitness-per-protein-per-codon.csv", row.names = F)

################################################################################
# Adding the amino acid for the NSPs:
genome <- readDNAStringSet("../DB/ensembl_ASM985889v3/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
tnsp <- read.csv("../DB/ensembl_ASM985889v3/nsp.csv")
tf <- tnsr
tf$NSP <- NA
tf$refAA <- NA
for(i in (1:nrow(tnsp))){
  print(i)
  nsp <- tnsp$nsp[i]
  start <- tnsp$start[i]
  end <- tnsp$end[i]
  geneSeq <- genome[[1]][start:end]
  protSeq <- translate(geneSeq)
  w <- which(tf$gPos>=start & tf$gPos<=end)
  if( length(w) > 0 ){
    tf$NSP[w] <- nsp
    tf$protLoc[w] <- ceiling( ((tf$gPos[w]-start)+1)/3 )
    tf$refAA[w] <- unlist(strsplit(as.character(protSeq[tf$protLoc[w]]) , split = ""))
    
    tfr <- tf[which(tf$NSP == nsp) , ]
    fName <- paste(nsp, "-fitness.csv", sep = "")
    write.csv(tfr, file = fName, row.names = F)
  }
}

################################################################################
# Get value per codon for 3D figure:

tt <- tmrc[which(tmrc$consequence=="nonsynonymous") , c("uid", "mut", "geneID", "proteinLoc", "consequence", "w") , ]

geneID <- "ENSSASG00005000004"

ttg <- tt[which(tt$geneID == geneID) , ]
vMax <- by(ttg$w, ttg$proteinLoc, max)
tmax <- data.frame(
  proteinPos = as.numeric(names(vMax)),
  w = as.numeric(vMax)
)
tap <- unique(tmrc[which(tmrc$geneID == geneID) , c("proteinLoc", "aa")])
tm <- merge(tmax, tap, by.x = "proteinPos", by.y = "proteinLoc")
fName <- paste(geneID, "-w-per-position.csv", sep = "")
write.csv(tmax, file = fName)

maxPos <- 1147
vv <- rep(-1.0, maxPos)
#vv <- round(rnorm(n = maxPos, mean = 0.5, sd = 0.2), 2)
w <- which(tmax$proteinPos<=maxPos)
vv[tmax$proteinPos[w]] <- round(tmax$w[w], 2)/2
write(vv, file = "beta-values.txt", sep = "\n")


################################################################################
# All fitness values (possibly multiple per codon) for the spike protein:

tExp <- read.csv("fitness-mapped-to-protein.csv")
tExp$puid <- paste(tExp$geneID, tExp$position, tExp$proteinLoc, sep = "_")
tns <- tExp[which(tExp$consequence == "nonsynonymous") , ]

tr <- mcols(res)[ , c("GENEID", "uid", "REF", "varAllele", "REFAA", "VARAA", "geneName")]
trs <- tr[which(tr$geneName == "S") , ]

tm <- merge(tns, trs, by = "uid")
dim(tm)
length(unique(tm$position))

tmr <- tm[ , c("uid", "position", "w", "proteinLoc", "REFAA", "VARAA")]
write.csv(tmr, file = "fitness-spike-all-mutations.csv", row.names = F)
