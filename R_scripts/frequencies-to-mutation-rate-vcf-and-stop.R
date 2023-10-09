################################################################################
#
# Using pre-calculated list of lethal mutations to obtain mutation rates
# (= mutation frequency for lethal mutations)
#

#options(width=Sys.getenv("COLUMNS"))

setwd("/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/RESULTS_BY_VARIANT")
source("/mnt/scratch/goutlab/RollingCircle/R_files/seqan-functions-v3.R")
source("/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/R_files/covid-functions.R")


BD_DATA = "/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/RESULTS_BY_VARIANT/DELTA/B"
DB_DIR = "/mnt/scratch/goutlab/RollingCircle/SARS_CoV_2/DB/ensembl_ASM985889v3/"

library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(genomation)
library(readxl)

vBases = getBases()
vMuts = getMutsFromBases(vBases)

NB_MAX_BEFORE_HS = 1e12
MIN_COV = 10
MAX_DIV = 10

fName = paste("lRes-MaxDiv", MAX_DIV, "-MinCov-", MIN_COV, "-HS", NB_MAX_BEFORE_HS, ".bin", sep="")
load(fName)


################################################################################
# These are variants never observed in the alignment of multiple genomes
vcf <- read.table(paste(DB_DIR, "/variants/vcf_usher_ensembl.tab", sep="/"), h=T)
load("mutations-track-all.bin") # <- tfull
lethal <- vcf$uid[which(vcf$count==0 & vcf$ensembl==0)]

# These are the variants that create stop codons
load(paste(DB_DIR, "res_lethal_stop.bin", sep="/"))
stop_really_lethal <- res_lethal_stop[which(res_lethal_stop$GENEID=="ENSSASG00005000002" & as.numeric(res_lethal_stop$PROTEINLOC)<4000)]
v_lethal_stops <- paste(start(stop_really_lethal), stop_really_lethal$REF, stop_really_lethal$varAllele, sep="_")
length(v_lethal_stops)
length(which(is.element(v_lethal_stops, lethal)==T))

lethal <- unique(lethal, v_lethal_stops)         

tfull_lethal <- tfull[which(is.element(tfull$uid, lethal)==T) , ]
t_non_lethal <- tfull[which(is.element(tfull$uid, lethal)==F) , ]

tlr <- unique(tfull_lethal[ , c("uid", "consequence")])
tnr <- unique(t_non_lethal[ , c("uid", "consequence")])

v1 <- c(length(which(tlr$consequence=="synonymous")), length(which(tnr$consequence=="synonymous")) )
v2 <- c( nrow(tlr), nrow(tnr) )
prop.test(v1, v2)


tfr <- tfull_lethal

tfr <- tfull_lethal[which(tfull_lethal$cov>0 & tfull_lethal$rate<0.01 & tfull_lethal$nb < 10) , ]
tfn <- tfull_lethal[which(tfull_lethal$cov>0 & tfull_lethal$rate<0.01 & tfull_lethal$nb < 10) , ]



tfr$sampleName <- paste(tfr$variant, tfr$replicate, tfr$passage, sep="_")
tfr$sampleGroup <- paste(tfr$variant, tfr$replicate, sep="_")
tfr$mut <- substr(tfr$uid, nchar(tfr$uid)-2, nchar(tfr$uid))

vMuts <- unique(tfr$mut)
vSampleGroups <- unique(tfr$sampleGroup)
vSampleNames <- unique(tfr$sampleName)

tm = data.frame(mutation = character(),
                   nbMut = numeric(),
                   nbObs = numeric(),
                   rate = numeric(),
                   variant = character(),
                   replicate = character(),
                   passage = character(),
                   sampleName = character(),
                   sampleGroup = character()
)

tmg <- data.frame(
  mutation = character(),
  nbMut = numeric(),
  nbObs = numeric(),
  sampleGroup = character()
)

for(sampleGroup in vSampleGroups){
  ts <- tfr[which(tfr$sampleGroup == sampleGroup) , ]
  for(mutation in vMuts){
    tsm <- ts[which(ts$mut == mutation) , ]
    nbMut <- sum(tsm$nb)
    nbObs <- sum(tsm$cov)
    tmp <- data.frame(
      mutation = mutation,
      nbMut = nbMut,
      nbObs = nbObs,
      sampleGroup = sampleGroup
    )
    tmg <- rbind(tmg, tmp)
  }
}

tmg$rate <- tmg$nbMut / tmg$nbObs

tr <- tmg[which(is.element(tmg$sampleGroup, c("USA_A", "USA_B", "ALPHA_A", "ALPHA_B", "DELTA_A", "DELTA_B", "OMICRON_A"))==T) , ]
tr <- tmg[which(is.element(tmg$sampleGroup, c("DELTA-ALI_B", "DELTA-CALU_B"))==F) , ]
tr <- tmg[which(is.element(tmg$sampleGroup, c("USA_A", "USA_B", "ALPHA_A", "ALPHA_B", "DELTA_A", "DELTA_B"))==T) , ]


tr <- tr[which(tr$mutation!="C_T"),]

tr <- tmg

pdf("mutation-rates-all-samples-all-mutations.pdf", width=16)
gp <- ggplot(data = tr, aes(x = mutation, y = rate, fill = sampleGroup) ) + 
  geom_bar(stat="identity", position=position_dodge())
plot(gp)
dev.off(dev.cur())







for(sampleName in vSampleNames){
  ts <- tfr[which(tfr$sampleName == sampleName) , ]
  variant <- unlist(strsplit(sampleName, "_"))[1]
  replicate <- unlist(strsplit(sampleName, "_"))[2]
  passage <- unlist(strsplit(sampleName, "_"))[3]
  sampleGroup <- paste(variant, replicate, sep="_")
  for(mutation in vMuts){
    tsm <- ts[which(ts$mut == mutation) , ]
    nbMut <- sum(tsm$nb)
    nbObs <- sum(tsm$cov)
    tmp <- data.frame(mutation = mutation,
                      nbMut = nbMut,
                      nbObs = nbObs,
                      rate = nbMut/nbObs,
                      variant = variant,
                      replicate = replicate,
                      passage = passage,
                      sampleName = sampleName,
                      sampleGroup = sampleGroup
                      )
    tm <- rbind(tm, tmp)
  }
}




tmr = tm[which(tm$nbMut>10 & is.element(tm$variant, c("ALPHA", "USA", "DELTA", "OMICRON", "GAMMA"))==T) , ]
tmr = tm[which(tm$nbMut>10 & is.element(tm$sampleGroup, c("BETA_A", "DELTA-ALI_B", "DELTA-CALU_B"))==F) , ]

pdf(file="mutations-rates.pdf", width=20)
for(mutation in unique(tm$mutation)){
  tmp = tmr[which(tmr$mutation==mutation) , ]
  boxplot(tmp$rate~tmp$variant, main=mutation)
}
dev.off(dev.cur())  



tda <- tfull[which(tfull$mutation!="C_T") , ]
tda <- tda[which( is.element(tda$sampleGroup, c("ALPHA_A", "ALPHA_B", "USA_A", "USA_B", "DELTA_A", "DELTA_B"))==T) ,  ]
pdf(file="test.pdf", width=24)
gp <- ggplot(data= tda, aes(x = mutation, y = rate, fill = sampleGroup)) +
    geom_boxplot()
plot(gp)
dev.off(dev.cur())





