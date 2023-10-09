BD_DATA <- ""
DB_DIR <- ""

library(ggplot2)

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

fname <- paste(BD_DATA, "all-frequencies-combined-and-lethal.tab.gz", sep="/")
tmuts <- read.table(fname, h=T)

lethal_uids <- tmuts$uid[which(tmuts$lethal==T)]

tExp <- read.csv("fitness-mapped-to-protein.csv")
tw <- tExp
tw$lethal <- F
tw$lethal[which(is.element(tw$uid, lethal_uids)==T)] <- T


tvcf <- read.table("../DB/ensembl_ASM985889v3/variants/vcf_usher_ensembl.tab", h = T)
tvcf <- read.table("../UShER/2023-05-23/vcf_usher_ensembl.tab", h = T)
tvcf$total <- log(tvcf$count + tvcf$ensembl + 0.1)

tvcf$lethal <- F
tvcf$lethal[which(is.element(tvcf$uid, lethal_uids)==T)] <- T

twr <- tw[ , c("uid", "w")]
tm <- merge(tvcf, twr, by = "uid", all.x = T)

res <- readRDS("../DB/protein_predict.rds")
tm$consequence <- "nonsynonymous"
tm$consequence[which(is.element(tm$uid, res$uid[which(res$CONSEQUENCE=="synonymous")]))] <- "synonymous"
tm$consequence[which(is.element(tm$uid, res$uid[which(res$CONSEQUENCE=="nonsense")]))] <- "nonsense"


tmr <- tm[which(is.na(tm$w)==F) , ]
dim(tmr)
tmr$newLethal <- F
tmr$newLethal[which(tmr$w==0 & tmr$lethal==F)] <- T

length(which(tmr$w==0 & tmr$lethal==F))/length(which(tmr$w==0))

tmn <- tmr[which(tmr$lethal==F) , ]
# Comparing the number of UShER-aligned genomes in which a mutation appears 
# for those classified as lethal vs the rest of mutations.
t.test(tmn$total[which(tmn$w==0)], tmn$total[which(tmn$w>0)])
t.test(tmn$count[which(tmn$w==0)], tmn$count[which(tmn$w>0)])

gp <- ggplot(tmn, aes(x = total, fill = newLethal, colour = newLethal)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 20)
plot(gp)

# Same thing, but specifically for synonymous mutations.
tmn <- tmr[which(tmr$lethal==F & tmr$consequence=="synonymous") , ]
t.test(tmn$total[which(tmn$w==0)], tmn$total[which(tmn$w>0)])
t.test(tmn$count[which(tmn$w==0)], tmn$count[which(tmn$w>0)])

# Looking only at synonymous C-to-U mutations
tmn <- tmr[which(tmr$lethal==F & tmr$consequence=="synonymous" & tmr$from=="C" & tmr$to=="T") , ]
t.test(tmn$total[which(tmn$w==0)], tmn$total[which(tmn$w>0)])
t.test(tmn$count[which(tmn$w==0)], tmn$count[which(tmn$w>0)])


tmn <- tmr[which(tmr$lethal==F & tmr$consequence=="nonsynonymous") , ]
t.test(tmn$total[which(tmn$w==0)], tmn$total[which(tmn$w>0)])
t.test(tmn$count[which(tmn$w==0)], tmn$count[which(tmn$w>0)])



################################################################################
# Split by consequence

library(tidyr)

res <- readRDS("../DB/protein_predict.rds")
rr <- mcols(res)[ , c("uid", "CONSEQUENCE")]

tt <- merge(tmr, rr, by = "uid")
ttr <- tt[which(tt$CONSEQUENCE=="synonymous" & tt$lethal==F)  , ]
#ttr <- tt[which(tt$CONSEQUENCE=="nonsynonymous" & tt$lethal==F)  , ]

t.test(ttr$total[which(ttr$w==0)] , ttr$total[which(ttr$w>0)])
# Even for synonymous mutations, those with w = 0 are significantly more rare among the genome alignments.

t.test(ttr$count[which(ttr$w==0)] , ttr$count[which(ttr$w>0)])
