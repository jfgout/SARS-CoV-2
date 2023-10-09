################################################################################
#
# Searching for mutations that show large changes in frequency between passages.
# This might reveal adaptation to Vero E6 cells or other adaptations.
#


setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")

library(tidyr)
library(rtracklayer)

# Loading the file with the consequences of every possible mutation (inside coding sequences)
res <- readRDS("../DB/protein_predict.rds")
# Adding gene names:
#res$geneName <- NA
##gff <- import.gff("../DB/ensembl_ASM985889v3/Sars_cov_2.ASM985889v3.101.gff3")
#genes <- gff[which(gff$type == "gene")]
#for(i in (1:length(genes))){
#  geneID <- genes$gene_id[i]
#  geneName <- genes$Name[i]
#  res$geneName[which(res$GENEID==geneID)] <- geneName
#}
res$mutation <- paste(res$Name, ":", res$REFAA, res$PROTEINLOC, res$VARAA, sep = "")

tf <- readRDS("all-obs-and-frequencies-combined-and-lethal.rds")
vSeries <- c("ALPHA_A", "ALPHA_B", "USA_A", "USA_B", "DELTA_A", "DELTA_B")

MIN_DELTA_FREQ <- 0.01
MIN_COVERAGE <- 1000 # Minimum number of reads covering the position


tf$min <- NA
tf$max <- NA
tf$delta <- NA
tf$keep <- F # This column will be set to TRUE for mutations with large frequency changes in at least one series.
tf$keepSeries <- "" # This column will keep track of which series were responsible for a large frequency change.

vDeltas <- c()

for(series in vSeries){
  for(passage in seq(from = 1, to = 7, by = 1)){
    cname_freq <- paste(series, passage, "rate", sep = "_")
    cname_cov <- paste(series, passage, "cov", sep = "_")
    
    wLowCoverage <- which(tf[ , cname_cov] < MIN_COVERAGE)
    
    if( length(wLowCoverage) > 0 ){
      tf[wLowCoverage , cname_freq] <- NA
    }
    
  }
  cnames <- paste(series, seq(from=1, to = 7, by=1), "rate", sep = "_")
  tf$min <- apply(tf[ , cnames], 1, min, na.rm = T)
  tf$max <- apply(tf[ , cnames], 1, max, na.rm = T)
  tf$delta <- tf$max - tf$min
  wk <- which(tf$delta >= MIN_DELTA_FREQ & tf$min!=Inf & tf$max != Inf) # The few positions that have only NAs would return Inf
  tf$keep[wk] <- TRUE
  tf$keepSeries[wk] <- paste(tf$keepSeries[wk], series, sep = ",")
  
  tfok <- tf[which(tf$min != Inf & tf$max != Inf) , ]
  vtmp <- tfok$delta
  names(vtmp) <- paste(series, tfok$uid, sep = "-")
  vDeltas <- c(vDeltas, vtmp)
  #vDeltas <- c(vDeltas, tf$delta[which(tf$min!=Inf & tf$max!=Inf)])
}

tfr <- tf[which(tf$keep == TRUE) , ]
dim(tfr)
hist(tfr$position, nclass = 20)

vp <- vDeltas[which(vDeltas >= MIN_DELTA_FREQ )]
vpu <- unlist(strsplit(names(vp), "-"))
td <- data.frame(
  uid = vpu[seq(from=2, to = length(vpu), by = 2)],
  series = vpu[seq(from=1, to = length(vpu)-1, by = 2)],
  delta = as.numeric(vp)
)

vu <- unlist(strsplit(td$series, "_"))
td$variant <- vu[seq(from=1, to=length(vu)-1, by = 2)]
td$replicate <- vu[seq(from=2, to=length(vu), by = 2)]

vu <- unlist(strsplit(td$uid, "_"))
td$gPos <- vu[seq(from=1, to = length(vu)-2, by = 3)]
td$bFrom <- vu[seq(from=2, to = length(vu)-1, by = 3)]
td$bTo <- vu[seq(from=3, to = length(vu), by = 3)]

rr <- as.data.frame(mcols(res)[ , c("uid", "mutation", "Name", "REFAA", "VARAA")])
tdm <- merge(td, rr, by = "uid", all.x = T)

tdm$consequence <- "nonsynonymous"
tdm$consequence[which(tdm$REFAA != "*" & tdm$VARAA=="*")] <- "nonsense"
tdm$consequence[which(tdm$REFAA == tdm$VARAA)] <- "synonymous"

write.csv(tdm, file = "large-frequency-mutations.csv", row.names = F)

################################################################################
# Is there an excess of non-synonymous mutations among the ones with large 
# frequency change ?

tdu <- unique(tdm[which(is.na(tdm$Name)==F) , c("uid", "consequence")])

v1 <- c( length(which(tdu$consequence=="nonsynonymous")) , length(which(res$CONSEQUENCE=="nonsynonymous")) )
v2 <- c(  nrow(tdu) , length(res) )
prop.test(v1, v2)

v1 <- c( length(which(tdu$consequence=="synonymous")) , length(which(res$CONSEQUENCE=="synonymous")) )
v2 <- c(  length(which(tdu$consequence=="nonsynonymous")) , length(which(res$CONSEQUENCE=="nonsynonymous")) )
prop.test(v1, v2)
# No Excess of synonymous substitutions!

v1 <- c( length(which(tdu$consequence=="nonsense")) , length(which(res$CONSEQUENCE=="nonsense")) )
v2 <- c(   length(vu) , length(res) )
prop.test(v1, v2)
# Significant depletion of nonsense mutations


################################################################################
# Is there an excess of mutations shared with the defining mutations from 
# covariants (https://covariants.org/variants/21L.Omicron)

tdf <- read.csv("../DB/defining-mutations.csv")
tdfr <- tdf[which(substr(tdf$mutation, nchar(tdf$mutation), nchar(tdf$mutation)) != "-") , ]
nrow(tdf)
nrow(tdfr)

nbCom <- length(which(is.element(tdm$mutation, tdfr$mutation)==T))
nbCom
NB <- nrow(tdm)
res_n <- res[which(res$CONSEQUENCE != "synonymous")]
NB_RANDS <- 1000
randRes <- rep(NA, NB_RANDS)
for(i in (1:NB_RANDS)){
  vMuts <- sample(res_n$mutation, NB)
  randRes[i] <- length(which(is.element(vMuts, tdfr$mutation)==T))
}
hist(randRes)
nbCom
# Yes, large excess. But this resource provides information for non-synonymous mutations only.


################################################################################
# Is there an excess overlap with the large fitness mutations from pyrocov?

tdm <- read.csv("large-frequency-mutations.csv")

tdu <- unique(tdm[which(is.na(tdm$Name)==F) , c("uid", "mutation", "consequence")])
tns <- tdu[which(tdu$consequence=="nonsynonymous") , ]

tpy <- read.table("../PyroCov/mutations.tsv", h = T, sep = "\t")
NB_OVERLAP <- length(which(is.element(tns$mutation, tpy$mutation)==T))
NB_OVERLAP

res <- readRDS("../DB/protein_predict-ORF1ab.rds")
res$mutation <- paste(res$Name, ":", res$REFAA, res$PROTEINLOC, res$VARAA, sep = "")
res_n <- res[which(res$CONSEQUENCE == "nonsynonymous")]

NB <- nrow(tns)
NB_RAND <- 1000
resRand <- rep(NA, NB_RAND)
for(i in (1:NB_RAND)){
  vMuts <- sample(res_n$mutation, NB)
  resRand[i] <- length(which(is.element(vMuts, tpy$mutation)==T))
}
hist(resRand)
median(resRand)
NB_OVERLAP

# Yes, large excess. But this resource provides information for non-synonymous mutations only.



#######################################
# Comparison with the JBloom data:
tlf <- read.csv("large-frequency-mutations.csv")
res <- readRDS("../DB/protein_predict.rds")
res$Name[which(start(res)<21555)] <- "ORF1ab"
w <- which(res$Name=="ORF1ab")

res$mutation <- paste(res$Name, ":", res$REFAA, res$PROTEINLOC, res$VARAA, sep = "")
dr <- as.data.frame(mcols(res)[ , c("uid", "mutation")])

tlfr <- tlf[ , which(colnames(tlf)!="mutation")]
tlf <- merge(tlfr, dr, by = "uid")
dim(tlf)
tdu <- unique(tlf[which(is.na(tlf$Name)==F) , c("uid", "mutation", "consequence")])

tdu_s <- tdu[which(tdu$consequence=="synonymous") , ]
tmiss <- tdu_s[which(is.element(tdu_s$mutation, tj$mutation)==F) , ]
tmiss

tj <- read.csv("../jbloomlab/aamut_fitness_all.csv")
tj$SYN <- F
tj$SYN[which(tj$clade_founder_aa==tj$mutant_aa)] <- T
tj$mutation <- paste(tj$gene, ":", tj$aa_mutation, sep = "")

tj$LF <- F
tj$LF[which(is.element(tj$mutation, tdu$mutation)==T)] <- T
length(which(tj$LF==T))

tjs <- tj[which(tj$SYN==T) , ]
length(which(tjs$LF==T))

gp <- ggplot(data = tjs, aes(x = delta_fitness, fill = LF, colour = LF)) +
  geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'dodge', bins = 15)
plot(gp)
t.test( tjs$delta_fitness[which(tjs$LF==F)] , tjs$delta_fitness[which(tjs$LF==T)] )

# --> Looking at genomes alignments.

vcf <- read.table("../DB/ensembl_ASM985889v3/variants/vcf_usher_ensembl.tab", h = T) # This file contains a summary of the frequency of each mutation in the multiple (~6 million) genomes alignment
vcf <- read.table("../UShER/2023-05-23/vcf_usher_ensembl.tab", h = T)
vcf$total <- vcf$count + vcf$ensembl
vcf$total <- log2(vcf$count + vcf$ensembl + 1)
vcf$LF <- F # This column will indicate if a mutation is among those with large frequency change or not
w <- which(is.element(vcf$uid, tdm$uid)==T)
vcf$LF[w] <- T

t.test(vcf$total[which(vcf$LF==T)], vcf$total[which(vcf$LF==F)])

# Looking specifically at synonymous mutations:
vcf$SYN <- F
syn_uids <- res$uid[which(res$CONSEQUENCE=="synonymous")]
w <- which(is.element(vcf$uid, syn_uids)==T)
vcf$SYN[w] <- T

vcf_syn <- vcf[which(vcf$SYN==T) , ]
dim(vcf_syn)
t.test(vcf_syn$total[which(vcf_syn$LF==T)], vcf_syn$total[which(vcf_syn$LF==F)])
t.test(vcf_syn$ensembl[which(vcf_syn$LF==T)], vcf_syn$ensembl[which(vcf_syn$LF==F)])
t.test(vcf_syn$count[which(vcf_syn$LF==T)], vcf_syn$count[which(vcf_syn$LF==F)])



################################################################################
# What is the fitness distribution of mutations with large frequency changes?

tlf <- read.csv("large-frequency-mutations.csv")

tf <- read.csv("../ALL_DATA/fitness-mapped-to-protein.csv")

tf$largeChange <- F
tf$largeChange[which(is.element(tf$uid, tlf$uid)==T)] <- T

t.test(tf$w[which(tf$largeChange==F)] , tf$w[which(tf$largeChange==T)])