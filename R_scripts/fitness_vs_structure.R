# Testing the hypothesis that mutations are more deletrious if they arise in 
# a region of the genome that is part of a double-stranded structure.

setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
library(ggplot2)
library(tidyverse)

tf <- read.csv("fitness-mapped-to-protein.csv")
tf$w[which(tf$w<0)] <- 0
tf$w[which(tf$w>1.5)] <- 1.5

ts <- read.table("../Structure/2022_NatCom_Lan/structure.ct", sep = "\t", h = F)
colnames(ts) <- c("position", "base", "prev_index", "next_index", "nbp", "nm")
pos_single <- ts$position[which(ts$nbp==0)]
length(pos_single)

tf$structure <- "Paired"
tf$structure[which(is.element(tf$position, pos_single)==T)] <- "Unpaired"

gp <- ggplot(data = tf[which(tf$consequence=="synonymous" & tf$mut=="C->T") , ], aes(x = w, fill = structure, colour = structure)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15) +
  ggtitle("Synonymous C-to-U mutations") + 
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.background = element_rect(colour="white", fill="white"),
    legend.position=c(0.8,0.8),
    legend.box.background=element_rect(),legend.box.margin=margin(5,5,5,5),
    plot.title = element_text(hjust = 0.5)
  )
plot(gp)

tfr <- tf[ , c("uid", "structure", "w")]
write.csv(tfr, file = "fitness-vs-structure-all.csv", row.names = F)

tfr <- tf[ which(tf$consequence == "synonymous" & tf$mut == "C->T") , ]
write.csv(tfr[ , c("uid", "structure", "w")], file = "fitness-vs-structure-synonymous_C-to-U-only.csv", row.names = F)

t.test(tfr$w[which(tfr$structure == "Paired")] , tfr$w[which(tfr$structure=="Unpaired")])


tfr <- tf[ which(tf$consequence == "nonsynonymous" & tf$mut == "C->T") , ]
t.test(tfr$w[which(tfr$structure == "Paired")] , tfr$w[which(tfr$structure=="Unpaired")])


gp <- ggplot(data = tf[which(tf$consequence=="nonsynonymous" & tf$mut=="C->T") , ], aes(x = w, fill = structure, colour = structure)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15)
plot(gp)


gp <- ggplot(data = tf[which(tf$consequence=="synonymous" & tf$mut!="C->T") , ], aes(x = w, fill = structure, colour = structure)) +
  geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15) +
  ggtitle("Synonymous other mutations")
plot(gp)


for(mut in unique(tf$mut)){
  gp <- ggplot(data = tf[which(tf$consequence=="synonymous" & tf$mut==mut) , ], aes(x = w, fill = single, colour = single)) +
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'dodge', bins = 15) +
    ggtitle(mut)
  plot(gp)
  readline(prompt="Press [enter] to continue")
}


tf <- read.csv("fitness-mapped-to-protein-ORF1ab.csv")
tf$structure <- "Paired"
tf$structure[which(is.element(tf$position, pos_single)==T)] <- "Unpaired"
tf$mutation <- paste(tf$Name, ":", tf$refAA, tf$proteinLoc, tf$varAA, sep = "")

CLADE <- "21J"
tp <- read.csv("../jbloomlab/aamut_fitness_by_clade.csv")
tp$mutation <- paste(tp$gene, ":", tp$clade_founder_aa, tp$aa_site, tp$mutant_aa, sep = "")
tpr <- tp[which(tp$clade==CLADE) , c("mutation", "delta_fitness", "expected_count", "count_terminal")]
tj <- tpr

#tj <- read.csv("../jbloomlab/aamut_fitness_all.csv")
#tj$mutation <- paste(tj$gene, ":", tj$aa_mutation, sep = "")
tm <- merge(tf, tj, by = "mutation")
dim(tm)

gp <- ggplot(data = tm[which(tm$consequence=="synonymous" & tm$mut!="C->T") , ], aes(x = delta_fitness, fill = structure, colour = structure)) +
  geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'dodge', bins = 15)
plot(gp)

gp <- ggplot(data = tm[which(tm$consequence=="nonsynonymous" & tm$mut=="C->T") , ], aes(x = delta_fitness, fill = structure, colour = structure)) +
  geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'dodge', bins = 15)
plot(gp)

gp <- ggplot(data = tm[which(tm$consequence=="synonymous" & tm$mut!="C->T") , ], aes(x = delta_fitness, fill = structure, colour = structure)) +
  geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'dodge', bins = 15)
plot(gp)


tmsct <- tm[which(tm$mut == "C->T" & tm$consequence == "synonymous") , ]
dim(tmsct)
boxplot(tmsct$delta_fitness~tmsct$structure)
t.test(tmsct$delta_fitness[which(tmsct$structure=="Paired")] , tmsct$delta_fitness[which(tmsct$structure=="Unpaired")])


res <- readRDS("../DB/protein_predict-ORF1ab.rds")
res$gPos <- start(res)
res$mutation <- paste(res$Name, ":", res$REFAA, res$PROTEINLOC, res$VARAA, sep = "")
res$mut <- paste(res$REF, "->", res$varAllele, sep = "")
rd <-  as.data.frame(mcols(res)[ , c("gPos", "mutation", "CONSEQUENCE", "mut")])

tm <- merge(tpr, rd, by = "mutation")
dim(tm)
tm$structure <- "Paired"
tm$structure[which(is.element(tm$gPos, pos_single)==T)] <- "Unpaired"

tmp <- tm[which(tm$CONSEQUENCE=="synonymous" & tm$mut == "C->T") , ]
tmp <- tm[which(tm$CONSEQUENCE=="synonymous" & tm$mut != "C->T") , ]
for(mut in unique(tm$mut)){
  tmp <- tm[which(tm$CONSEQUENCE=="synonymous" & tm$mut == mut) , ]
  gp <- ggplot(data = tmp, aes(x = delta_fitness, fill = structure, colour = structure)) + 
    geom_histogram(alpha = 0.3, aes(y = after_stat(density)), position = 'dodge', bins = 15) +
    ggtitle(mut)
  plot(gp)
  tt <- t.test(tmp$delta_fitness[which(tmp$structure=="Paired")] , tmp$delta_fitness[which(tmp$structure=="Unpaired")])
  print(tt)
  readline(prompt="Press [enter] to continue")
}

tms <- tm[which(tm$CONSEQUENCE=="synonymous") , ]
tms$rp <- round(tms$gPos/3000, 0)
boxplot(tms$delta_fitness~tms$rp)


################################################################################
# Using SHAPE scores from Sun et al. 2021

ts <- read.csv("../Structure/2021_Sun_SHAPE/shape_scores.csv")
ts$shapeScore <- as.numeric(as.character(ts$shapeScore))
tms <- merge(tm, ts, by.x = "gPos", by.y = "position")
dim(tms)
tmsru <- unique(tms[ , c("gPos", "structure", "shapeScore")])
boxplot(tmsru$shapeScore~tmsru$structure)

tmp <- tms[which(tms$mut == "C->T" & tms$CONSEQUENCE == "synonymous") , ]
cor(tmp$delta_fitness, tmp$shapeScore)
plot(tmp$delta_fitness, tmp$shapeScore)

# Comparison between SHAPE reactivity score and fitness value from our study:
tsf <- merge(tf, ts, by = "position")
dim(tsf)
tmp <- tsf[which(tsf$mut == "C->T" & tsf$consequence == "synonymous") , ]
cor(tmp$w,tmp$shapeScore)
summary(lm(tmp$w~tmp$shapeScore))
