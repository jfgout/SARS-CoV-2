# Looking at the distribution of fitness value on a gene by gene basis.
# Taking into account the RNA secondary structure.

setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
library(ggplot2)
library(tidyverse)

# Loading the fitness values
tf <- read.csv("fitness-mapped-to-protein.csv")

# Truncating extreme fitness values
tf$w[which(tf$w<0)] <- 0
tf$w[which(tf$w>1.5)] <- 1.5

# Loading the RNA secondary structure information.
ts <- read.table("../Structure/2022_NatCom_Lan/structure.ct", sep = "\t", h = F)
colnames(ts) <- c("position", "base", "prev_index", "next_index", "nbp", "nm")
pos_single <- ts$position[which(ts$nbp==0)]

# Adding the secondary structure information to the fitness table.
tf$paired <- T
tf$paired[which(is.element(tf$position, pos_single)==T)] <- F
length(which(tf$paired == T))
length(which(tf$paired == F))


# Ordering the gene names by position
vPos <- by(tf$position, tf$Name, min)
vPosSorted <- vPos[order(vPos)]
tf$Name <- factor(tf$Name, levels = names(vPosSorted))
tf$nb <- 1

tfr <- tf[which(tf$consequence != "nonsense") , ]
tfr <- tf[which(tf$consequence != "nonsense" & tf$paired == F) , ]
dim(tfr)

# Keep only genes with at least "nbMin" data points
nbMin <- 20

vnb <- by(tfr$nb, tfr$Name, sum)
genesOk <- names(vnb[which(vnb>nbMin)])

tfr <- tfr[which(is.element(tfr$Name, genesOk)==T) , ]

gp <- ggplot(data = tfr, aes(x=Name, y=w, fill=factor(consequence))) +
  geom_boxplot() + 
  labs(fill = "Consequence")
plot(gp)


tgenes <- read.csv("../DB/ensembl_ASM985889v3/genes-location.csv")
tgenes$mean <- NA
tgenes$median <- NA
tobs <- readRDS("../ALL_DATA/all-obs-and-frequencies-combined-and-lethal.rds")
for(i in (1:nrow(tgenes))){
  geneName <- tgenes$geneName[i]
  geneStart <- tgenes$start[i]
  geneEnd <- tgenes$end[i]
  w <- which(tobs$position>=geneStart & tobs$position<=geneEnd)
  tgenes$mean[i] <- mean(tobs$DELTA_B_cov[w])
  tgenes$median[i] <- median(tobs$DELTA_B_cov[w])
}
plot(tgenes$median~factor(tgenes$geneName))

hist(tgenes$median)
rug(x = tgenes$median[which(tgenes$geneName=="S")], col = "red", ticksize = 1.0, lwd = 2)
rug(x = tgenes$median[which(tgenes$geneName=="M")], col = "green", ticksize = 1.0, lwd = 2)
rug(x = tgenes$median[which(tgenes$geneName=="N")], col = "blue", ticksize = 1.0, lwd = 2)


tmp <- tfr[which(tfr$consequence == "nonsynonymous") , ]
vnb <- by(tmp$nb, tmp$Name, sum)
genesOk <- names(vnb[which(vnb>nbMin)])
tmp <- tmp[which(is.element(tmp$Name, genesOk)==T) , ]

t.test(tmp$w[which(tmp$Name == "N")] , tmp$w[which(tmp$Name == "S")])
t.test(tmp$w[which(tmp$Name == "N")] , tmp$w[which(tmp$Name == "M")])
t.test(tmp$w[which(tmp$Name == "M")] , tmp$w[which(tmp$Name == "S")])



tfr <- tf[which(tf$consequence == "nonsynonymous") , ]
vnb <- by(tfr$nb, tfr$Name, sum)
genesOk <- names(vnb[which(vnb>nbMin)])
tfr <- tfr[which(is.element(tfr$Name, genesOk)==T) , ]

gp <- ggplot(data = tfr, aes(x=Name, y=w, fill=factor(paired))) +
  geom_boxplot() + 
  labs(fill = "Paired")
plot(gp)


tmp <- tf[which(tf$Name != "S" & tf$consequence == "nonsynonymous") , ]
tmp <- tf[which(tf$Name != "S" & tf$consequence == "synonymous") , ]
dim(tmp)
t.test(tmp$w[which(tmp$paired==T)] , tmp$w[which(tmp$paired==F)])
