################################################################################
# Statistics about lethal vs non-lethal mutations

BD_DATA <- ""
DB_DIR <- ""

if( .Platform$OS.type == "windows" ){
  setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
  source("C:/lab/TR_ERRORS/R_files/seqan-test/seqan-functions-v3.R")
  source("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/R_files/getLethals.R")
  source("R_files/covid-functions-v2.R")
  
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

# This file contains the unfiltered mutation frequencies for every possible mutation in every sample
fname <- paste(BD_DATA, "all-frequencies-combined-and-lethal.tab.gz", sep="/")
tmuts <- read.table(fname, h=T)

tr <- tmuts[ , c("position", "uid", "ref", "bTo", "lethal")]
tr$mutation <- paste(tr$ref, "->", tr$bTo, sep = "")
tr$nb <- 1

vnb_all <- by(tr$nb, tr$mutation, sum)

trl <- tr[which(tr$lethal == T) , ]
vnb_l <- by(trl$nb, trl$mutation, sum)

t1 <- data.frame(
  mutation = names(vnb_all),
  nb_total = as.numeric(vnb_all)
)
t2 <- data.frame(
  mutation = names(vnb_l),
  nb_lethal = as.numeric(vnb_l)
)

tres <- merge(t1, t2, by = "mutation", all.x = T)
tres$ratio <- tres$nb_lethal / tres$nb_total
tres$bFrom <- substr(tres$mutation, 1, 1)
tres$bFrom <- factor(tres$bFrom, levels = c("A", "C", "G", "T"))

pdf( file = "ratio-lethal-non-lethal.pdf", width = 12, height = 8)
gp <- ggplot(data = tres, aes(x = mutation, y = ratio, fill = bFrom)) +
  geom_bar(stat="identity", position=position_dodge() )
plot(gp)
dev.off(dev.cur())


pdf( file = "nb-lethal-mutations-by-type.pdf", width = 12, height = 8)
gp <- ggplot(data = tres, aes(x = mutation, y = nb_lethal, fill = bFrom)) +
  geom_bar(stat="identity", position=position_dodge() )
plot(gp)
dev.off(dev.cur())
