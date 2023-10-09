################################################################################
#
# Generating a table with the coverage at each position of the genome for each 
# sample.

setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")

library(genomation)
library(Biostrings)
library(BSgenome)

library(readxl)

genome_file = "C:/lab/TR_ERRORS/DATA/Coronavirus/DB/ensembl_ASM985889v3/Sars_cov_2.ASM985889v3.dna.toplevel.fa"
gff_file = "C:/lab/TR_ERRORS/DATA/Coronavirus/DB/ensembl_ASM985889v3/Sars_cov_2.ASM985889v3.101.gff3"

BD_DATA = "C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/"

# This file contains the information about where the data from each sample is stored.
library_xls_file <- paste(BD_DATA, "ALL_COVID_LIBRARIES_PATH.xlsx", sep="/")
xls_libs <- read_excel(library_xls_file, sheet = 2)
tlibs <- as.data.frame(xls_libs)


genome = readDNAStringSet(genome_file)
genomeSeq <- genome[1]
ref_seq <- strsplit(as.character(genomeSeq), "")[[1]]

tres <- data.frame(pos = seq(from = 1, to = length(ref_seq), by = 1),
                   ref = ref_seq
                   )

sample_names <- tlibs$true_sampleName[order(tlibs$true_sampleName)]
#sample_names <- sample_names[which(sample_names != "ALPHA_A_1")] # <- Tempoorary fix
# I had been running the ALPHA_A_1 against the human genome to fish out reads from the host cells.
# Need to re-download the data from the actual sars-cov-2 genome...

for( sample_name in sample_names){
  cat("Working on ", sample_name, " ...")
  w <- which(tlibs$true_sampleName == sample_name)
  sample_dir <- tlibs$run[w]
  data_dir <- paste(BD_DATA, sample_dir, sep="/")
  prefix <- tlibs$prefix[w]
  obs_file <- paste(data_dir, "/", prefix, "-obs.tab.gz", sep="")
  tobs <- read.table(obs_file, h=F)
  colnames(tobs) <- c("chromosome", "start", "ref", "strand", "A", "T", "C", "G")
  tobs$cov <- tobs$A + tobs$T + tobs$G + tobs$C
  w <- tobs$start
  cname <- sample_name
  tres[ , cname] <- 0
  tres[w,cname] <- tobs$cov
  cat(" done.\n")
}

write.csv(tres, file = "coverage-all-samples.csv")

library(rtracklayer)
library(GenomicFeatures)
library(Gviz)


options(ucscChromosomeNames=FALSE)
txdb = makeTxDbFromGFF(gff_file)
grt <- GeneRegionTrack(txdb)
chr <- "MN908947.3"
gen <- "MN908947.3"


GENERATE_PDFS = FALSE
if( GENERATE_PDFS == TRUE ){

  for( sample_name in sample_names){
    
    axisTrack <- GenomeAxisTrack()
  
    dtrack <- DataTrack(data = tres[ , sample_name], 
                         start = tres$pos, 
                         end = tres$pos, 
                         chromosome = chr, 
                         genome = gen, 
                         name = "coverage")
  
    plot_file <- paste("coverage-graphs/", sample_name, "-coverage.pdf", sep="")  
    pdf(plot_file, width=16, height = 8)
    plotTracks( 
      list(axisTrack, grt, dtrack), 
      from = 0,
      to = max(tres$pos),
      chromosome = chr, transcriptAnnotation = "symbol", type="polygon",
      main = sample_name)
    dev.off(dev.cur())
  }
  
  
  pdf("coverage-graphs/all-coverages.pdf", width=16, height=8)
  for( sample_name in sample_names){
    
    axisTrack <- GenomeAxisTrack()
    
    dtrack <- DataTrack(data = tres[ , sample_name ], 
                        start = tres$pos, 
                        end = tres$pos, 
                        chromosome = chr, 
                        genome = gen, 
                        name = "coverage")
    
    plotTracks( 
      list(axisTrack, dtrack), 
      from = 0,
      to = max(tres$pos),
      chromosome = chr, type="polygon",
      main = sample_name)
  }
  dev.off(dev.cur())
}

for(sample_name in sample_names){
  nb <- as.numeric(sum(as.numeric(tres[ , sample_name])))
  tres[ , sample_name] <- (tres[ , sample_name]*100) / nb
}

mcor <- matrix(data = 1, ncol = length(sample_names), nrow = length(sample_names))
tcor <- as.data.frame(mcor)
colnames(tcor) <- sample_names
rownames(tcor) <- sample_names

for(i in (1:(length(sample_names)-1)) ){
  sample_name <- sample_names[i]
  for(j in (i+1):length(sample_names)){
    second_sample <- sample_names[j]
    cc <- cor(tres[ , sample_name], tres[ , second_sample])
    tcor[sample_name,second_sample] <- as.numeric(cc)
    tcor[second_sample,sample_name] <- as.numeric(cc)
  }
}

library(ggplot2)
library(reshape2)
library(heatmaply)

library(gplots)

pdf(file="coverage-correlation.pdf", width=12, height=12)
#heatmaply_cor(as.matrix(tcor))
heatmap.2(x=as.matrix(tcor))
dev.off(dev.cur())


melted_cormat <- melt(as.matrix(tcor))
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill = value)) + 
  geom_tile()






sample_1 <- "DELTA_B_1"
sample_2 <- "OMICRON_A_1"
tt <- tres[which(tres[ , sample_1]>0 & tres[ , sample_2]>0) , ]
v <- log(tt[ , sample_1]/tt[ , sample_2])
vabs <- abs(v)
vbig <- vabs[which(vabs>2)]
sum(vbig)

axisTrack <- GenomeAxisTrack()
dtrack <- DataTrack(data = v, 
                    start = tt$pos, 
                    end = tt$pos, 
                    chromosome = chr, 
                    genome = gen, 
                    name = "coverage")

plotTracks( 
  list(axisTrack, grt, dtrack), 
  from = 0,
  to = max(tres$pos),
  chromosome = chr, type="polygon",
  transcriptAnnotation = "symbol", type="polygon",
  main = paste(sample_1, " vs ", sample_2, sep="")
)

