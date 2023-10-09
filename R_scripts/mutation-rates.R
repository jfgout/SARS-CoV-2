MIN_COV <- 1000
MAX_RATE <- 1/MIN_COV
MAX_BEFORE_HS = 1000

USE_ONLY_LETHAL <- F
EXCLUDE_C_TO_T <- F

####################################################################################################
# First part: location of files:
setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")
source("../R_files/seqan-functions.R")
source("../R_files/covid-functions.R")

FILE_PROTEIN_PREDICT <- "../DB/protein_predict.rds" # R binary file with the results from predicting mutation consequence (synonymous, non-synonymous, ...)
FILE_ALL_FREQUENCIES <- "all-frequencies.rds"



library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(genomation)
library(Biostrings)

library(readxl)

vBases = getBases()
vMuts = getMutsFromBases(vBases)

# This file contains the unfiltered mutation frequencies for every possible mutation in every sample
tmuts <- readRDS("all-frequencies-combined-and-lethal.rds")
if( EXCLUDE_C_TO_T == T ){
  vMuts <- vMuts[which(vMuts != "C->T")]
  w <- which(tmuts$ref=="C" & tmuts$bTo=="T")
  if( length(w) > 0 ){
    tmuts <- tmuts[-w , ]
  }
}

# This file contains the "high frequency" mutations to remove from the list of lethal ones
fName <- paste("Mutation_higher_than_", MAX_RATE, ".tab", sep = "")
th <- read.table(fName, h = F, as.is = T)
vh <- th$V1

length(which(tmuts$lethal == T))
tmuts$lethal[which( is.element(tmuts$uid, vh) == T )] = FALSE
length(which(tmuts$lethal == T))

# Subsetting only lethal mutations
tl <- tmuts[which(tmuts$lethal == T) , ]
tnl <- tmuts[which(tmuts$lethal == F) , ]
dim(tl)
dim(tnl)

# Extracting the list of sample names from the columns names
vcn <- colnames(tmuts)
wmut <- grepl("_muts", vcn, fixed = TRUE)
vcr <- vcn[which(wmut==T)]
sampleNames <- substr(vcr, 1, nchar(vcr)-5)
#sampleNames <- sampleNames[which(sampleNames != "total")]

tSpecs_lethal <- computeSpectrumFromMuts(
  tmut = tl,
  sampleNames = sampleNames,
  MIN_COV = 1,
  MAX_RATE = 0.1,
  MAX_BEFORE_HS = MAX_BEFORE_HS,
  bFrom = "ref",
  bTo = "bTo",
  vMuts = vMuts,
)

tSpecs_non_lethal <- computeSpectrumFromMuts(
  tmut = tnl,
  sampleNames = sampleNames,
  MIN_COV = MIN_COV,
  MAX_RATE = MAX_RATE,
  MAX_BEFORE_HS = MAX_BEFORE_HS,
  bFrom = "ref",
  bTo = "bTo",
  vMuts = vMuts,
)

tSpecs_lethal$lethal <- TRUE
tSpecs_non_lethal$lethal <- FALSE
tSpecs <- tSpecs_lethal
if( USE_ONLY_LETHAL == F ){
  tSpecs <- rbind(tSpecs_lethal, tSpecs_non_lethal)
}

tSpecs <- decomposeVariantReplicatePassage(tSpecs)
tSpecs$series <- paste(tSpecs$variant, tSpecs$replicate, sep = "_")
wFullVariant <- which(is.na(tSpecs$replicate)==T)
tSpecs$series[wFullVariant] <- tSpecs$variant[wFullVariant]

unique(tSpecs$series)

tres <- tSpecs
tres$rate <- tres$nbMut / tres$nbCalls

# Computing the confidence intervals:
for(i in (1:nrow(tres))){
  
  if( is.na(tres$passage[i]) == F ){
    # If I'm looking at a unique sample, I compute the CI with the confidence interval of the fraction
    nb <- tres$nbMut[i]
    cov <- tres$nbCalls[i]
    pp <- prop.test(x=nb, n=cov)
    tres$ci_low[i] <- as.numeric(pp$conf.int[1])
    tres$ci_high[i] <- as.numeric(pp$conf.int[2])
  } else {
    # I'm looking at a group of samples either a series or a full variant -> using sd for the confidence interval
    ww <- NA
    if( is.na(tres$replicate[i]) == T ){
      # Looking at the full variant (= replicates combined)
      ww <- which(tres$variant == tres$variant[i] 
                  & is.na(tres$replicate) == F 
                  & tres$mutation == tres$mutation[i]
                  & tres$lethal == tres$lethal[i]
                    )
    } else {
      # I'm looking at a series of passages (= only one replicate).
      ww <- which(tres$variant == tres$variant[i] 
                  & is.na(tres$replicate) == F 
                  & tres$replicate == tres$replicate[i]
                  & tres$mutation == tres$mutation[i]
                  & tres$lethal == tres$lethal[i]
      )
    }
    
    # Calculation is the same for both:
    sd_val <- sd(tres$rate[ww])
    tres$ci_high[i] <- tres$rate[i] + sd_val/2
    tres$ci_low[i] <- tres$rate[i] - sd_val/2
    
  }

}


fName <- paste("spectrum-lethal-COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
if( EXCLUDE_C_TO_T == T ){
  fName <- paste("spectrum-lethal-no_C_to_T-COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
}
if( USE_ONLY_LETHAL == F ){
  fName <- gsub("lethal", "non-lethal", fName)
}

write.csv(tres, file = fName, row.names = F)


################################################################################
# Run this second part to generate the graphs after having already done all the 
# calculations in the previous half of the script.

fname <- paste("spectrum-COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
if( EXCLUDE_C_TO_T == T ){
  fName <- paste("spectrum-COV-no_C_to_T-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
}

tres <- read.csv(fname)
tres$mutation <- factor(tres$mutation, levels = c(vMuts, "N->N"))

PLOT_GRAPHS <- T

fpdf <- paste("spectrum-lethal_COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".pdf", sep="")
if( EXCLUDE_C_TO_T == T ){
  fpdf <- paste("spectrum-lethal_COV_no_C_to_T-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".pdf", sep="")
}
LETHAL_STATUS_TO_PLOT <- TRUE
if( LETHAL_STATUS_TO_PLOT == FALSE ){
  fpdf <- gsub("lethal", "non-lethal", fpdf)
}


if( PLOT_GRAPHS == T ){ pdf(file = fpdf, width = 12, height = 8) }

# Let's start by plotting the graph with all the populous series:
tr <- tres[which(is.element(tres$variant, c("USA", "ALPHA", "DELTA"))==T & is.na(tres$passage)==T & is.na(tres$replicate)==F & tres$lethal==LETHAL_STATUS_TO_PLOT) , ]
gp <- ggplot( data = tr, aes(x = mutation, y = rate, fill = series) ) + 
  geom_bar(stat="identity", position=position_dodge() ) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.2, position=position_dodge(.9), color = "grey") +
  ggtitle("Mutation spectrum") +
  theme(plot.title = element_text(hjust = 0.5))
plot(gp)


# Next step is to plot all the series individually (= all the passages of USA_A in one graph, all the passage of DELTA_B, ...)
for(series in unique(tres$series)){
  tr <- tres[which(tres$series==series & is.na(tres$passage)==F & tres$lethal == LETHAL_STATUS_TO_PLOT) , ]
  tr$passage <- factor(tr$passage)
  gp <- ggplot( data = tr, aes(x = mutation, y = rate, fill = passage) ) + 
    geom_bar(stat="identity", position=position_dodge() ) +
    geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.2, position=position_dodge(.9), color = "grey") +
    ggtitle(series) +
    theme(plot.title = element_text(hjust = 0.5))
  plot(gp)
  
}

# Then, I can plot each "sample" (including combined series/variants) separately:
for(sampleName in unique(tres$sampleName)){
  tr <- tres[which(tres$sampleName == sampleName & tres$lethal == LETHAL_STATUS_TO_PLOT) , ]
  
  gp <- ggplot( data = tr, aes(x = mutation, y = rate) ) + 
      geom_bar(stat="identity", position=position_dodge() ) +
      geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.2, position=position_dodge(.9), color = "grey") +
      ggtitle(sampleName) +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot(gp)
}
if( PLOT_GRAPHS == T ){ dev.off(dev.cur()) }

png(file = "spectrum-populous-series.png", width = 800, height = 400)
tr <- tres[which(is.element(tres$variant, c("USA", "ALPHA", "DELTA"))==T & is.na(tres$passage)==T & is.na(tres$replicate)==F & tres$lethal==LETHAL_STATUS_TO_PLOT) , ]
gp <- ggplot( data = tr, aes(x = mutation, y = rate, fill = series) ) + 
  geom_bar(stat="identity", position=position_dodge() ) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.2, position=position_dodge(.9), color = "grey") +
  ggtitle("Mutation spectrum") +
  theme(plot.title = element_text(hjust = 0.5))
plot(gp)
dev.off(dev.cur())




################################################################################
# Comparing the frequencies of lethal and non-lethal mutations:

fname_lethal <- paste("spectrum-lethal_COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
fname_non_lethal <- paste("spectrum-non-lethal_COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")

tsl <- read.csv(fname_lethal)
tsn <- read.csv(fname_non_lethal)

tsl$lethal = T
tsn$lethal = F
tsf <- rbind(tsn, tsl)
tsf$lethal <- factor(tsf$lethal)

tsf$series <- substr(tsf$sampleName, 1, nchar(tsf$sampleName)-2)
tsf$passage <- substr(tsf$sampleName, nchar(tsf$sampleName), nchar(tsf$sampleName))
ts <- tsf[which(tsf$series == "DELTA_B") , ]

fname <- paste("spectrum-lethal_vs_nonlethal_COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".csv", sep="")
write.csv(ts, file = fname)

fname_pdf <- paste("spectrum-lethal-vs-non-lethal_DELTA_B_COV-", MIN_COV, "_RATE-", MAX_RATE, "_HS-", MAX_BEFORE_HS, ".pdf", sep="")
#fname_pdf <- "tmp.pdf"
pdf(file = fname_pdf, width = 12, height = 8)

for(passage in unique(ts$passage)) {
  tsr <- ts[which(ts$passage == passage) , ]
  title <- paste("Mutation frequency - lethal vs non-lethal - passage ", passage, sep = "")
  gp <- ggplot( data = tsr, aes(x = mutation, y = rate, fill = lethal) ) + 
    geom_bar(stat="identity", position=position_dodge() ) +
    #geom_errorbar(aes(ymin=ci_min, ymax=ci_max), width=.2, position=position_dodge(.9), color = "grey") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot(gp)
  
}
dev.off( dev.cur() )
