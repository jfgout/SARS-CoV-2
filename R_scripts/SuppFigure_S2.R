# Generating the mutation spectrum by passage for lethal and non-lethal mutations

setwd("C:/lab/TR_ERRORS/DATA/Coronavirus/ALL_DATA/")

library(ggplot2)

tsl <- read.csv("spectrum-lethal_COV-1000_RATE-0.001_HS-10.csv")
tsn <- read.csv("spectrum-non-lethal_COV-1000_RATE-0.001_HS-1e+12.csv")

tsl <- tsl[which(tsl$variant=="DELTA" & tsl$replicate=="B") , ]
tsn <- tsn[which(tsn$variant=="DELTA" & tsn$replicate=="B") , ]

tsn$passage = factor(tsn$passage)

gp <- ggplot(data = tsn, aes(x = mutation, y = rate, by = passage)) +
  geom_bar(stat="identity", position=position_dodge(), aes(fill = passage) )
plot(gp)
