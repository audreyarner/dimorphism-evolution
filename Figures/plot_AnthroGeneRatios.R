#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used to create our ratio graph for sex differentiation gene enrichment
#------------------------------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

#Ratio table pulling from all anthropometric phenotypes (Supplementary Figure 2)
#--table consists of Pdiff level and ratio 
Ratios <- read.csv("EnrichmentWindow7-28.csv", header = TRUE)

#red 95th percentile and observed gene point
Plot <- ggplot(Ratios, aes(x=pdiff, y = ratio, color = Window, group = Window)) +
  geom_point(size = 8) + 
  geom_line(size = 2) + 
  theme_classic(base_size = 20) + 
  geom_hline(yintercept = 1, color = "#98c658", linetype = "dashed", size = 2) + 
  scale_x_discrete(limits = c("Pdiff4", "Pdiff3", "Pdiff2", "Pdiff1"), labels = c(".001", ".005", ".01", ".05"), name = "P-diff adjusted FDR threshold") +
  scale_color_manual(values=c("black"), labels = c("2000 bp")) + 
  theme(legend.position = "none", text=element_text(size=50)) + 
  scale_y_continuous(name="Ratio of observed proportion of sexually differentiated /ngenes for all anthropometric phenotypes to /nobserved proportion of sexually differentiated /ngenes genome wide") +
  ylim(0,3)

Plot

ggsave(Plot, file="AnthroRatioEnrichment7-28.pdf", height=8.5, width =10, useDingbats=FALSE)

