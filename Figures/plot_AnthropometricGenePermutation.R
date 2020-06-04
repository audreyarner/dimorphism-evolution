#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used to create our permutation graph for sex differentiation gene enrichment
#------------------------------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

#Permutation table pulling from all anthropometric phenotypes
AnthroPermutationPdiff4 <- read.table("Pdiff4AnthroPermutationTable.tsv", header =T)

#95th percentile of observed permutations
highpercentile <- quantile(AnthroPermutationPdiff4$NumberGenes, .95)

#red 95th percentile and observed gene point
AnthroPermutationPlot <- ggplot(AnthroPermutationPdiff4, aes(x=NumberGenes)) + 
  theme(text=element_text(size=50)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(name = "Number of simulated unique sexual differentiation genes\n per permutation") + 
  scale_y_continuous(name="Count per permutation") + 
  geom_vline(xintercept = highpercentile, size = 3, color = "#92c648", linetype = "dashed") + 
  geom_point(x=9, y=0, size = 15, color = "#92c648") +
  theme_classic() + 
  theme(text = element_text(size=20))

ggsave(AnthroPermutationPlot, file="Pdiff4AnthropometricGenePermutation.pdf", height=7, width =12, useDingbats=FALSE)

