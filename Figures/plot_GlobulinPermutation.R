#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used to create our permutation graph for globulin phenotype enrichment
#------------------------------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)

#Permutation table
GlobulinPermutation <- read.table("PhenotypeGlobulinPerm.txt", header =T)

#95th percentile of observed permutations
highpercentile <- quantile(GlobulinPermutation$x, .95)

#95th percentile and observed gene point
AnthroPermutationPlot <- ggplot(GlobulinPermutation, aes(x=x)) + 
  theme(text=element_text(size=50)) + 
  geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = highpercentile, size = 3, color = "#92c648", linetype = "dashed") + 
  geom_point(x=3, y=0, size = 8, color = "#92c648") +
  theme_classic() + 
  theme(text = element_text(size=20))

ggsave(AnthroPermutationPlot, file="GlobulinPermutation.pdf", height=7, width =7, useDingbats=FALSE)

