#!/usr/bin/env Rscript

#===================================================================================
#--Plot for gene permutation Supplementary Figure 1
#===================================================================================

AnthroPermutationPdiff4 <- read.table("Pdiff4AnthroPermutationTable.tsv", header =T)

#95th percentile of observed permutations
highpercentile <- quantile(AnthroPermutationPdiff4$NumberGenes, .95)

#red 95th percentile and observed gene point
AnthroPermutationPlot <- ggplot(AnthroPermutationPdiff4, aes(x=NumberGenes)) + 
  theme(text=element_text(size=50)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(name = "Number of simulated unique sexual differentiation genes\n per permutation") + 
  scale_y_continuous(name="Count per permutation") + 
  geom_vline(xintercept = highpercentile, size = 2, color = "#92c648", linetype = "dashed") + 
  geom_point(x=9, y=0, size = 10, color = "#92c648") +
  theme_classic() + 
  theme(text = element_text(size=20))

AnthroPermutationPlot

ggsave(AnthroPermutationPlot, file="Pdiff4AnthropometricGenePermutation.pdf", height=7, width =10)
################
AnthroPermutationPdiff3 <- read.table("Pdiff3AnthroPermutationTable.tsv", header =T)

#95th percentile of observed permutations
highpercentile <- quantile(AnthroPermutationPdiff3$NumberGenes, .95)

#red 95th percentile and observed gene point
AnthroPermutationPlot <- ggplot(AnthroPermutationPdiff3, aes(x=NumberGenes)) + 
  theme(text=element_text(size=50)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(name = "Number of simulated unique sexual differentiation genes\n per permutation") + 
  scale_y_continuous(name="Count per permutation") + 
  geom_vline(xintercept = highpercentile, size = 2, color = "#92c648", linetype = "dashed") + 
  geom_point(x=13, y=0, size = 10, color = "#92c648") +
  theme_classic() + 
  theme(text = element_text(size=20))

AnthroPermutationPlot

ggsave(AnthroPermutationPlot, file="Pdiff3AnthropometricGenePermutation.pdf", height=7, width =10)
###############
AnthroPermutationPdiff2 <- read.table("Pdiff2AnthroPermutationTable.tsv", header =T)

#95th percentile of observed permutations
highpercentile <- quantile(AnthroPermutationPdiff2$NumberGenes, .95)

#red 95th percentile and observed gene point
AnthroPermutationPlot <- ggplot(AnthroPermutationPdiff2, aes(x=NumberGenes)) + 
  theme(text=element_text(size=50)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(name = "Number of simulated unique sexual differentiation genes\n per permutation") + 
  scale_y_continuous(name="Count per permutation") + 
  geom_vline(xintercept = highpercentile, size = 2, color = "#92c648", linetype = "dashed") + 
  geom_point(x=16, y=0, size = 10, color = "#92c648") +
  theme_classic() + 
  theme(text = element_text(size=20))

AnthroPermutationPlot

ggsave(AnthroPermutationPlot, file="Pdiff2AnthropometricGenePermutation.pdf", height=7, width =10)
############
AnthroPermutationPdiff1 <- read.table("Pdiff1AnthroPermutationTable.tsv", header =T)

#95th percentile of observed permutations
highpercentile <- quantile(AnthroPermutationPdiff1$NumberGenes, .95)

#red 95th percentile and observed gene point
AnthroPermutationPlot <- ggplot(AnthroPermutationPdiff1, aes(x=NumberGenes)) + 
  theme(text=element_text(size=50)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(name = "Number of simulated unique sexual differentiation genes\n per permutation") + 
  scale_y_continuous(name="Count per permutation") + 
  geom_vline(xintercept = highpercentile, size = 2, color = "#92c648", linetype = "dashed") + 
  geom_point(x=24, y=0, size = 10, color = "#92c648") +
  theme_classic() + 
  theme(text = element_text(size=20))

AnthroPermutationPlot

ggsave(AnthroPermutationPlot, file="Pdiff1AnthropometricGenePermutation.pdf", height=7, width =10)
