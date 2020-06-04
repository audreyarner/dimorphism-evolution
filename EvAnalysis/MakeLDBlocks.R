library(dplyr)
LDBlocks <- read.table("NEWall_fourier_ls.bed")

LDBlocks <- LDBlocks[order((LDBlocks$V1)),]

LDBlocks <- LDBlocks %>% mutate(Index = row_number())

names(LDBlocks)[1] <- "CHR"
names(LDBlocks)[2] <- "Block_Start"
names(LDBlocks)[3] <- "Block_End"

write.table(LDBlocks, "LDBlocks.txt", sep="\t", row.names = F, quote = F)