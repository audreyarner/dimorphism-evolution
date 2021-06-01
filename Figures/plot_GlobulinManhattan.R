library(dplyr)
library(ggplot2)

get_chrom_lengths=function(build=c('hg18','hg19','hg38')){
  build=match.arg(build)
  
  chrom_lengths_hg18=c(chr1=247249719,chr2=242951149,chr3=199501827,
                       chr4=191273063,chr5=180857866,chr6=170899992,
                       chr7=158821424,chr8=146274826,chr9=140273252,
                       chr10=135374737,chr11=134452384,chr12=132349534,
                       chr13=114142980,chr14=106368585,chr15=100338915,
                       chr16=88827254,chr17=78774742,chr18=76117153,
                       chr19=63811651,chr20=62435964,chr21=46944323,
                       chr22=49691432)
  
  
  chrom_lengths_hg19=c(chr1=249250621,chr2=243199373,chr3=198022430,
                       chr4=191154276,chr5=180915260,chr6=171115067,
                       chr7=159138663,chr8=146364022,chr9=141213431,
                       chr10=135534747,chr11=135006516,chr12=133851895,
                       chr13=115169878,chr14=107349540,chr15=102531392,
                       chr16=90354753,chr17=81195210,chr18=78077248,
                       chr19=59128983,chr20=63025520,chr21=48129895,
                       chr22=51304566)
  
  # https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
  chrom_lengths_hg38=c(chr1=248956422,chr2=242193529,chr3=198295559,
                       chr4=190214555,chr5=181538259,chr6=170805979,
                       chr7=159345973,chr8=145138636,chr9=138394717,
                       chr10=133797422,chr11=135086622,chr12=133275309,
                       chr13=114364328,chr14=107043718,chr15=101991189,
                       chr16=90338345,chr17=83257441,chr18=80373285,
                       chr19=58617616,chr20=64444167,chr21=46709983,
                       chr22=50818468)
  
  if (build=='hg18'){
    return(chrom_lengths_hg18)
  } else if (build=='hg19'){
    return(chrom_lengths_hg19)
  } else {
    return(chrom_lengths_hg38)
  }
}

get_cumulative_length=function(chrom_lengths){
  n_chrom=length(chrom_lengths)
  cumulative_length=c(0,cumsum(x = unname(chrom_lengths))[1:(n_chrom-1)])
  names(cumulative_length)=names(chrom_lengths)
  return(cumulative_length)
}

get_total_length=function(chrom_lengths){
  sum(chrom_lengths)
}

get_x_breaks=function(chrom_lengths){
  cumulative_length=get_cumulative_length(chrom_lengths)
  x_breaks=cumulative_length+round(chrom_lengths/2)
  names(x_breaks)=gsub('chr','',names(x_breaks))
  names(x_breaks)[20]=''
  names(x_breaks)[22]=''
  return(x_breaks)
}


add_cumulative_pos=function(data,build=c('hg18','hg19','hg38')){
  build=match.arg(build)
  
  if (!all(c('chrom','pos','y')%in%colnames(data))){
    stop('data must have chrom, pos and y columns')
  }
  
  chrom_lengths=get_chrom_lengths(build)
  cumulative_length=get_cumulative_length(chrom_lengths)
  
  calc_cumulative_pos=function(chrom,pos,cumulative_length){
    stopifnot(length(unique(chrom))==1)
    stopifnot(chrom%in%names(cumulative_length))
    pos+cumulative_length[chrom]
  }
  
  for (i in paste0('chr',seq(1,22))){
    data[data$chrom==i,'cumulative_pos']=data[data$chrom==i,'pos']+cumulative_length[i]
  }
  
  return(data)
}

add_color=function(data,color1='black',color2='grey'){
  if ('color'%in%colnames(data)){
    user_color=data$color
  } else {
    user_color=rep(NA,nrow(data))
  }
  
  data$color=ifelse(data$chrom%in%paste0('chr',seq(1,22,2)),color1,color2)
  data$color=ifelse(is.na(user_color),data$color,user_color)
  return(data)
}

add_shape = function(data,shape=16){
  if ('shape'%in% colnames(data)){
    user_shape = data$shape
  } else {
    user_shape = rep(NA, nrow(data))
  }
  
  data$shape = shape
  data$shape = ifelse(is.na(user_shape),data$shape,user_shape)
  return(data)
}

add_fill = function(data){
  if ('fill'%in%colnames(data)){
    NULL
  } else {
    data$fill = data$color
  } 
  return(data)
}

#' Make a Manhattan plot
#' @param gwas A GWAS study. Must have chrom, pos, and y columns
#' @param build Genomic build. Currently supports hg18, hg19, and hg38
#' @param color1 Color for odd-numbered chromosomes
#' @param color2 Color for even-numbered chromosomes
#' @example
#' data(cad_gwas)
#' cad_gwas$y=-log10(cad_gwas$pval)
#' manhattan(cad_gwas,build='hg18')
#' @return a ggplot object that makes a Manhattan plot
#' @details
#' \code{manhattan} is a wrapper around \code{ggplot}. It uses a few tricks to transform a genomic axis to a scatterplot axis. For instance, chr2:1 would be the length of chromosome 1 plus 1, chr3:1 would be chromosome 1 plus chromosome 2 plus 1, so on and so forth. It is important to specify the genomic build (e.g. hg19) so that `manhattan` can make the correct transformation. It positions the chromosome labels on the x-axis according to these transformations.
manhattan=function(gwas,build=c('hg18','hg19','hg38'),color1='black',color2='grey'){
  data=gwas
  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  data=add_color(data,color1 = color1,color2 = color2)
  data=add_shape(data,shape=16)
  data=add_fill(data)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)
  
  color_map=unique(data$color)
  names(color_map)=unique(data$color)
  
  ggplot2::ggplot(data,aes(x=cumulative_pos,y=y,color=color,shape=shape,fill=fill))+
    geom_point()+
    theme_classic()+
    scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                       labels=names(x_breaks),name='Chromosome')+
    scale_y_continuous(expand=c(0.01,0),name=expression('-log10(P-value)'))+
    scale_color_manual(values=color_map,guide='none')+
    scale_fill_manual(values = color_map, guide = 'none')+
    scale_shape_identity()
}
#------------------------------------------------------------------------------------------------------
#End of code from GitHub
#------------------------------------------------------------------------------------------------------

#LD Blocks for pruning SNPs
LDBlocks <- read.table("LDBlocks.txt", header =T)

#read both SexDiff and phenotype-associated pruned SNPs for each fo the five phenotypes
Height <- read.table("HeightPdiffPrunedSNPs.txt", header =T)
HeightAll <- read.table("HeightPrunedSNPs.txt", header =T)
Height <- Height %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)
HeightAll <- HeightAll %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)

HeightCombine <- merge(Height,HeightAll, all=TRUE)
names(HeightCombine)[4] <- "pdiff.Height"
names(HeightCombine)[5] <- "pdiffFDR.Height"
HeightCombine$pdiff.Height <- as.numeric(HeightCombine$pdiff.Height)

BodyFatPercent <- read.table("BodyFatPercentPdiffPrunedSNPs.txt", header =T)
BodyFatPercentAll <- read.table("BodyFatPercentPrunedSNPs.txt", header =T)
BodyFatPercent <- BodyFatPercent %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)
BodyFatPercentAll <- BodyFatPercentAll %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)

BodyFatPercentCombine <- merge(BodyFatPercent,BodyFatPercentAll, all=TRUE)
names(BodyFatPercentCombine)[4] <- "pdiff.BodyFatPercent"
names(BodyFatPercentCombine)[5] <- "pdiffFDR.BodyFatPercent"
BodyFatPercentCombine$pdiff.BodyFatPercent <- as.numeric(BodyFatPercentCombine$pdiff.BodyFatPercent)

BodyMass <- read.table("BodyMassPdiffPrunedSNPs.txt", header =T)
BodyMassAll <- read.table("BodyMassPrunedSNPs.txt", header =T)
BodyMass <- BodyMass %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)
BodyMassAll <- BodyMassAll %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)

BodyMassCombine <- merge(BodyMass,BodyMassAll, all=TRUE)
names(BodyMassCombine)[4] <- "pdiff.BodyMass"
names(BodyMassCombine)[5] <- "pdiffFDR.BodyMass"
BodyMassCombine$pdiff.BodyMass <- as.numeric(BodyMassCombine$pdiff.BodyMass)


HipCirc <- read.table("HipCircPdiffPrunedSNPs.txt", header =T)
HipCircAll <- read.table("HipCircPrunedSNPs.txt", header =T)
HipCirc <- HipCirc %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)
HipCircAll <- HipCircAll %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)

HipCircCombine <- merge(HipCirc,HipCircAll, all=TRUE)
names(HipCircCombine)[4] <- "pdiff.HipCirc"
names(HipCircCombine)[5] <- "pdiffFDR.HipCirc"
HipCircCombine$pdiff.HipCirc <- as.numeric(HipCircCombine$pdiff.HipCirc)

WaistCirc <- read.table("WaistCircPdiffPrunedSNPs.txt", header =T)
WaistCircAll <- read.table("WaistCircPrunedSNPs.txt", header =T)
WaistCirc <- WaistCirc %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)
WaistCircAll <- WaistCircAll %>% select(rsid, CHR, POS, pdiff, pdiff.FDR)

WaistCircCombine <- merge(WaistCirc,WaistCircAll, all=TRUE)
names(WaistCircCombine)[4] <- "pdiff.WaistCirc"
names(WaistCircCombine)[5] <- "pdiffFDR.WaistCirc"
WaistCircCombine$pdiff.WaistCirc <- as.numeric(WaistCircCombine$pdiff.WaistCirc)

#merge all lists of SexDiff and phenotype-associated SNPs
All <- merge(HeightCombine,BodyMassCombine,by=c("rsid","CHR","POS"),all=TRUE)
All <- merge(All,WaistCircCombine,by=c("rsid","CHR","POS"),all=TRUE)
All <- merge(All,HipCircCombine,by=c("rsid","CHR","POS"),all=TRUE)
All <- merge(All,BodyFatPercentCombine,by=c("rsid","CHR","POS"),all=TRUE)
All[is.na(All)] <- 1

#Identify the most significant P-SexDiff value among all SexDiff-associated variants
All$Pdiff <- NA
for (i in 1:nrow(All)){
  if(All$pdiff.Height[i]<All$pdiff.BodyMass[i]){
    All$Pdiff[i] = All$pdiff.Height[i]
  } else {All$Pdiff[i] = All$pdiff.BodyMass[i]}
}

for (i in 1:nrow(All)){
  if(All$pdiff.BodyFatPercent[i]<All$Pdiff[i]){
    All$Pdiff[i] = All$pdiff.BodyFatPercent[i]
  }
}

for (i in 1:nrow(All)){
  if(All$pdiff.HipCirc[i]<All$Pdiff[i]){
    All$Pdiff[i] = All$pdiff.HipCirc[i]
  }
}

for (i in 1:nrow(All)){
  if(All$pdiff.WaistCirc[i]<All$Pdiff[i]){
    All$Pdiff[i] = All$pdiff.WaistCirc[i]
  }
}

#Identify the most significant P-SexDiff FDR value among all SexDiff-associated variants
All$pdiffFDR <- NA
for (i in 1:nrow(All)){
  if(All$pdiffFDR.Height[i]<All$pdiffFDR.BodyMass[i]){
    All$pdiffFDR[i] = All$pdiffFDR.Height[i]
  } else {All$pdiffFDR[i] = All$pdiffFDR.BodyMass[i]}
}

for (i in 1:nrow(All)){
  if(All$pdiffFDR.BodyFatPercent[i]<All$pdiffFDR[i]){
    All$pdiffFDR[i] = All$pdiffFDR.BodyFatPercent[i]
  }
}

for (i in 1:nrow(All)){
  if(All$pdiffFDR.HipCirc[i]<All$pdiffFDR[i]){
    All$pdiffFDR[i] = All$pdiffFDR.HipCirc[i]
  }
}

for (i in 1:nrow(All)){
  if(All$pdiffFDR.WaistCirc[i]<All$pdiffFDR[i]){
    All$pdiffFDR[i] = All$pdiffFDR.WaistCirc[i]
  }
}

#Put the desired color depending on whether SNP is significant (P<0.001) or not
All$color <- NA
for (i in 1:nrow(All)){
  if(All$pdiffFDR[i]<0.001){
    All$color[i] = "#F6D000"
  } else {All$color[i] = "#F58EDD"}
}

All$mer <- NA
for (i in 1:nrow(All)){
  if(All$pdiffFDR[i]<0.001){
    All$mer[i] = 1
  } else {All$mer[i] = 2}
}

All$CHR <- as.numeric(All$CHR)
LDBlocks$CHR <- as.numeric(LDBlocks$CHR)

All$POS <- as.numeric(All$POS)
LDBlocks$Block_Start <- as.numeric(LDBlocks$Block_Start)
LDBlocks$Block_End <- as.numeric(LDBlocks$Block_End)

#Prune SNPs so we only get one per LD block, which will be the SNP with the highest P-SexDiff FDR
PrunedSNPsOutput <- rep('a', nrow(LDBlocks))

for(i in 1:nrow(LDBlocks)){
  BlockSNPs <- All %>% filter((LDBlocks$CHR[i] == All$CHR) & (LDBlocks$Block_Start[i] <= All$POS) & (LDBlocks$Block_End[i] >= All$POS))
  OrderbyP <- BlockSNPs[order(BlockSNPs$pdiffFDR),]
  OrderbyP$rsid <- as.character(OrderbyP$rsid)
  PrunedSNPsOutput[i] <- OrderbyP$rsid[1]
}

PrunedSigSNPs <- All[All$rsid %in% PrunedSNPsOutput,]

PrunedSigSNPs[,"Index"] <- NA

for (i in 1:nrow(PrunedSigSNPs)){
  SNP.Overlap = which((LDBlocks$CHR == PrunedSigSNPs$CHR[i]) & (LDBlocks$Block_Start <= PrunedSigSNPs$POS[i]) &(LDBlocks$Block_End>=PrunedSigSNPs$POS[i]))
  if(length(SNP.Overlap>=1)){
    PrunedSigSNPs$Index[i] = LDBlocks$Index[SNP.Overlap]
  } else {SNPs$GeneOverlap[i] = 0}
} 

#read in the SNP lists that have the different GWAS Catalog phenotypes each SNP is associated with
Phenotype <- read.delim("PhenotypeSNPsPleiotropy.txt", header =T, sep="\t")
SexDiff <- read.delim("SexDiffSNPsPleiotropy.txt", header =T, sep="\t")

#identify SexDiff-associated SNPs associated with Sex hormone-binding globulin
for(i in 1:nrow(SexDiff)){ 
  value <- "Sex hormone-binding globulin"
  pleiotropy <- SexDiff$Phenotypes[i]
  if(grepl(value,pleiotropy,fixed=TRUE)){
    SexDiff$Globulin[i] = "Y"
  } else {SexDiff$Globulin[i] = "N"}
}

#if SNP is associated with Sex hormone-binding globulin, it will be a different color in the plot
SexDiff <- SexDiff %>% filter(Globulin == "Y")

for(i in 1:nrow(PrunedSigSNPs)){
  SNP.Overlap = which(PrunedSigSNPs$rsid[i]==SexDiff$rsid)
  if(length(SNP.Overlap>=1)){
    PrunedSigSNPs$color[i] = "#000000"
  }
}

#repeat for phenotype-associated SNPs
for(i in 1:nrow(Phenotype)){ 
  value <- "Sex hormone-binding globulin"
  pleiotropy <- Phenotype$Phenotypes[i]
  if(grepl(value,pleiotropy,fixed=TRUE)){
    Phenotype$Globulin[i] = "Y"
  } else {Phenotype$Globulin[i] = "N"}
}

Phenotype <- Phenotype %>% filter(Globulin == "Y")

for(i in 1:nrow(PrunedSigSNPs)){
  SNP.Overlap = which(PrunedSigSNPs$rsid[i]==Phenotype$rsid)
  if(length(SNP.Overlap>=1)){
    PrunedSigSNPs$color[i] = "#000000"
  }
}

#Have new column "mer" where SNPs associated with Sex hormone-binding globulin will be plotted last
for(i in 1:nrow(PrunedSigSNPs)){
  SNP.Overlap = which(PrunedSigSNPs$rsid[i]==Phenotype$rsid)
  if(length(SNP.Overlap>=1)){
    PrunedSigSNPs$mer[i] = 3
  }
}

PrunedSNPsOutput$mer <- as.numeric(PrunedSNPsOutput$mer)

#change names of columns so runs properly with GitHub Manhattan plot code
names(PrunedSigSNPs)[1] <- "rsid"
names(PrunedSigSNPs)[2] <- "chrom"
names(PrunedSigSNPs)[3] <- "pos"

#formatting data frame to fit within Github Manhattan plot code parameters
PrunedSigSNPs$chrom <- sub("^", "chr", PrunedSigSNPs$chrom)

PrunedSigSNPs$y=-log10(PrunedSigSNPs$Pdiff)

PrunedSigSNPs$pos <- as.numeric(as.character(PrunedSigSNPs$pos))  

#order data so SNPs associated with Sex hormone-binding globulin are plotted last
PrunedSigSNPs2 <- PrunedSigSNPs[with(PrunedSigSNPs, order(mer)),]

p <- manhattan(PrunedSigSNPs2, build='hg19') + 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title=element_text(size=18)) +
  geom_point(size=3.5) +
  ylim(0,16)

p

ggsave(p, file="GlobulinManhattan3-3.pdf", width=10, height=6, useDingbats=FALSE)

