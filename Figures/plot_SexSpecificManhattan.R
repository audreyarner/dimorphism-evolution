#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used to create our female-specific and male-specific manhattan plots
#------------------------------------------------------------------------------------------------------

#code taken from https://github.com/boxiangliu/manhattan/blob/master/R/manhattan.R
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

#my code to make a manhattan plot
#commonly used genome-wide significance line
GenomeWideSig <- -log10(5e-8)

#Manhattan plot for height female
ManhattanWomen <- read.table("HeightFemaleManhattan.tsv", header =T)

#changes names for function
names(ManhattanWomen)[1] <- "rsid"
names(ManhattanWomen)[2] <- "chrom"
names(ManhattanWomen)[3] <- "pos"
names(ManhattanWomen)[4] <- "pval"

#function requires chr1, chr2, etc
ManhattanWomen$chrom <- sub("^", "chr", ManhattanWomen$chrom)

#transform -log10(pval)
ManhattanWomen$y=-log10(ManhattanWomen$pval)

ManhattanWomen$pos <- as.numeric(as.character(ManhattanWomen$pos))  

#make plot
plot <- manhattan(ManhattanWomen, build='hg19', color1 = "#c94f7c", color2 = "#f3cbd3") +
  geom_hline(yintercept = GenomeWideSig, size = 1.6, color = "black") + 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title=element_text(size=18))

ggsave(plot, "FemaleManhattan.pdf", width =10, height=6, useDingbats=FALSE)

#Manhattan plot for height Male  
ManhattanMen <- read.table("HeightMaleManhattan.tsv", header =T)

#changes names for function
names(ManhattanMen)[1] <- "rsid"
names(ManhattanMen)[2] <- "chrom"
names(ManhattanMen)[3] <- "pos"
names(ManhattanMen)[4] <- "pval"
ManhattanMen$chrom <- sub("^", "chr", ManhattanMen$chrom)

ManhattanMen$y=-log10(ManhattanMen$pval)

ManhattanMen$pos <- as.numeric(as.character(ManhattanMen$pos))  

plot <- manhattan(ManhattanMen, build='hg19', color1 = "darkolivegreen4", color2 = "darkolivegreen3") +
  geom_hline(yintercept = GenomeWideSig, size = 1.6, color = "black") + 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title=element_text(size=18))

ggsave(plot, "MaleManhattan.pdf", width =10, height=6, useDingbats=FALSE)

  