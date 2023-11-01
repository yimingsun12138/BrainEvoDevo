#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: calculate TF binding affinity influenced by human SNC           ##
## Data: 2022.11.08                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#general setting
setwd('/content/data/sunym/project/Brain/')
.libPaths('/content/data/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scibet)
library(Matrix)
library(tidyverse)
library(cowplot)
library(viridis)
library(ComplexHeatmap)
library(parallel)
library(ggsignif)
library(RColorBrewer)
library(ggsci)
library(scales)
library(patchwork)
library(ggpointdensity)
library(latex2exp)
library(ArchR)
library(scales)
library(circlize)
library(ggpubr)
library(ggtext)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(universalmotif)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

#notice:
#Remember that the sequence filter in step_72 is wrong!

# load motif ----------------------------------------------------------
meta_data <- read.csv(file = './data/reference/motif_DB/CISBP/human_221108/TF_information.csv')

#keep only direct motif
table(meta_data$TF_Status)
meta_data <- meta_data[meta_data$TF_Status == 'D',]
motif_list <- unique(meta_data$Motif_ID)

#read in motif files
motif_file <- list()
for (i in motif_list) {
  temp <- list.files(path = './data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
  temp <- temp[grep(pattern = i,x = temp,fixed = TRUE)]
  if(length(temp) != 1){
    print('motif num error!')
    break
  }
  temp <- paste('./data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',temp,sep = '/')
  char <- paste0('cat ',temp,' | wc -l')
  char <- system(command = char,intern = TRUE)
  char <- as.integer(x = char)
  if(char > 1){
    temp <- read_cisbp(file = temp)
    temp@name <- i
    if(length(motif_file) == 0){
      motif_file <- list(temp)
    }else{
      motif_file <- append(x = motif_file,values = temp)
    }
  }
}

#save meme data
write_meme(motifs = motif_file,file = './res/step_77_fig_221108/CISBP_human.meme',overwrite = TRUE)

# extent human sequence ---------------------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_human_SNC.rds')

extented_human_SNC <- human_SNC$human_SNC
extented_human_SNC <- base::do.call(what = rbind,args = strsplit(x = extented_human_SNC,split = '-'))
extented_human_SNC <- as.data.frame(extented_human_SNC)
colnames(extented_human_SNC) <- c('chrom','start')

#extent
extented_human_SNC$start <- as.numeric(extented_human_SNC$start)
extented_human_SNC$end <- extented_human_SNC$start + 30
extented_human_SNC$start <- extented_human_SNC$start - 30
extented_human_SNC[extented_human_SNC$start < 0,"start"] <- 0

extented_human_SNC <- as(extented_human_SNC,'GRanges')
table(extented_human_SNC@ranges@width)

#get sequence
human_seq <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,names = extented_human_SNC,as.character = FALSE)
names(human_seq) <- human_SNC$human_SNC

# extent macaque seq ------------------------------------------------------
#load data
macaque_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_macaque_SNC.rds')
names(macaque_SNC) <- macaque_SNC$human_SNC
macaque_SNC <- macaque_SNC[human_SNC$human_SNC]

extented_macaque_SNC <- macaque_SNC$macaque_SNC
extented_macaque_SNC <- base::do.call(what = rbind,args = strsplit(x = extented_macaque_SNC,split = '-'))
extented_macaque_SNC <- as.data.frame(extented_macaque_SNC)
colnames(extented_macaque_SNC) <- c('chrom','start')

#extent
extented_macaque_SNC$start <- as.numeric(extented_macaque_SNC$start)
extented_macaque_SNC$end <- extented_macaque_SNC$start + 30
extented_macaque_SNC$start <- extented_macaque_SNC$start - 30
extented_macaque_SNC[extented_macaque_SNC$start < 0,"start"] <- 0

extented_macaque_SNC <- as(extented_macaque_SNC,'GRanges')
table(extented_macaque_SNC@ranges@width)

#get seq
macaque_seq <- BSgenome::getSeq(x = BSgenome.Mmulatta.UCSC.rheMac10,names = extented_macaque_SNC,as.character = FALSE)
names(macaque_SNC) <- macaque_SNC$human_SNC

#compare sequence similarity
#create NSM
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

pairwiseAlignment(pattern = human_seq[3],subject = reverseComplement(macaque_seq[3]),type="global",
                  substitutionMatrix = NSM,gapOpening=10, gapExtension=4,
                  scoreOnly = FALSE)

# get paired human macaque sequence ---------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC/extented_human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC/extented_macaque_SNC.rds')

names(human_SNC) <- human_SNC$human_SNC
names(macaque_SNC) <- macaque_SNC$human_SNC

#create NSM
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq_list <- c()
macaque_seq_list <- c()
alignment_score_list <- c()

pb <- txtProgressBar(style = 3)
for (i in human_SNC$human_SNC) {
  #get human seq
  extented_human_SNC <- i
  extented_human_SNC <- strsplit(x = extented_human_SNC,split = '-')
  temp_start <- as.integer(extented_human_SNC[[1]][2]) - 30
  if(temp_start < 0){
    temp_start <- 0
  }
  temp_end <- as.integer(extented_human_SNC[[1]][2]) + 30
  extented_human_SNC <- paste0(extented_human_SNC[[1]][1],':',temp_start,'-',temp_end)
  extented_human_SNC <- as(extented_human_SNC,'GRanges')
  human_seq <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,names = extented_human_SNC,as.character = FALSE)
  
  if(length(human_seq_list) == 0){
    human_seq_list <- human_seq
  }else{
    human_seq_list <- append(x = human_seq_list,values = human_seq)
  }
  
  #get macaque seq
  extented_macaque_SNC <- macaque_SNC[i]$macaque_SNC
  extented_macaque_SNC <- strsplit(x = extented_macaque_SNC,split = '-')
  temp_start <- as.integer(extented_macaque_SNC[[1]][2]) - 30
  if(temp_start < 0){
    temp_start <- 0
  }
  temp_end <- as.integer(extented_macaque_SNC[[1]][2]) + 30
  extented_macaque_SNC <- paste0(extented_macaque_SNC[[1]][1],':',temp_start,'-',temp_end)
  extented_macaque_SNC <- as(extented_macaque_SNC,'GRanges')
  
  macaque_seq <- BSgenome::getSeq(x = BSgenome.Mmulatta.UCSC.rheMac10,names = extented_macaque_SNC)
  
  #find the perfect match
  macaque_seq <- list(macaque_seq,reverse(macaque_seq),complement(macaque_seq),reverseComplement(macaque_seq))
  temp_norm <- pairwiseAlignment(pattern = human_seq,subject = macaque_seq[[1]],type="global",
                                 substitutionMatrix = NSM,gapOpening=10, gapExtension=4,
                                 scoreOnly = TRUE)
  temp_rev <- pairwiseAlignment(pattern = human_seq,subject = macaque_seq[[2]],type="global",
                                substitutionMatrix = NSM,gapOpening=10, gapExtension=4,
                                scoreOnly = TRUE)
  temp_com <- pairwiseAlignment(pattern = human_seq,subject = macaque_seq[[3]],type="global",
                                substitutionMatrix = NSM,gapOpening=10, gapExtension=4,
                                scoreOnly = TRUE)
  temp_comrev <- pairwiseAlignment(pattern = human_seq,subject = macaque_seq[[4]],type="global",
                                   substitutionMatrix = NSM,gapOpening=10, gapExtension=4,
                                   scoreOnly = TRUE)
  
  idx <- which.max(c(temp_norm,temp_rev,temp_com,temp_comrev))
  
  #return
  if(length(macaque_seq_list) == 0){
    macaque_seq_list <- macaque_seq[[idx]]
  }else{
    macaque_seq_list <- append(x = macaque_seq_list,values = macaque_seq[[idx]])
  }
  
  if(length(alignment_score_list) == 0){
    alignment_score_list <- max(c(temp_norm,temp_rev,temp_com,temp_comrev))
  }else{
    alignment_score_list <- append(x = alignment_score_list,values = max(c(temp_norm,temp_rev,temp_com,temp_comrev)))
  }
  setTxtProgressBar(pb,which(human_SNC$human_SNC == i)/length(human_SNC))
}
close(pb)
my_send_sms('get seq done!')

#rename and save it
names(human_seq_list) <- human_SNC$human_SNC
names(macaque_seq_list) <- human_SNC$human_SNC
names(alignment_score_list) <- human_SNC$human_SNC

temp <- SimpleList(human_seq = human_seq_list,
                   macaque_seq = macaque_seq_list,
                   alignment_score = alignment_score_list)

saveRDS(object = temp,file = './res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')

#filter sequence and save it
sequence_list <- readRDS(file = './res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')

temp_human <- unlist(base::lapply(X = human_SNC$human_SNC,FUN = function(x){
  temp <- as.character(sequence_list$human_seq[[x]][31])
  return(temp)
}))

temp_macaque <- unlist(base::lapply(X = human_SNC$human_SNC,FUN = function(x){
  temp <- as.character(sequence_list$macaque_seq[[x]][31])
  return(temp)
}))

idx <- human_SNC$human_SNC[which(temp_human != temp_macaque)]

human_seq_list <- sequence_list$human_seq[idx]
macaque_seq_list <- sequence_list$macaque_seq[idx]

writeXStringSet(x = human_seq_list,filepath = './res/step_77_fig_221108/human_SNC_seq.fasta')
writeXStringSet(x = macaque_seq_list,filepath = './res/step_77_fig_221108/macaque_SNC_seq.fasta')

# see what we get from fimo -----------------------------------------------
#load fimo_out
human_fimo_out <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo_out <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)

human_fimo_out <- human_fimo_out[human_fimo_out$q.value < 0.05,]
macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$q.value < 0.05,]

#filter fimo out
idx <- which(human_fimo_out$start <= 31 & human_fimo_out$stop >= 31)
human_fimo_out <- human_fimo_out[idx,]

idx <- which(macaque_fimo_out$start <= 31 & macaque_fimo_out$stop >= 31)
macaque_fimo_out <- macaque_fimo_out[idx,]

#try IP
IP_DF_list <- readRDS(file = './res/step_72_fig_221101/DF_list/IP_DF_list.rds')
IP_DF_list <- IP_DF_list[IP_DF_list$group == 'human_specific',"human_SNC"]

#subset fimo out to only include DF human SNC
subset_human_fimo_out <- human_fimo_out[human_fimo_out$sequence_name %in% IP_DF_list,]
subset_macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$sequence_name %in% IP_DF_list,]

#subset fiomo out to only include the highest score
subset_human_fimo_out <- subset_human_fimo_out[order(subset_human_fimo_out$score,decreasing = TRUE),]
temp <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '_')
temp <- which(!duplicated(temp))
subset_human_fimo_out <- subset_human_fimo_out[temp,]

subset_macaque_fimo_out <- subset_macaque_fimo_out[order(subset_macaque_fimo_out$score,decreasing = TRUE),]
temp <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '_')
temp <- which(!duplicated(temp))
subset_macaque_fimo_out <- subset_macaque_fimo_out[temp,]

#create full set
temp_human <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '#')
temp_macaque <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '#')

rownames(subset_human_fimo_out) <- temp_human
rownames(subset_macaque_fimo_out) <- temp_macaque

full_set <- unique(c(temp_macaque,temp_human))
full_set <- data.frame(set = full_set,human_score = NA,macaque_score = NA)

full_set$human_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
  if(x %in% rownames(subset_human_fimo_out)){
    return(subset_human_fimo_out[x,"score"])
  }else{
    temp <- strsplit(x = x,split = '#')
    temp <- temp[[1]][1]
    if(temp %in% human_fimo_out$motif_id){
      temp <- human_fimo_out[human_fimo_out$motif_id == temp,"score"]
      return(min(temp))
    }else{
      return(NA)
    }
  }
}))

full_set$macaque_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
  if(x %in% rownames(subset_macaque_fimo_out)){
    return(subset_macaque_fimo_out[x,"score"])
  }else{
    temp <- strsplit(x = x,split = '#')
    temp <- temp[[1]][1]
    if(temp %in% macaque_fimo_out$motif_id){
      temp <- macaque_fimo_out[macaque_fimo_out$motif_id == temp,"score"]
      return(min(temp))
    }else{
      return(NA)
    }
  }
}))

full_set <- na.omit(full_set)
full_set$difference <- full_set$human_score - full_set$macaque_score
summary(full_set$difference)

#seems meaning less

# for loop access the motif affinity difference ---------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get fimo out
human_fimo_out <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo_out <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)

#back up fimo out
human_fimo_out_backup <- human_fimo_out
macaque_fimo_out_backup <- macaque_fimo_out
human_fimo_out_backup$idx <- paste(human_fimo_out_backup$motif_id,human_fimo_out_backup$sequence_name,human_fimo_out_backup$strand,sep = '#')
macaque_fimo_out_backup$idx <- paste(macaque_fimo_out_backup$motif_id,macaque_fimo_out_backup$sequence_name,macaque_fimo_out_backup$strand,sep = '#')

#filter fimo out
human_fimo_out <- human_fimo_out[human_fimo_out$q.value < 0.05,]
macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$q.value < 0.05,]

idx <- which(human_fimo_out$start <= 31 & human_fimo_out$stop >= 31)
human_fimo_out <- human_fimo_out[idx,]

idx <- which(macaque_fimo_out$start <= 31 & macaque_fimo_out$stop >= 31)
macaque_fimo_out <- macaque_fimo_out[idx,]

#get cell type list 
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

#loop
affinity_difference <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  
  #get cell type
  cell_type <- x
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load DF list
  DF_list <- list.files(path = './res/step_72_fig_221101/DF_list/')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_72_fig_221101/DF_list',DF_list,sep = '/')
  DF_list <- readRDS(file = DF_list)
  
  #human specific
  i <- 'human_specific'
  SNC_list <- DF_list[DF_list$group == i,"human_SNC"]
  subset_human_fimo_out <- human_fimo_out[human_fimo_out$sequence_name %in% SNC_list,]
  subset_macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$sequence_name %in% SNC_list,]
  
  subset_human_fimo_out <- subset_human_fimo_out[order(subset_human_fimo_out$score,decreasing = TRUE),]
  temp <- subset_human_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_human_fimo_out <- subset_human_fimo_out[temp,]
  
  subset_macaque_fimo_out <- subset_macaque_fimo_out[order(subset_macaque_fimo_out$score,decreasing = TRUE),]
  temp <- subset_macaque_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_macaque_fimo_out <- subset_macaque_fimo_out[temp,]
  
  temp_human <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '#')
  temp_macaque <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '#')
  
  rownames(subset_human_fimo_out) <- temp_human
  rownames(subset_macaque_fimo_out) <- temp_macaque
  
  full_set <- unique(c(temp_macaque,temp_human))
  full_set <- data.frame(set = full_set,human_score = NA,macaque_score = NA)
  
  full_set$human_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% human_fimo_out_backup$idx){
      temp <- human_fimo_out_backup[human_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% human_fimo_out_backup$motif_id){
        temp <- human_fimo_out_backup[human_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set$macaque_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% macaque_fimo_out_backup$idx){
      temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% macaque_fimo_out_backup$motif_id){
        temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set <- na.omit(full_set)
  full_set$difference <- full_set$human_score - full_set$macaque_score
  full_set$cell_type <- cell_type
  full_set$group <- i
  
  human_specific <- full_set
  
  #macaque specific
  i <- 'macaque_specific'
  SNC_list <- DF_list[DF_list$group == i,"human_SNC"]
  subset_human_fimo_out <- human_fimo_out[human_fimo_out$sequence_name %in% SNC_list,]
  subset_macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$sequence_name %in% SNC_list,]
  
  subset_human_fimo_out <- subset_human_fimo_out[order(subset_human_fimo_out$score,decreasing = TRUE),]
  temp <- subset_human_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_human_fimo_out <- subset_human_fimo_out[temp,]
  
  subset_macaque_fimo_out <- subset_macaque_fimo_out[order(subset_macaque_fimo_out$score,decreasing = TRUE),]
  temp <- subset_macaque_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_macaque_fimo_out <- subset_macaque_fimo_out[temp,]
  
  temp_human <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '#')
  temp_macaque <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '#')
  
  rownames(subset_human_fimo_out) <- temp_human
  rownames(subset_macaque_fimo_out) <- temp_macaque
  
  full_set <- unique(c(temp_macaque,temp_human))
  full_set <- data.frame(set = full_set,human_score = NA,macaque_score = NA)
  
  full_set$human_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% human_fimo_out_backup$idx){
      temp <- human_fimo_out_backup[human_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% human_fimo_out_backup$motif_id){
        temp <- human_fimo_out_backup[human_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set$macaque_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% macaque_fimo_out_backup$idx){
      temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% macaque_fimo_out_backup$motif_id){
        temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set <- na.omit(full_set)
  full_set$difference <- full_set$human_score - full_set$macaque_score
  full_set$cell_type <- cell_type
  full_set$group <- i
  
  macaque_specific <- full_set
  
  #species conserved
  i <- 'species_conserved'
  SNC_list <- DF_list[DF_list$group == i,"human_SNC"]
  subset_human_fimo_out <- human_fimo_out[human_fimo_out$sequence_name %in% SNC_list,]
  subset_macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$sequence_name %in% SNC_list,]
  
  subset_human_fimo_out <- subset_human_fimo_out[order(subset_human_fimo_out$score,decreasing = TRUE),]
  temp <- subset_human_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_human_fimo_out <- subset_human_fimo_out[temp,]
  
  subset_macaque_fimo_out <- subset_macaque_fimo_out[order(subset_macaque_fimo_out$score,decreasing = TRUE),]
  temp <- subset_macaque_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_macaque_fimo_out <- subset_macaque_fimo_out[temp,]
  
  temp_human <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '#')
  temp_macaque <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '#')
  
  rownames(subset_human_fimo_out) <- temp_human
  rownames(subset_macaque_fimo_out) <- temp_macaque
  
  full_set <- unique(c(temp_macaque,temp_human))
  full_set <- data.frame(set = full_set,human_score = NA,macaque_score = NA)
  
  full_set$human_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% human_fimo_out_backup$idx){
      temp <- human_fimo_out_backup[human_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% human_fimo_out_backup$motif_id){
        temp <- human_fimo_out_backup[human_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set$macaque_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% macaque_fimo_out_backup$idx){
      temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% macaque_fimo_out_backup$motif_id){
        temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set <- na.omit(full_set)
  full_set$difference <- full_set$human_score - full_set$macaque_score
  full_set$cell_type <- cell_type
  full_set$group <- i
  
  species_conserved <- full_set
  
  #return
  temp <- rbind(human_specific,macaque_specific,species_conserved)
  return(temp)
  
}))

my_send_sms('calculate done!')

#plot
ggplot(data = affinity_difference,aes(x = group,y = difference,fill = group)) + 
  geom_boxplot() + facet_wrap(~ cell_type)


# for loop access the motif affinity difference after new DF list ---------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get fimo out
human_fimo_out <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo_out <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)

#back up fimo out
human_fimo_out_backup <- human_fimo_out
macaque_fimo_out_backup <- macaque_fimo_out
human_fimo_out_backup$idx <- paste(human_fimo_out_backup$motif_id,human_fimo_out_backup$sequence_name,human_fimo_out_backup$strand,sep = '#')
macaque_fimo_out_backup$idx <- paste(macaque_fimo_out_backup$motif_id,macaque_fimo_out_backup$sequence_name,macaque_fimo_out_backup$strand,sep = '#')

#filter fimo out
human_fimo_out <- human_fimo_out[human_fimo_out$q.value < 0.05,]
macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$q.value < 0.05,]

idx <- which(human_fimo_out$start <= 31 & human_fimo_out$stop >= 31)
human_fimo_out <- human_fimo_out[idx,]

idx <- which(macaque_fimo_out$start <= 31 & macaque_fimo_out$stop >= 31)
macaque_fimo_out <- macaque_fimo_out[idx,]

#get cell type list 
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

#loop
affinity_difference <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  
  #get cell type
  cell_type <- x
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load DF list
  DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
  DF_list <- readRDS(file = DF_list)
  
  #human specific
  i <- 'human_specific'
  SNC_list <- DF_list[DF_list$group == i,"human_SNC"]
  subset_human_fimo_out <- human_fimo_out[human_fimo_out$sequence_name %in% SNC_list,]
  subset_macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$sequence_name %in% SNC_list,]
  
  subset_human_fimo_out <- subset_human_fimo_out[order(subset_human_fimo_out$score,decreasing = TRUE),]
  temp <- subset_human_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_human_fimo_out <- subset_human_fimo_out[temp,]
  
  subset_macaque_fimo_out <- subset_macaque_fimo_out[order(subset_macaque_fimo_out$score,decreasing = TRUE),]
  temp <- subset_macaque_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_macaque_fimo_out <- subset_macaque_fimo_out[temp,]
  
  temp_human <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '#')
  temp_macaque <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '#')
  
  rownames(subset_human_fimo_out) <- temp_human
  rownames(subset_macaque_fimo_out) <- temp_macaque
  
  full_set <- unique(c(temp_macaque,temp_human))
  full_set <- data.frame(set = full_set,human_score = NA,macaque_score = NA)
  
  full_set$human_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% human_fimo_out_backup$idx){
      temp <- human_fimo_out_backup[human_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% human_fimo_out_backup$motif_id){
        temp <- human_fimo_out_backup[human_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set$macaque_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% macaque_fimo_out_backup$idx){
      temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% macaque_fimo_out_backup$motif_id){
        temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set <- na.omit(full_set)
  full_set$difference <- full_set$human_score - full_set$macaque_score
  full_set$cell_type <- cell_type
  full_set$group <- i
  
  human_specific <- full_set
  
  #macaque specific
  i <- 'macaque_specific'
  SNC_list <- DF_list[DF_list$group == i,"human_SNC"]
  subset_human_fimo_out <- human_fimo_out[human_fimo_out$sequence_name %in% SNC_list,]
  subset_macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$sequence_name %in% SNC_list,]
  
  subset_human_fimo_out <- subset_human_fimo_out[order(subset_human_fimo_out$score,decreasing = TRUE),]
  temp <- subset_human_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_human_fimo_out <- subset_human_fimo_out[temp,]
  
  subset_macaque_fimo_out <- subset_macaque_fimo_out[order(subset_macaque_fimo_out$score,decreasing = TRUE),]
  temp <- subset_macaque_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_macaque_fimo_out <- subset_macaque_fimo_out[temp,]
  
  temp_human <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '#')
  temp_macaque <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '#')
  
  rownames(subset_human_fimo_out) <- temp_human
  rownames(subset_macaque_fimo_out) <- temp_macaque
  
  full_set <- unique(c(temp_macaque,temp_human))
  full_set <- data.frame(set = full_set,human_score = NA,macaque_score = NA)
  
  full_set$human_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% human_fimo_out_backup$idx){
      temp <- human_fimo_out_backup[human_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% human_fimo_out_backup$motif_id){
        temp <- human_fimo_out_backup[human_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set$macaque_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% macaque_fimo_out_backup$idx){
      temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% macaque_fimo_out_backup$motif_id){
        temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set <- na.omit(full_set)
  full_set$difference <- full_set$human_score - full_set$macaque_score
  full_set$cell_type <- cell_type
  full_set$group <- i
  
  macaque_specific <- full_set
  
  #species conserved
  i <- 'species_conserved'
  SNC_list <- DF_list[DF_list$group == i,"human_SNC"]
  subset_human_fimo_out <- human_fimo_out[human_fimo_out$sequence_name %in% SNC_list,]
  subset_macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$sequence_name %in% SNC_list,]
  
  subset_human_fimo_out <- subset_human_fimo_out[order(subset_human_fimo_out$score,decreasing = TRUE),]
  temp <- subset_human_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_human_fimo_out <- subset_human_fimo_out[temp,]
  
  subset_macaque_fimo_out <- subset_macaque_fimo_out[order(subset_macaque_fimo_out$score,decreasing = TRUE),]
  temp <- subset_macaque_fimo_out$sequence_name
  temp <- which(!duplicated(temp))
  subset_macaque_fimo_out <- subset_macaque_fimo_out[temp,]
  
  temp_human <- paste(subset_human_fimo_out$motif_id,subset_human_fimo_out$sequence_name,subset_human_fimo_out$strand,sep = '#')
  temp_macaque <- paste(subset_macaque_fimo_out$motif_id,subset_macaque_fimo_out$sequence_name,subset_macaque_fimo_out$strand,sep = '#')
  
  rownames(subset_human_fimo_out) <- temp_human
  rownames(subset_macaque_fimo_out) <- temp_macaque
  
  full_set <- unique(c(temp_macaque,temp_human))
  full_set <- data.frame(set = full_set,human_score = NA,macaque_score = NA)
  
  full_set$human_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% human_fimo_out_backup$idx){
      temp <- human_fimo_out_backup[human_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% human_fimo_out_backup$motif_id){
        temp <- human_fimo_out_backup[human_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set$macaque_score <- unlist(base::lapply(X = full_set$set,FUN = function(x){
    if(x %in% macaque_fimo_out_backup$idx){
      temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$idx == x,"score"]
      return(max(temp))
    }else{
      temp <- strsplit(x = x,split = '#')
      temp <- temp[[1]][1]
      if(temp %in% macaque_fimo_out_backup$motif_id){
        temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$motif_id == temp,"score"]
        return(min(temp))
      }else{
        return(NA)
      }
    }
  }))
  
  full_set <- na.omit(full_set)
  full_set$difference <- full_set$human_score - full_set$macaque_score
  full_set$cell_type <- cell_type
  full_set$group <- i
  
  species_conserved <- full_set
  
  #return
  temp <- rbind(human_specific,macaque_specific,species_conserved)
  return(temp)
  
}))

my_send_sms('calculate done!')

#save data
saveRDS(object = affinity_difference,file = './res/step_77_fig_221108/affinity_difference_data_frame.rds')

#plot
compare_group <- list(c('human_specific','macaque_specific'),
                      c('macaque_specific','species_conserved'),
                      c('human_specific','species_conserved'))
affinity_difference$cell_type <- factor(affinity_difference$cell_type,levels = cell_type_list)
affinity_difference$group <- factor(affinity_difference$group,levels = c('human_specific','macaque_specific','species_conserved'))

pdf(file = './res/step_77_fig_221108/affinity_difference_boxplot.pdf',width = 15,height = 8)
ggplot(data = affinity_difference,aes(x = group,y = difference,fill = group)) + 
  geom_boxplot(outlier.alpha = 0) + facet_wrap(~ cell_type,ncol = 5) + 
  theme_bw() + theme(aspect.ratio = 1) + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  stat_compare_means(comparisons = compare_group,paired = FALSE,method = 'wilcox.test',vjust = 1.5) + 
  scale_fill_manual(values = as.character(color_param$species))
dev.off()