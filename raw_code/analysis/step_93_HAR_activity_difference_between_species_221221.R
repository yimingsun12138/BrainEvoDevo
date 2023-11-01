#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: HAR activity difference between species                         ##
## Data: 2022.12.21                                                                ##
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
library(BSgenome.Mmusculus.UCSC.mm10)
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(universalmotif)
library(topGO)
library(future.apply)
library(transPlotR)
library(aplot)
library(Biostrings)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# HAR sequence difference between species ---------------------------------
#load data
HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')

#get human seq
human_seq <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,names = HAR_list$human)
names(human_seq) <- names(HAR_list$human)

macaque_seq <- BSgenome::getSeq(x = BSgenome.Mmulatta.UCSC.rheMac10,names = HAR_list$macaque)
names(macaque_seq) <- HAR_list$macaque[names(macaque_seq)]$ori_peak

mouse_seq <- BSgenome::getSeq(x = BSgenome.Mmusculus.UCSC.mm10,names = HAR_list$mouse)
names(mouse_seq) <- HAR_list$mouse[names(mouse_seq)]$ori_peak

#create score matrix
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

#find the correct order
#macaque
macaque_seq <- base::lapply(X = names(human_seq),FUN = function(x){
  #get seq
  temp_human <- human_seq[[x]]
  temp_macaque <- macaque_seq[[x]]
  
  #normal alignment
  normal_score <- pairwiseAlignment(pattern = temp_human,subject = temp_macaque,type="global",
                                    substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #complementary alignment
  comp_score <- pairwiseAlignment(pattern = temp_human,subject = Biostrings::complement(temp_macaque),type="global",
                                  substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #reverse alignment
  rev_score <- pairwiseAlignment(pattern = temp_human,subject = Biostrings::reverse(temp_macaque),type="global",
                                 substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #reverse complementary alignment
  rev_comp_score <- pairwiseAlignment(pattern = temp_human,subject = Biostrings::reverseComplement(temp_macaque),type="global",
                                      substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #compare
  if(normal_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(temp_macaque)
  }else if(comp_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(Biostrings::complement(temp_macaque))
  }else if(rev_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(Biostrings::reverse(temp_macaque))
  }else if(rev_comp_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(Biostrings::reverseComplement(temp_macaque))
  }else{
    stop('something wrong!')
  }
})
macaque_seq <- DNAStringSet(macaque_seq)
names(macaque_seq) <- names(human_seq)


macaque_score <- unlist(base::lapply(X = names(human_seq),FUN = function(x){
  #get seq
  temp_human <- human_seq[[x]]
  temp_macaque <- macaque_seq[[x]]
  
  #score
  score_num <- pairwiseAlignment(pattern = temp_human,subject = temp_macaque,type="global",
                                 substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #return
  return(score_num)
}))

#mouse
mouse_seq <- base::lapply(X = names(human_seq),FUN = function(x){
  #get seq
  temp_human <- human_seq[[x]]
  temp_mouse <- mouse_seq[[x]]
  
  #normal alignment
  normal_score <- pairwiseAlignment(pattern = temp_human,subject = temp_mouse,type="global",
                                    substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #complementary alignment
  comp_score <- pairwiseAlignment(pattern = temp_human,subject = Biostrings::complement(temp_mouse),type="global",
                                  substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #reverse alignment
  rev_score <- pairwiseAlignment(pattern = temp_human,subject = Biostrings::reverse(temp_mouse),type="global",
                                 substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #reverse complementary alignment
  rev_comp_score <- pairwiseAlignment(pattern = temp_human,subject = Biostrings::reverseComplement(temp_mouse),type="global",
                                      substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #compare
  if(normal_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(temp_mouse)
  }else if(comp_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(Biostrings::complement(temp_mouse))
  }else if(rev_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(Biostrings::reverse(temp_mouse))
  }else if(rev_comp_score == max(normal_score,comp_score,rev_score,rev_comp_score)){
    return(Biostrings::reverseComplement(temp_mouse))
  }else{
    stop('something wrong!')
  }
})
mouse_seq <- DNAStringSet(mouse_seq)
names(mouse_seq) <- names(human_seq)


mouse_score <- unlist(base::lapply(X = names(human_seq),FUN = function(x){
  #get seq
  temp_human <- human_seq[[x]]
  temp_mouse <- mouse_seq[[x]]
  
  #score
  score_num <- pairwiseAlignment(pattern = temp_human,subject = temp_mouse,type="global",
                                 substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #return
  return(score_num)
}))

#macaque_vs_mouse_score
macaque_vs_mouse_score <- unlist(base::lapply(X = names(macaque_seq),FUN = function(x){
  #get seq
  temp_macaque <- macaque_seq[[x]]
  temp_mouse <- mouse_seq[[x]]
  
  #score
  score_num <- pairwiseAlignment(pattern = temp_macaque,subject = temp_mouse,type="global",
                                 substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  
  #return
  return(score_num)
}))

#alignment score
alignment_score <- data.frame(HAR = names(human_seq),human_vs_macaque = macaque_score,human_vs_mouse = mouse_score,macaque_vs_mouse = macaque_vs_mouse_score)
rownames(alignment_score) <- alignment_score$HAR

#save data
HAR_seq <- SimpleList(human_seq = human_seq,
                      macaque_seq = macaque_seq,
                      mouse_seq = mouse_seq,
                      alignment_score = alignment_score)
saveRDS(object = HAR_seq,file = './res/step_93_fig_221216/HAR_seq.rds')

#compare human vs macaque and macaque vs mouse
hist(alignment_score$human_vs_macaque - alignment_score$macaque_vs_mouse)

# # HAR group number --------------------------------------------------------
# #load data
# HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')
# 
# Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
# macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
# mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')
# 
# color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')
# 
# #cell type list
# cell_type_list <- names(color_param$celltype)
# cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% mouse_multiome_ArchR$Gex_macaque_cell_type]
# cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]
# 
# #lapply count HAR group number
# count_matrix <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
#   #cell type
#   cell_type <- x
#   cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
#   
#   #load HAR_list
#   temp_HAR <- list.files(path = './res/step_93_fig_221212/DF_three_species')
#   temp_HAR <- temp_HAR[grep(pattern = cell_type_dot,x = temp_HAR,fixed = TRUE)]
#   temp_HAR <- paste('./res/step_93_fig_221212/DF_three_species',temp_HAR,sep = '/')
#   temp_HAR <- readRDS(file = temp_HAR)
#   
#   #count
#   group_list <- names(temp_HAR)
#   temp <- unlist(base::lapply(X = group_list,FUN = function(i){
#     temp <- temp_HAR[[i]]
#     return(length(temp))
#   }))
#   temp <- data.frame(group_list = group_list,number = temp,cell_type = cell_type)
#   return(temp)
# }))
# count_matrix$group_list <- factor(count_matrix$group_list,levels = c('human_specific','primate_conserved','human_loss','species_conserved'))

# HAR regulated genes( human specific ) -----------------------------------------------------

## calculate gene list affected by human specific HAR ----------------------

#load data
HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% mouse_multiome_ArchR$Gex_macaque_cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#lapply
gene_list <- base::lapply(X = cell_type_list,FUN = function(x){
  #cell type
  cell_type <- x
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load HAR list
  temp_HAR <- list.files(path = './res/step_93_fig_221212/DF_three_species')
  temp_HAR <- temp_HAR[grep(pattern = cell_type_dot,x = temp_HAR,fixed = TRUE)]
  temp_HAR <- paste('./res/step_93_fig_221212/DF_three_species',temp_HAR,sep = '/')
  temp_HAR <- readRDS(file = temp_HAR)
  
  temp_HAR <- temp_HAR$human_specific
  temp_HAR <- HAR_list$human[temp_HAR]
  
  #get peakset
  peakset <- p2g@metadata$peakSet
  names(peakset) <- paste(peakset@seqnames,as.character(peakset@ranges),sep = '-')
  
  #count overlap
  temp <- countOverlaps(query = peakset,subject = temp_HAR)
  temp <- names(temp)[temp > 0]
  peak_idx <- which(names(peakset) %in% temp)
  
  #idx
  idx <- which(p2g$idxATAC %in% peak_idx)
  gene_idx <- p2g$idxRNA[idx]
  gene_name <- p2g@metadata$geneSet$name[gene_idx]
  gene_name <- unique(gene_name)
  
  #return
  return(gene_name)
})
names(gene_list) <- cell_type_list

#save gene list
saveRDS(object = gene_list,file = './res/step_93_fig_221216/gene_list_affected_by_human_specific_HAR.rds')

## find cases --------------------------------------------------------------
#plot function
plot_track <- function(gene_list = gene_list,species = 'human'){
  #load data
  HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')
  macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
  Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
  mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')
  color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')
  human_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
  human_anno <- rtracklayer::as.data.frame(x = human_anno)
  macaque_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
  macaque_anno <- rtracklayer::as.data.frame(x = macaque_anno)
  mouse_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Mus_musculus.GRCm38.102.gtf',format = 'gtf')
  mouse_anno <- rtracklayer::as.data.frame(x = mouse_anno)
  
  #set parameter
  cell_type_list <- c('RG-1','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE')
  cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)
  extend_size <- 100000
  resolution <- 1000
  ylim_track <- 0.03
  
  #load bw
  if(species == 'human'){
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type/')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }else if(species == 'macaque'){
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }else{
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }
  
  #for loop
  for (i in gene_list) {
    #set gene name
    gene_name <- i
    
    #get p2g
    p2g <- getPeak2GeneLinks(
      ArchRProj = Greenleaf_ATAC_ArchR,
      corCutOff = 0.45,
      resolution = 1,
      returnLoops = FALSE
    )
    
    #get HAR list that may affect the gene express
    temp <- which(p2g@metadata$geneSet$name == gene_name)
    temp <- which(p2g$idxRNA %in% temp)
    temp <- p2g$idxATAC[temp]
    temp <- p2g@metadata$peakSet[temp]
    
    human_HAR <- HAR_list$human
    macaque_HAR <- HAR_list$macaque
    mouse_HAR <- HAR_list$mouse
    names(macaque_HAR) <- macaque_HAR$ori_peak
    names(mouse_HAR) <- mouse_HAR$ori_peak
    
    temp <- countOverlaps(query = HAR_list$human,subject = temp)
    temp_HAR <- names(temp)[temp > 0]
    if(species == 'human'){
      temp_HAR <- human_HAR[temp_HAR]
    }else if(species == 'macaque'){
      temp_HAR <- macaque_HAR[temp_HAR]
    }else{
      temp_HAR <- mouse_HAR[temp_HAR]
    }
    
    #get target site
    if(species == 'human'){
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = Greenleaf_ATAC_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }else if(species == 'macaque'){
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = macaque_multiome_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }else{
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = mouse_multiome_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }
    
    if(species == 'mouse'){
      mouse_gene_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCh38_to_GRCm39.csv')
      mouse_gene_anno <- mouse_gene_anno[which(mouse_gene_anno[,2] == gene_name),,drop = FALSE]
      gene_name <- unique(mouse_gene_anno[,4])[1]
    }else if(species == 'macaque'){
      macaque_gene_anno <- readRDS(file = '/home/sunym/temp/macaque_anno.rds')
      macaque_gene_anno <- macaque_gene_anno[which(macaque_gene_anno[,2] == gene_name),,drop = FALSE]
      gene_name <- unique(macaque_gene_anno[,4])[1]
    }
    TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]
    if(length(TSS_site) > 1){
      stop('TSS number error!')
    }
    chrom <- as.character(TSS_site@seqnames)
    start_site <- TSS_site@ranges@start - extend_size
    end_site <- TSS_site@ranges@start + extend_size
    
    #modify extension
    temp <- rtracklayer::as.data.frame(x = temp_HAR)
    start_range <- min(temp$start,temp$end)
    end_range <- max(temp$start,temp$end)
    
    if(start_site >= start_range){
      start_site <- start_range - 10000
    }
    if(end_site <= end_range){
      end_site <- end_range + 10000
    }
    
    #plot gene annotation
    if(species == 'human'){
      p_gene_anno <- trancriptVis(gtfFile = human_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }else if(species == 'macaque'){
      p_gene_anno <- trancriptVis(gtfFile = macaque_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }else{
      p_gene_anno <- trancriptVis(gtfFile = mouse_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }
    
    #plot peak
    if(species == 'human'){
      temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
      temp <- temp_HAR
      if(length(temp) == 0){
        temp <- human_HAR
      }
      rtracklayer::export(object = temp,con = '/home/sunym/temp/temp_HAR.bed',format = 'bed')
      rtracklayer::export(object = human_HAR,con = '/home/sunym/temp/HAR_list.bed',format = 'bed')
    }else if(species == 'macaque'){
      temp <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
      temp <- temp_HAR
      if(length(temp) == 0){
        temp <- macaque_HAR
      }
      rtracklayer::export(object = temp,con = '/home/sunym/temp/temp_HAR.bed',format = 'bed')
      rtracklayer::export(object = macaque_HAR,con = '/home/sunym/temp/HAR_list.bed',format = 'bed')
    }else{
      temp <- getPeakSet(ArchRProj = mouse_multiome_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
      temp <- temp_HAR
      if(length(temp) == 0){
        temp <- mouse_HAR
      }
      rtracklayer::export(object = temp,con = '/home/sunym/temp/temp_HAR.bed',format = 'bed')
      rtracklayer::export(object = mouse_HAR,con = '/home/sunym/temp/HAR_list.bed',format = 'bed')
    }
    
    p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed','/home/sunym/temp/HAR_list.bed','/home/sunym/temp/temp_HAR.bed'),
                          chr = chrom,region.min = start_site - 10000,region.max = end_site + 10000,
                          show.legend = FALSE,fill = ggsci::pal_d3()(3),track.width = 0.3)
    
    #plot track
    if(species == 'human'){
      temp_anno <- human_anno
    }else if(species == 'macaque'){
      temp_anno <- macaque_anno
    }else{
      temp_anno <- mouse_anno
    }
    
    p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = temp_anno,
                             chr = chrom,region.min = start_site,region.max = end_site,
                             sample.order = cell_type_list,space.y = 0,
                             y.max = ylim_track,color = color_param$celltype[cell_type_list],
                             theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)
    
    #plot link
    if(species == 'human'){
      p2g <- getPeak2GeneLinks(
        ArchRProj = Greenleaf_ATAC_ArchR,
        corCutOff = 0.45,
        resolution = resolution,
        returnLoops = FALSE
      )
    }else if(species == 'macaque'){
      p2g <- getPeak2GeneLinks(
        ArchRProj = macaque_multiome_ArchR,
        corCutOff = 0.45,
        resolution = resolution,
        returnLoops = FALSE
      )
    }else{
      p2g <- getPeak2GeneLinks(
        ArchRProj = mouse_multiome_ArchR,
        corCutOff = 0.45,
        resolution = resolution,
        returnLoops = FALSE
      )
    }
    
    temp <- which(p2g@metadata$geneSet$name == gene_name)
    temp <- which(p2g$idxRNA %in% temp)
    if(length(temp) == 0){
      links_meta_data <- data.frame(chr = chrom,
                                    start = c(start_site - 3000,start_site - 6000),
                                    end = c(end_site + 3000,end_site + 6000),
                                    cor = c(1,0),group = 'group1')
    }else{
      idx_ATAC <- p2g$idxATAC[temp]
      idx_ATAC <- p2g@metadata$peakSet[idx_ATAC]
      idx_ATAC <- rtracklayer::as.data.frame(idx_ATAC)
      idx_ATAC <- round((idx_ATAC$start + idx_ATAC$end) / 2)
      links_meta_data <- data.frame(chr = chrom,start = NA,end = NA,cor = p2g$Correlation[temp],group = 'group1')
      for (i in 1:length(idx_ATAC)) {
        if(TSS_site@ranges@start <= idx_ATAC[i]){
          links_meta_data[i,'start'] <- TSS_site@ranges@start
          links_meta_data[i,"end"] <- idx_ATAC[i]
        }else{
          links_meta_data[i,"start"] <- idx_ATAC[i]
          links_meta_data[i,'end'] <- TSS_site@ranges@start
        }
      }
      temp <- data.frame(chr = chrom,
                         start = c(start_site - 3000,start_site - 6000),
                         end = c(end_site + 3000,end_site + 6000),
                         cor = c(1,0),group = 'group1')
      links_meta_data <- rbind(links_meta_data,temp)
    }
    p_p2g <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                     facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),
                     xAixs.info = FALSE,curvature = 0.3)
    
    #save plot
    char <- paste0('./res/step_93_fig_221216/',species,'_',gene_name,'_trackplot.pdf')
    pdf(file = char,width = 8,height = 8)
    print(p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.05) %>% insert_bottom(plot = p_p2g,height = 0.1) %>% insert_bottom(plot = p_gene_anno,height = 0.05))
    dev.off()
  }
  return(NULL)
}

#gene_list
gene_list <- readRDS(file = './res/step_93_fig_221216/gene_list_affected_by_human_specific_HAR.rds')
gene_list <- unlist(gene_list)
names(gene_list) <- NULL
gene_list <- unique(gene_list)

#modify gene list
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
macaque_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCh38_to_Mmul10.csv')
macaque_anno <- macaque_anno[macaque_anno[,2] %in% gene_list,]
table(duplicated(macaque_anno))
temp <- gene_list[!(gene_list %in% macaque_anno[,2])]
table(temp %in% rownames(macaque_multiome_Seurat@assays$RNA@counts))
temp <- temp[!(temp %in% rownames(macaque_multiome_Seurat@assays$RNA@counts))]
gene_list <- gene_list[!(gene_list %in% temp)]

#modify macaque anno
macaque_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCh38_to_Mmul10.csv')
macaque_anno <- macaque_anno[macaque_anno[,2] %in% gene_list,]
table(duplicated(macaque_anno))

temp <- gene_list[!(gene_list %in% macaque_anno[,2])]
table(temp %in% rownames(macaque_multiome_Seurat@assays$RNA@counts))
temp <- data.frame(Gene.stable.ID = NA,Gene.name = temp,Gene.stable.ID.1 = NA,Gene.name.1 = temp)
macaque_anno <- rbind(macaque_anno,temp)

table(duplicated(macaque_anno[,2]))
table(gene_list %in% macaque_anno[,2])
saveRDS(object = macaque_anno,file = '/home/sunym/temp/macaque_anno.rds')

#plot
plot_track(gene_list = gene_list,species = 'macaque')
plot_track(gene_list = gene_list,species = 'mouse')
plot_track(gene_list = gene_list,species = 'human')
