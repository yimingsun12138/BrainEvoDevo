#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: HAR activity difference between species                         ##
## Data: 2022.12.12                                                                ##
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

# liftover to macaque and mouse -------------------------------------------
#load data
human_HAR <- readRDS(file = './res/step_92_fig_221209/processed_HAR_region.rds')
names(human_HAR) <- paste(human_HAR@seqnames,as.character(human_HAR@ranges),sep = '-')

#liftover to macaque
macaque_HAR <- my_rtracklayer_liftOver(ori_GRanges = human_HAR,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac10.over.chain',merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
table(macaque_HAR@seqnames)
summary((macaque_HAR@ranges@width - human_HAR[macaque_HAR$ori_peak]@ranges@width)/human_HAR[macaque_HAR$ori_peak]@ranges@width)
macaque_HAR <- macaque_HAR[abs(macaque_HAR@ranges@width - human_HAR[macaque_HAR$ori_peak]@ranges@width)/human_HAR[macaque_HAR$ori_peak]@ranges@width <= 10]

#liftover to mouse
mouse_HAR <- my_rtracklayer_liftOver(ori_GRanges = human_HAR,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToMm10.over.chain',merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
table(mouse_HAR@seqnames)
summary((mouse_HAR@ranges@width - human_HAR[mouse_HAR$ori_peak]@ranges@width)/human_HAR[mouse_HAR$ori_peak]@ranges@width)
#mouse_HAR <- mouse_HAR[abs(mouse_HAR@ranges@width - human_HAR[mouse_HAR$ori_peak]@ranges@width)/human_HAR[mouse_HAR$ori_peak]@ranges@width <= 10]

#HAR_list
HAR_list <- dplyr::intersect(x = macaque_HAR$ori_peak,y = mouse_HAR$ori_peak)

#generate final har list
human_HAR <- readRDS(file = './res/step_92_fig_221209/processed_HAR_region.rds')
names(human_HAR) <- paste(human_HAR@seqnames,as.character(human_HAR@ranges),sep = '-')
human_HAR <- human_HAR[HAR_list]

macaque_HAR <- my_rtracklayer_liftOver(ori_GRanges = human_HAR,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac10.over.chain',merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
table(macaque_HAR$ori_peak == HAR_list)
names(macaque_HAR) <- paste(macaque_HAR@seqnames,as.character(macaque_HAR@ranges),sep = '-')

mouse_HAR <- my_rtracklayer_liftOver(ori_GRanges = human_HAR,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToMm10.over.chain',merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
table(mouse_HAR$ori_peak == HAR_list)
names(mouse_HAR) <- paste(mouse_HAR@seqnames,as.character(mouse_HAR@ranges),sep = '-')

HAR_list <- SimpleList(human = human_HAR,macaque = macaque_HAR,mouse = mouse_HAR)
saveRDS(object = HAR_list,file = './res/step_93_fig_221212/HAR_list.rds')

# generate peak matrix ----------------------------------------------------
#load data
HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = '/home/sunym/temp/merged_Greenleaf_ATAC_ArchR_221018/')
macaque_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/mouse_multiome_ArchR_221009/')

#human peak matrix
Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = HAR_list$human,
                                   genomeAnnotation = getGenomeAnnotation(ArchRProj = Greenleaf_ATAC_ArchR),
                                   force = TRUE)
Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR,verbose = TRUE,force = TRUE)
human_peak_matrix <- getMatrixFromProject(ArchRProj = Greenleaf_ATAC_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)
gene_list <- human_peak_matrix@rowRanges
gene_list <- paste(gene_list@seqnames,as.character(gene_list@ranges),sep = '-')
human_peak_matrix <- human_peak_matrix@assays@data$PeakMatrix
rownames(human_peak_matrix) <- gene_list

saveRDS(object = human_peak_matrix,file = './res/step_93_fig_221212/human_peak_matrix.rds')

#macaque peak matrix
macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = HAR_list$macaque,
                                     genomeAnnotation = getGenomeAnnotation(ArchRProj = macaque_multiome_ArchR),
                                     force = TRUE)
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)
macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)
gene_list <- macaque_peak_matrix@rowRanges
gene_list <- paste(gene_list@seqnames,as.character(gene_list@ranges),sep = '-')
macaque_peak_matrix <- macaque_peak_matrix@assays@data$PeakMatrix
rownames(macaque_peak_matrix) <- gene_list

saveRDS(object = macaque_peak_matrix,file = './res/step_93_fig_221212/macaque_peak_matrix.rds')

#mouse peak matrix
mouse_multiome_ArchR <- addPeakSet(ArchRProj = mouse_multiome_ArchR,peakSet = HAR_list$mouse,
                                   genomeAnnotation = getGenomeAnnotation(ArchRProj = mouse_multiome_ArchR),
                                   force = TRUE)
mouse_multiome_ArchR <- addPeakMatrix(ArchRProj = mouse_multiome_ArchR,verbose = TRUE,force = TRUE)
mouse_peak_matrix <- getMatrixFromProject(ArchRProj = mouse_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)
gene_list <- mouse_peak_matrix@rowRanges
gene_list <- paste(gene_list@seqnames,as.character(gene_list@ranges),sep = '-')
mouse_peak_matrix <- mouse_peak_matrix@assays@data$PeakMatrix
rownames(mouse_peak_matrix) <- gene_list

saveRDS(object = mouse_peak_matrix,file = './res/step_93_fig_221212/mouse_peak_matrix.rds')

# normalize peak matrix ---------------------------------------------------

#load data
HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

human_peak_matrix <- readRDS(file = './res/step_93_fig_221212/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_93_fig_221212/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_93_fig_221212/mouse_peak_matrix.rds')

#norm factor
depth_factor <- c(Greenleaf_ATAC_ArchR$ReadsInPeaks,macaque_multiome_ArchR$ReadsInPeaks,mouse_multiome_ArchR$ReadsInPeaks)
depth_factor <- median(depth_factor)

length_factor <- c(HAR_list$human@ranges@width,HAR_list$macaque@ranges@width,HAR_list$mouse@ranges@width)
length_factor <- median(length_factor)

#normalize matrix depth
gene_list <- rownames(human_peak_matrix)
cell_list <- colnames(human_peak_matrix)
human_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = cell_list,FUN = function(x){
  temp <- human_peak_matrix[,x]
  temp <- temp/Greenleaf_ATAC_ArchR@cellColData[x,"ReadsInPeaks"]*depth_factor*10^4
  return(temp)
}))
rownames(human_peak_matrix) <- gene_list
colnames(human_peak_matrix) <- cell_list

gene_list <- rownames(macaque_peak_matrix)
cell_list <- colnames(macaque_peak_matrix)
macaque_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = cell_list,FUN = function(x){
  temp <- macaque_peak_matrix[,x]
  temp <- temp/macaque_multiome_ArchR@cellColData[x,"ReadsInPeaks"]*depth_factor*10^4
  return(temp)
}))
rownames(macaque_peak_matrix) <- gene_list
colnames(macaque_peak_matrix) <- cell_list

gene_list <- rownames(mouse_peak_matrix)
cell_list <- colnames(mouse_peak_matrix)
mouse_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = cell_list,FUN = function(x){
  temp <- mouse_peak_matrix[,x]
  temp <- temp/mouse_multiome_ArchR@cellColData[x,"ReadsInPeaks"]*depth_factor*10^4
  return(temp)
}))
rownames(mouse_peak_matrix) <- gene_list
colnames(mouse_peak_matrix) <- cell_list

#normalize HAR length
gene_list <- rownames(human_peak_matrix)
cell_list <- colnames(human_peak_matrix)
human_HAR <- HAR_list$human
names(human_HAR) <- paste(human_HAR@seqnames,as.character(human_HAR@ranges),sep = '-')
human_peak_matrix <- base::do.call(what = rbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- human_peak_matrix[x,]
  temp <- temp/human_HAR[x]@ranges@width*length_factor
  return(temp)
}))
rownames(human_peak_matrix) <- gene_list
colnames(human_peak_matrix) <- cell_list

gene_list <- rownames(macaque_peak_matrix)
cell_list <- colnames(macaque_peak_matrix)
macaque_HAR <- HAR_list$macaque
names(macaque_HAR) <- paste(macaque_HAR@seqnames,as.character(macaque_HAR@ranges),sep = '-')
macaque_peak_matrix <- base::do.call(what = rbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- macaque_peak_matrix[x,]
  temp <- temp/macaque_HAR[x]@ranges@width*length_factor
  return(temp)
}))
rownames(macaque_peak_matrix) <- gene_list
colnames(macaque_peak_matrix) <- cell_list

gene_list <- rownames(mouse_peak_matrix)
cell_list <- colnames(mouse_peak_matrix)
mouse_HAR <- HAR_list$mouse
names(mouse_HAR) <- paste(mouse_HAR@seqnames,as.character(mouse_HAR@ranges),sep = '-')
mouse_peak_matrix <- base::do.call(what = rbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- mouse_peak_matrix[x,]
  temp <- temp/mouse_HAR[x]@ranges@width*length_factor
  return(temp)
}))
rownames(mouse_peak_matrix) <- gene_list
colnames(mouse_peak_matrix) <- cell_list

#save data
peak_matrix <- SimpleList(human = human_peak_matrix,macaque = macaque_peak_matrix,mouse = mouse_peak_matrix)
saveRDS(object = peak_matrix,file = './res/step_93_fig_221212/normalized_peak_matrix.rds')

# DF analysis -------------------------------------------------------------
#load data
HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')
peak_matrix <- readRDS(file = './res/step_93_fig_221212/normalized_peak_matrix.rds')

human_peak_matrix <- peak_matrix$human
macaque_peak_matrix <- peak_matrix$macaque
mouse_peak_matrix <- peak_matrix$mouse

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#unify har names
human_HAR <- HAR_list$human
macaque_HAR <- HAR_list$macaque
mouse_HAR <- HAR_list$mouse

table(rownames(macaque_peak_matrix) %in% names(macaque_HAR))
rownames(macaque_peak_matrix) <- macaque_HAR[rownames(macaque_peak_matrix)]$ori_peak

table(rownames(mouse_peak_matrix) %in% names(mouse_HAR))
rownames(mouse_peak_matrix) <- mouse_HAR[rownames(mouse_peak_matrix)]$ori_peak

#cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% mouse_multiome_ArchR$Gex_macaque_cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

#generate cell type HAR accessibility difference
for (i in cell_type_list) {
  #cell type
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  lfc_cutoff <- 1
  fdr_cutoff <- 0.01
  
  #gene list and cell list
  gene_list <- names(human_HAR)
  human_cell_list <- rownames(Greenleaf_ATAC_ArchR@cellColData)[Greenleaf_ATAC_ArchR$cell_type == cell_type]
  macaque_cell_list <- rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$cell_type == cell_type]
  mouse_cell_list <- rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Gex_macaque_cell_type == cell_type]
  
  #do wilcox
  DF_human_macaque <- my_DF_wilcox_test(mat1 = human_peak_matrix[gene_list,human_cell_list],mat2 = macaque_peak_matrix[gene_list,macaque_cell_list],alternative = 'two.sided',paired = FALSE,workers = 5,future.globals.maxSize = 200*(1024^3))
  DF_human_mouse <- my_DF_wilcox_test(mat1 = human_peak_matrix[gene_list,human_cell_list],mat2 = mouse_peak_matrix[gene_list,mouse_cell_list],alternative = 'two.sided',paired = FALSE,workers = 5,future.globals.maxSize = 200*(1024^3))
  DF_macaque_mouse <- my_DF_wilcox_test(mat1 = macaque_peak_matrix[gene_list,macaque_cell_list],mat2 = mouse_peak_matrix[gene_list,mouse_cell_list],alternative = 'two.sided',paired = FALSE,workers = 5,future.globals.maxSize = 200*(1024^3))
  
  #load cell type peakset
  human_peakset <- list.files(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls')
  human_peakset <- human_peakset[grep(pattern = cell_type_dot,x = human_peakset,fixed = TRUE)]
  human_peakset <- paste('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls',human_peakset,sep = '/')
  human_peakset <- readRDS(file = human_peakset)
  
  macaque_peakset <- list.files(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls')
  macaque_peakset <- macaque_peakset[grep(pattern = cell_type_dot,x = macaque_peakset,fixed = TRUE)]
  macaque_peakset <- paste('./processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls',macaque_peakset,sep = '/')
  macaque_peakset <- readRDS(file = macaque_peakset)
  
  mouse_peakset <- list.files(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls')
  mouse_peakset <- mouse_peakset[grep(pattern = cell_type_dot,x = mouse_peakset,fixed = TRUE)]
  mouse_peakset <- paste('./processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls',mouse_peakset,sep = '/')
  mouse_peakset <- readRDS(file = mouse_peakset)
  
  #human specific
  temp_1 <- rownames(DF_human_macaque)[DF_human_macaque$log2FC > lfc_cutoff & DF_human_macaque$fdr < fdr_cutoff]
  temp_2 <- rownames(DF_human_mouse)[DF_human_mouse$log2FC > lfc_cutoff & DF_human_mouse$fdr < fdr_cutoff]
  human_specific <- dplyr::intersect(x = temp_1,y = temp_2)
  
  temp <- human_HAR[human_specific]
  temp <- countOverlaps(query = temp,subject = human_peakset)
  human_specific <- names(temp)[temp > 0]
  
  #primate conserved
  temp_1 <- rownames(DF_human_macaque)[abs(DF_human_macaque$log2FC) <= lfc_cutoff]
  temp_2 <- rownames(DF_human_mouse)[DF_human_mouse$log2FC > lfc_cutoff & DF_human_mouse$fdr < fdr_cutoff]
  temp_3 <- rownames(DF_macaque_mouse)[DF_macaque_mouse$log2FC > lfc_cutoff & DF_macaque_mouse$fdr < fdr_cutoff]
  primate_conserved <- dplyr::intersect(x = temp_1,y = temp_2)
  primate_conserved <- dplyr::intersect(x = primate_conserved,y = temp_3)
  
  temp <- human_HAR[primate_conserved]
  temp <- countOverlaps(query = temp,subject = human_peakset)
  primate_conserved <- names(temp)[temp > 0]
  
  names(macaque_HAR) <- macaque_HAR$ori_peak
  temp <- macaque_HAR[primate_conserved]
  temp <- countOverlaps(query = temp,subject = macaque_peakset)
  primate_conserved <- names(temp)[temp > 0]
  
  #human loss
  temp_1 <- rownames(DF_human_macaque)[DF_human_macaque$log2FC < -1*lfc_cutoff & DF_human_macaque$fdr < fdr_cutoff]
  temp_2 <- rownames(DF_human_mouse)[DF_human_mouse$log2FC < -1*lfc_cutoff & DF_human_mouse$fdr < fdr_cutoff]
  human_loss <- dplyr::intersect(x = temp_1,temp_2)
  
  names(macaque_HAR) <- macaque_HAR$ori_peak
  temp <- macaque_HAR[human_loss]
  temp <- countOverlaps(query = temp,subject = macaque_peakset)
  human_loss <- names(temp)[temp > 0]
  
  names(mouse_HAR) <- mouse_HAR$ori_peak
  temp <- mouse_HAR[human_loss]
  temp <- countOverlaps(query = temp,subject = mouse_peakset)
  human_loss <- names(temp)[temp > 0]
  
  #species conserved
  temp_1 <- rownames(DF_human_macaque)[abs(DF_human_macaque$log2FC) <= lfc_cutoff]
  temp_2 <- rownames(DF_human_mouse)[abs(DF_human_mouse$log2FC) <= lfc_cutoff]
  species_conserved <- dplyr::intersect(x = temp_1,y = temp_2)
  
  temp <- human_HAR[species_conserved]
  temp <- countOverlaps(query = temp,subject = human_peakset)
  species_conserved <- names(temp)[temp > 0]
  
  names(macaque_HAR) <- macaque_HAR$ori_peak
  temp <- macaque_HAR[species_conserved]
  temp <- countOverlaps(query = temp,subject = macaque_peakset)
  species_conserved <- names(temp)[temp > 0]
  
  names(mouse_HAR) <- mouse_HAR$ori_peak
  temp <- mouse_HAR[species_conserved]
  temp <- countOverlaps(query = temp,subject = mouse_peakset)
  species_conserved <- names(temp)[temp > 0]
  
  #return
  HAR_list <- SimpleList(human_specific = human_specific,
                         primate_conserved = primate_conserved,
                         human_loss = human_loss,
                         species_conserved = species_conserved)
  char <- paste0('./res/step_93_fig_221212/DF_three_species/',cell_type_dot,'_DF_list.rds')
  saveRDS(object = HAR_list,file = char)
  print(paste(cell_type,'done!',sep = ' '))
}

# enrich heatmap plot -----------------------------------------------------
#load data
HAR_list <- readRDS(file = './res/step_93_fig_221212/HAR_list.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#modify HAR list
human_HAR <- HAR_list$human
macaque_HAR <- HAR_list$macaque
mouse_HAR <- HAR_list$mouse

names(macaque_HAR) <- macaque_HAR$ori_peak
names(mouse_HAR) <- mouse_HAR$ori_peak

#cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% mouse_multiome_ArchR$Gex_macaque_cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

#for loop generate the enrich heatmap
for (i in cell_type_list) {
  #cell type
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load DF_list
  DF_list <- list.files(path = './res/step_93_fig_221212/DF_three_species')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_93_fig_221212/DF_three_species',DF_list,sep = '/')
  DF_list <- readRDS(file = DF_list)
  
  temp <- c(rep('human_specific',times = length(DF_list$human_specific)),
            rep('primate_conserved',times = length(DF_list$primate_conserved)),
            rep('human_loss',times = length(DF_list$human_loss)),
            rep('species_conserved',times = length(DF_list$species_conserved)))
  names(temp) <- c(DF_list$human_specific,DF_list$primate_conserved,DF_list$human_loss,DF_list$species_conserved)
  DF_list <- temp
  
  #human
  target_site <- human_HAR[names(DF_list)]
  target_site <- paste(target_site@seqnames,as.character(target_site@ranges),sep = '-')
  temp <- unlist(base::lapply(X = target_site,FUN = function(x){
    temp <- strsplit(x = x,split = '-')[[1]]
    mid_point <- round((as.numeric(temp[2]) + as.numeric(temp[3]))/2)
    temp <- paste0(temp[1],':',mid_point,'-',mid_point)
    return(temp)
  }))
  target_site <- as(temp,'GRanges')
  
  signal_coverage <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type')
  signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
  signal_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type',signal_coverage,sep = '/')
  signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
  mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p1 <- EnrichedHeatmap(mat = mat,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_loss','species_conserved')),use_raster = TRUE)
  
  #macaque
  target_site <- macaque_HAR[names(DF_list)]
  target_site <- paste(target_site@seqnames,as.character(target_site@ranges),sep = '-')
  temp <- unlist(base::lapply(X = target_site,FUN = function(x){
    temp <- strsplit(x = x,split = '-')[[1]]
    mid_point <- round((as.numeric(temp[2]) + as.numeric(temp[3]))/2)
    temp <- paste0(temp[1],':',mid_point,'-',mid_point)
    return(temp)
  }))
  target_site <- as(temp,'GRanges')
  
  signal_coverage <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type')
  signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
  signal_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type',signal_coverage,sep = '/')
  signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
  mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p2 <- EnrichedHeatmap(mat = mat,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_loss','species_conserved')),use_raster = TRUE)
  
  #mouse
  target_site <- mouse_HAR[names(DF_list)]
  target_site <- paste(target_site@seqnames,as.character(target_site@ranges),sep = '-')
  temp <- unlist(base::lapply(X = target_site,FUN = function(x){
    temp <- strsplit(x = x,split = '-')[[1]]
    mid_point <- round((as.numeric(temp[2]) + as.numeric(temp[3]))/2)
    temp <- paste0(temp[1],':',mid_point,'-',mid_point)
    return(temp)
  }))
  target_site <- as(temp,'GRanges')
  
  signal_coverage <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type')
  signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
  signal_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type',signal_coverage,sep = '/')
  signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
  mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p3 <- EnrichedHeatmap(mat = mat,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_loss','species_conserved')),use_raster = TRUE)
  
  #plot
  col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
  char <- paste0('./res/step_93_fig_221212/Enrich_heatmap_plot/',cell_type_dot,'_enrich_heatmap.pdf')
  pdf(file = char,width = 3.6,height = 6.5)
  print(EnrichedHeatmap(mat = p1@matrix,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_loss','species_conserved')),
                        use_raster = TRUE,col = col_fun,name = 'insertion',pos_line = FALSE) + 
          EnrichedHeatmap(mat = p2@matrix,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_loss','species_conserved')),
                          use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order) + 
          EnrichedHeatmap(mat = p3@matrix,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_loss','species_conserved')),
                          use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order))
  dev.off()
  
  #report
  char <- paste(cell_type_dot,'done!',sep = ' ')
  print(char)
}

# calculate HAR activity percentage ---------------------------------------
#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% mouse_multiome_ArchR$Gex_macaque_cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

#get percentage
Proportion_table <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  #set cell type
  cell_type <- x
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load DF_list
  DF_list <- list.files(path = './res/step_93_fig_221212/DF_three_species')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_93_fig_221212/DF_three_species',DF_list,sep = '/')
  DF_list <- readRDS(file = DF_list)
  
  #get proportion
  temp <- data.frame(group = names(DF_list),Proportion = NA,cell_type = cell_type)
  temp$Proportion <- unlist(base::lapply(X = temp$group,FUN = function(y){
    return(length(DF_list[[y]]))
  }))
  temp$Proportion <- temp$Proportion/sum(temp$Proportion)
  
  #return
  return(temp)
}))

#plot
Proportion_table$cell_type <- factor(Proportion_table$cell_type,levels = rev(cell_type_list))
Proportion_table$group <- factor(Proportion_table$group,levels = rev(names(DF_list)))
ggplot(data = Proportion_table,aes(x = cell_type,y = Proportion,fill = group)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7) + 
  theme_bw() + coord_flip() + 
  scale_fill_manual(values = pal_npg()(4)) + 
  theme(aspect.ratio = 0.5)
