#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: DAP analysis across species                                     ##
## Data: 2023.04.22                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.1/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(dplyr)
library(Seurat)
library(ggplot2)
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
library(riverplot)
library(CellChat)
library(OpenAI4R)
library(paletteer)
library(ggpattern)
library(ggrepel)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0314')

# re-name consensus peak set and peak matrix ------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_73_fig_221102/Brain_ATAC_peak.rds')
human_peak_matrix <- readRDS(file = './res/step_73_fig_221102/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_73_fig_221102/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_73_fig_221102/mouse_peak_matrix.rds')

#re-create consensus peakset
consensus_peakset <- as.character(Brain_ATAC_peakset$human)
names(consensus_peakset) <- NULL
consensus_peakset <- as(consensus_peakset,"GRanges")

names(consensus_peakset) <- as.character(consensus_peakset)
consensus_peakset@elementMetadata$human_coord <- as.character(consensus_peakset)

temp <- Brain_ATAC_peakset$macaque
table(names(temp) == temp$name)
temp <- temp[gsub(pattern = ':',replacement = '-',x = consensus_peakset$human_coord,fixed = TRUE)]
consensus_peakset@elementMetadata$macaque_coord <- as.character(temp)

temp <- Brain_ATAC_peakset$mouse
table(names(temp) == temp$name)
temp <- temp[gsub(pattern = ':',replacement = '-',x = consensus_peakset$human_coord,fixed = TRUE)]
consensus_peakset@elementMetadata$mouse_coord <- as.character(temp)

saveRDS(object = consensus_peakset,file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')
temp <- rtracklayer::as.data.frame(consensus_peakset)
write.csv(x = temp,file = './res/step_118_fig_230422/consensus_peakset.csv')

#re-name peak matrix
temp <- rownames(human_peak_matrix)
temp <- base::lapply(X = temp,FUN = function(x){
  peak_name <- strsplit(x = x,split = '-')[[1]]
  peak_name <- paste0(peak_name[1],":",peak_name[2],"-",peak_name[3])
  return(peak_name)
})
rownames(human_peak_matrix) <- unlist(temp)
saveRDS(object = human_peak_matrix,file = './res/step_118_fig_230422/human_peak_matrix.rds')

temp <- rownames(macaque_peak_matrix)
temp <- base::lapply(X = temp,FUN = function(x){
  peak_name <- strsplit(x = x,split = '-')[[1]]
  peak_name <- paste0(peak_name[1],":",peak_name[2],"-",peak_name[3])
  return(peak_name)
})
rownames(macaque_peak_matrix) <- unlist(temp)
saveRDS(object = macaque_peak_matrix,file = './res/step_118_fig_230422/macaque_peak_matrix.rds')

temp <- rownames(mouse_peak_matrix)
temp <- base::lapply(X = temp,FUN = function(x){
  peak_name <- strsplit(x = x,split = '-')[[1]]
  peak_name <- paste0(peak_name[1],":",peak_name[2],"-",peak_name[3])
  return(peak_name)
})
rownames(mouse_peak_matrix) <- unlist(temp)
saveRDS(object = mouse_peak_matrix,file = './res/step_118_fig_230422/mouse_peak_matrix.rds')

# wilcox ------------------------------------------------------------------
#load data
consensus_peakset <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')
human_peak_matrix <- readRDS(file = './res/step_118_fig_230422/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_118_fig_230422/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_118_fig_230422/mouse_peak_matrix.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

cell_type_list <- base::Reduce(f = base::intersect,x = list(human = as.character(Greenleaf_ATAC_ArchR$cell_type),
                                                            macaque = as.character(macaque_multiome_ArchR$cell_type),
                                                            mouse = as.character(mouse_multiome_ArchR$Gex_macaque_cell_type)))

#for loop
for (cell_type in cell_type_list) {
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load cell type peaks
  human_peaks <- list.files(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls')
  human_peaks <- human_peaks[grep(pattern = cell_type_dot,x = human_peaks,fixed = TRUE)]
  human_peaks <- paste('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls',human_peaks,sep = '/')
  human_peaks <- readRDS(human_peaks)
  
  macaque_peaks <- list.files(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls')
  macaque_peaks <- macaque_peaks[grep(pattern = cell_type_dot,x = macaque_peaks,fixed = TRUE)]
  macaque_peaks <- paste('./processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls',macaque_peaks,sep = '/')
  macaque_peaks <- readRDS(macaque_peaks)
  
  mouse_peaks <- list.files(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls')
  mouse_peaks <- mouse_peaks[grep(pattern = cell_type_dot,x = mouse_peaks,fixed = TRUE)]
  mouse_peaks <- paste('./processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls',mouse_peaks,sep = '/')
  mouse_peaks <- readRDS(mouse_peaks)
  
  #get consensus peaks that overlap with original cell type peaks in at least one species
  human_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$human_coord,'GRanges'),subject = human_peaks))
  macaque_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$macaque_coord,'GRanges'),subject = macaque_peaks))
  mouse_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$human_coord,'GRanges'),subject = mouse_peaks))
  peak_list <- consensus_peakset$human_coord[which((human_overlap + macaque_overlap + mouse_overlap) > 0)]
  
  #get subsetted peak matrix
  subset_human_peak_matrix <- human_peak_matrix[,rownames(Greenleaf_ATAC_ArchR@cellColData)[Greenleaf_ATAC_ArchR$cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$human_coord
  rownames(subset_human_peak_matrix) <- temp[rownames(subset_human_peak_matrix)]$human_coord
  
  subset_macaque_peak_matrix <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$macaque_coord
  rownames(subset_macaque_peak_matrix) <- temp[rownames(subset_macaque_peak_matrix)]$human_coord
  
  subset_mouse_peak_matrix <- mouse_peak_matrix[,rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Gex_macaque_cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$mouse_coord
  rownames(subset_mouse_peak_matrix) <- temp[rownames(subset_mouse_peak_matrix)]$human_coord
  
  subset_human_peak_matrix <- subset_human_peak_matrix[peak_list,]
  subset_macaque_peak_matrix <- subset_macaque_peak_matrix[peak_list,]
  subset_mouse_peak_matrix <- subset_mouse_peak_matrix[peak_list,]
  gc()
  
  #normalize
  subset_human_peak_matrix <- t(t(subset_human_peak_matrix)*median(c(Greenleaf_ATAC_ArchR$ReadsInPeaks,macaque_multiome_ArchR$ReadsInPeaks,mouse_multiome_ArchR$ReadsInPeaks))*10000/Greenleaf_ATAC_ArchR@cellColData[colnames(subset_human_peak_matrix),"ReadsInPeaks"])
  subset_macaque_peak_matrix <- t(t(subset_macaque_peak_matrix)*median(c(Greenleaf_ATAC_ArchR$ReadsInPeaks,macaque_multiome_ArchR$ReadsInPeaks,mouse_multiome_ArchR$ReadsInPeaks))*10000/macaque_multiome_ArchR@cellColData[colnames(subset_macaque_peak_matrix),"ReadsInPeaks"])
  subset_mouse_peak_matrix <- t(t(subset_mouse_peak_matrix)*median(c(Greenleaf_ATAC_ArchR$ReadsInPeaks,macaque_multiome_ArchR$ReadsInPeaks,mouse_multiome_ArchR$ReadsInPeaks))*10000/mouse_multiome_ArchR@cellColData[colnames(subset_mouse_peak_matrix),"ReadsInPeaks"])
  gc()
  
  #run wilcox
  HR_df <- my_sparseMatWilcoxon(mat1 = subset_human_peak_matrix,mat2 = subset_macaque_peak_matrix)
  gc()
  HM_df <- my_sparseMatWilcoxon(mat1 = subset_human_peak_matrix,mat2 = subset_mouse_peak_matrix)
  gc()
  RM_df <- my_sparseMatWilcoxon(mat1 = subset_macaque_peak_matrix,mat2 = subset_mouse_peak_matrix)
  gc()
  
  #human specific
  res_HR <- rownames(HR_df)[HR_df$log2FC > 1 & HR_df$fdr < 0.01]
  res_HM <- rownames(HM_df)[HM_df$log2FC > 1 & HM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  HGARs <- peak_list
  
  #primate conserved
  res_HR <- rownames(HR_df)[abs(HR_df$log2FC) <= 1]
  res_HM <- rownames(HM_df)[HM_df$log2FC > 1 & HM_df$fdr < 0.01]
  res_RM <- rownames(RM_df)[RM_df$log2FC > 1 & RM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM,res_RM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  PCARs <- peak_list
  
  #human lost
  res_HR <- rownames(HR_df)[HR_df$log2FC < -1 & HR_df$fdr < 0.01]
  res_HM <- rownames(HM_df)[HM_df$log2FC < -1 & HM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM))
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$mouse_coord,'GRanges'),subject = mouse_peaks) > 0]
  HLARs <- peak_list
  
  #species conserved
  res_HR <- rownames(HR_df)[abs(HR_df$log2FC) <= 1]
  res_HM <- rownames(HM_df)[abs(HM_df$log2FC) <= 1]
  res_RM <- rownames(RM_df)[abs(RM_df$log2FC) <= 1]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM,res_RM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$mouse_coord,'GRanges'),subject = mouse_peaks) > 0]
  CARs <- peak_list
  
  #save list
  HR_df$peak <- rownames(HR_df)
  HR_df$compare <- 'HR'
  HM_df$peak <- rownames(HM_df)
  HM_df$compare <- 'HM'
  RM_df$peak <- rownames(RM_df)
  RM_df$compare <- 'RM'
  
  HR_df$group <- 'others'
  HR_df[which(rownames(HR_df) %in% HGARs),"group"] <- 'HGARs'
  HR_df[which(rownames(HR_df) %in% PCARs),"group"] <- 'PCARs'
  HR_df[which(rownames(HR_df) %in% HLARs),"group"] <- 'HLARs'
  HR_df[which(rownames(HR_df) %in% CARs),"group"] <- 'CARs'
  
  HM_df$group <- 'others'
  HM_df[which(rownames(HM_df) %in% HGARs),"group"] <- 'HGARs'
  HM_df[which(rownames(HM_df) %in% PCARs),"group"] <- 'PCARs'
  HM_df[which(rownames(HM_df) %in% HLARs),"group"] <- 'HLARs'
  HM_df[which(rownames(HM_df) %in% CARs),"group"] <- 'CARs'
  
  RM_df$group <- 'others'
  RM_df[which(rownames(RM_df) %in% HGARs),"group"] <- 'HGARs'
  RM_df[which(rownames(RM_df) %in% PCARs),"group"] <- 'PCARs'
  RM_df[which(rownames(RM_df) %in% HLARs),"group"] <- 'HLARs'
  RM_df[which(rownames(RM_df) %in% CARs),"group"] <- 'CARs'
  
  wilcox_df <- rbind(HR_df,HM_df,RM_df)
  rownames(wilcox_df) <- NULL
  
  char_cmd <- paste0('./res/step_118_fig_230422/Wilcox/',cell_type_dot,'_df_wilcox.rds')
  saveRDS(object = wilcox_df,file = char_cmd)
  gc()
  print(paste(cell_type,'done',sep = ' '))
}

# overlap with previous results -------------------------------------------
cell_type <- 'Ex-4'

cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
old_list <- list.files('./res/step_74_fig_221104')
old_list <- old_list[grep(pattern = 'wilcox',x = old_list,fixed = TRUE)]
old_list <- old_list[grep(pattern = cell_type_dot,x = old_list,fixed = TRUE)]
old_list <- paste('./res/step_74_fig_221104',old_list,sep = '/')
old_list <- readRDS(old_list)

new_list <- list.files('./res/step_118_fig_230422/Wilcox')
new_list <- new_list[grep(pattern = cell_type_dot,x = new_list,fixed = TRUE)]
new_list <- paste('./res/step_118_fig_230422/Wilcox',new_list,sep = '/')
new_list <- readRDS(new_list)

ggvenn(data = list(old = names(old_list)[old_list == 'human_specific'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'HGARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'human_macaque_conserved'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'PCARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'human_loss'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'HLARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'species_conserved'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'CARs']))))

# DESeq2 random ------------------------------------------------------------------
#load data
consensus_peakset <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')
human_peak_matrix <- readRDS(file = './res/step_118_fig_230422/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_118_fig_230422/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_118_fig_230422/mouse_peak_matrix.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

cell_type_list <- base::Reduce(f = base::intersect,x = list(human = as.character(Greenleaf_ATAC_ArchR$cell_type),
                                                            macaque = as.character(macaque_multiome_ArchR$cell_type),
                                                            mouse = as.character(mouse_multiome_ArchR$Gex_macaque_cell_type)))

#for loop
for (cell_type in cell_type_list) {
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load cell type peaks
  human_peaks <- list.files(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls')
  human_peaks <- human_peaks[grep(pattern = cell_type_dot,x = human_peaks,fixed = TRUE)]
  human_peaks <- paste('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls',human_peaks,sep = '/')
  human_peaks <- readRDS(human_peaks)
  
  macaque_peaks <- list.files(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls')
  macaque_peaks <- macaque_peaks[grep(pattern = cell_type_dot,x = macaque_peaks,fixed = TRUE)]
  macaque_peaks <- paste('./processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls',macaque_peaks,sep = '/')
  macaque_peaks <- readRDS(macaque_peaks)
  
  mouse_peaks <- list.files(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls')
  mouse_peaks <- mouse_peaks[grep(pattern = cell_type_dot,x = mouse_peaks,fixed = TRUE)]
  mouse_peaks <- paste('./processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls',mouse_peaks,sep = '/')
  mouse_peaks <- readRDS(mouse_peaks)
  
  #get consensus peaks that overlap with original cell type peaks in at least one species
  human_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$human_coord,'GRanges'),subject = human_peaks))
  macaque_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$macaque_coord,'GRanges'),subject = macaque_peaks))
  mouse_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$human_coord,'GRanges'),subject = mouse_peaks))
  peak_list <- consensus_peakset$human_coord[which((human_overlap + macaque_overlap + mouse_overlap) > 0)]
  
  #get subsetted peak matrix
  subset_human_peak_matrix <- human_peak_matrix[,rownames(Greenleaf_ATAC_ArchR@cellColData)[Greenleaf_ATAC_ArchR$cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$human_coord
  rownames(subset_human_peak_matrix) <- temp[rownames(subset_human_peak_matrix)]$human_coord
  
  subset_macaque_peak_matrix <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$macaque_coord
  rownames(subset_macaque_peak_matrix) <- temp[rownames(subset_macaque_peak_matrix)]$human_coord
  
  subset_mouse_peak_matrix <- mouse_peak_matrix[,rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Gex_macaque_cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$mouse_coord
  rownames(subset_mouse_peak_matrix) <- temp[rownames(subset_mouse_peak_matrix)]$human_coord
  
  subset_human_peak_matrix <- subset_human_peak_matrix[peak_list,]
  subset_macaque_peak_matrix <- subset_macaque_peak_matrix[peak_list,]
  subset_mouse_peak_matrix <- subset_mouse_peak_matrix[peak_list,]
  gc()
  
  # #sampled by donors
  # human_donor_list <- unique(Greenleaf_ATAC_ArchR$Sample)
  # temp <- rownames(subset_human_peak_matrix)
  # subset_human_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = human_donor_list,FUN = function(x){
  #   aggregated_peak_matrix <- rowSums(subset_human_peak_matrix[,which(Greenleaf_ATAC_ArchR@cellColData[colnames(subset_human_peak_matrix),"Sample"] == x)])
  #   return(aggregated_peak_matrix)
  # }))
  # colnames(subset_human_peak_matrix) <- human_donor_list
  # rownames(subset_human_peak_matrix) <- temp
  # subset_human_peak_matrix <- subset_human_peak_matrix[,colSums(subset_human_peak_matrix) > 0]
  # 
  # macaque_donor_list <- unique(macaque_multiome_ArchR$Sample)
  # temp <- rownames(subset_macaque_peak_matrix)
  # subset_macaque_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = macaque_donor_list,FUN = function(x){
  #   aggregated_peak_matrix <- rowSums(subset_macaque_peak_matrix[,which(macaque_multiome_ArchR@cellColData[colnames(subset_macaque_peak_matrix),"Sample"] == x)])
  #   return(aggregated_peak_matrix)
  # }))
  # colnames(subset_macaque_peak_matrix) <- macaque_donor_list
  # rownames(subset_macaque_peak_matrix) <- temp
  # subset_macaque_peak_matrix <- subset_macaque_peak_matrix[,colSums(subset_macaque_peak_matrix) > 0]
  # 
  # mouse_donor_list <- unique(mouse_multiome_ArchR$Sample)
  # temp <- rownames(subset_mouse_peak_matrix)
  # subset_mouse_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = mouse_donor_list,FUN = function(x){
  #   aggregated_peak_matrix <- rowSums(subset_mouse_peak_matrix[,which(mouse_multiome_ArchR@cellColData[colnames(subset_mouse_peak_matrix),"Sample"] == x)])
  #   return(aggregated_peak_matrix)
  # }))
  # colnames(subset_mouse_peak_matrix) <- mouse_donor_list
  # rownames(subset_mouse_peak_matrix) <- temp
  # subset_mouse_peak_matrix <- subset_mouse_peak_matrix[,colSums(subset_mouse_peak_matrix) > 0]
  
  # sampled randomly
  tag_list <- random_equal_split(x = colnames(subset_human_peak_matrix),n = 3)
  temp <- rownames(subset_human_peak_matrix)
  subset_human_peak_matrix <- base::do.call(cbind,args = base::lapply(X = names(tag_list),FUN = function(x){
    aggregated_peak_matrix <- rowSums(subset_human_peak_matrix[,tag_list[[x]]])
    return(aggregated_peak_matrix)
  }))
  rownames(subset_human_peak_matrix) <- temp
  colnames(subset_human_peak_matrix) <- paste('human',names(tag_list),sep = '_')
  
  tag_list <- random_equal_split(x = colnames(subset_macaque_peak_matrix),n = 3)
  temp <- rownames(subset_macaque_peak_matrix)
  subset_macaque_peak_matrix <- base::do.call(cbind,args = base::lapply(X = names(tag_list),FUN = function(x){
    aggregated_peak_matrix <- rowSums(subset_macaque_peak_matrix[,tag_list[[x]]])
    return(aggregated_peak_matrix)
  }))
  rownames(subset_macaque_peak_matrix) <- temp
  colnames(subset_macaque_peak_matrix) <- paste('macaque',names(tag_list),sep = '_')
  
  tag_list <- random_equal_split(x = colnames(subset_mouse_peak_matrix),n = 3)
  temp <- rownames(subset_mouse_peak_matrix)
  subset_mouse_peak_matrix <- base::do.call(cbind,args = base::lapply(X = names(tag_list),FUN = function(x){
    aggregated_peak_matrix <- rowSums(subset_mouse_peak_matrix[,tag_list[[x]]])
    return(aggregated_peak_matrix)
  }))
  rownames(subset_mouse_peak_matrix) <- temp
  colnames(subset_mouse_peak_matrix) <- paste('mouse',names(tag_list),sep = '_')
  
  #create DESeq2
  dds_HR <- cbind(subset_human_peak_matrix,subset_macaque_peak_matrix)
  meta_data <- data.frame(sample = colnames(dds_HR),species = c(rep('human',times = ncol(subset_human_peak_matrix)),
                                                                rep('macaque',times = ncol(subset_macaque_peak_matrix))))
  rownames(meta_data) <- meta_data$sample
  dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
  dds_HR <- DESeq(object = dds_HR)
  
  dds_HM <- cbind(subset_human_peak_matrix,subset_mouse_peak_matrix)
  meta_data <- data.frame(sample = colnames(dds_HM),species = c(rep('human',times = ncol(subset_human_peak_matrix)),
                                                                rep('mouse',times = ncol(subset_mouse_peak_matrix))))
  rownames(meta_data) <- meta_data$sample
  dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
  dds_HM <- DESeq(object = dds_HM)
  
  dds_RM <- cbind(subset_macaque_peak_matrix,subset_mouse_peak_matrix)
  meta_data <- data.frame(sample = colnames(dds_RM),species = c(rep('macaque',times = ncol(subset_macaque_peak_matrix)),
                                                                rep('mouse',times = ncol(subset_mouse_peak_matrix))))
  rownames(meta_data) <- meta_data$sample
  dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
  dds_RM <- DESeq(object = dds_RM)
  
  #DESeq2 df
  HR_df <- results(object = dds_HR,contrast = c('species','human','macaque'))
  HR_df <- na.omit(HR_df)
  colnames(HR_df) <- gsub(pattern = 'log2FoldChange',replacement = 'log2FC',x = colnames(HR_df),fixed = TRUE)
  colnames(HR_df) <- gsub(pattern = 'padj',replacement = 'fdr',x = colnames(HR_df),fixed = TRUE)
  
  HM_df <- results(object = dds_HM,contrast = c('species','human','mouse'))
  HM_df <- na.omit(HM_df)
  colnames(HM_df) <- gsub(pattern = 'log2FoldChange',replacement = 'log2FC',x = colnames(HM_df),fixed = TRUE)
  colnames(HM_df) <- gsub(pattern = 'padj',replacement = 'fdr',x = colnames(HM_df),fixed = TRUE)
  
  RM_df <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
  RM_df <- na.omit(RM_df)
  colnames(RM_df) <- gsub(pattern = 'log2FoldChange',replacement = 'log2FC',x = colnames(RM_df),fixed = TRUE)
  colnames(RM_df) <- gsub(pattern = 'padj',replacement = 'fdr',x = colnames(RM_df),fixed = TRUE)
  
  #human specific
  res_HR <- rownames(HR_df)[HR_df$log2FC > 1 & HR_df$fdr < 0.01]
  res_HM <- rownames(HM_df)[HM_df$log2FC > 1 & HM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  HGARs <- peak_list
  
  #primate conserved
  res_HR <- rownames(HR_df)[abs(HR_df$log2FC) <= 1]
  res_HM <- rownames(HM_df)[HM_df$log2FC > 1 & HM_df$fdr < 0.01]
  res_RM <- rownames(RM_df)[RM_df$log2FC > 1 & RM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM,res_RM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  PCARs <- peak_list
  
  #human lost
  res_HR <- rownames(HR_df)[HR_df$log2FC < -1 & HR_df$fdr < 0.01]
  res_HM <- rownames(HM_df)[HM_df$log2FC < -1 & HM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM))
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$mouse_coord,'GRanges'),subject = mouse_peaks) > 0]
  HLARs <- peak_list
  
  #species conserved
  res_HR <- rownames(HR_df)[abs(HR_df$log2FC) <= 1]
  res_HM <- rownames(HM_df)[abs(HM_df$log2FC) <= 1]
  res_RM <- rownames(RM_df)[abs(RM_df$log2FC) <= 1]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM,res_RM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$mouse_coord,'GRanges'),subject = mouse_peaks) > 0]
  CARs <- peak_list
  
  #save list
  HR_df <- as.data.frame(HR_df)
  HM_df <- as.data.frame(HM_df)
  RM_df <- as.data.frame(RM_df)
  
  HR_df$peak <- rownames(HR_df)
  HR_df$compare <- 'HR'
  HM_df$peak <- rownames(HM_df)
  HM_df$compare <- 'HM'
  RM_df$peak <- rownames(RM_df)
  RM_df$compare <- 'RM'
  
  HR_df$group <- 'others'
  HR_df[which(rownames(HR_df) %in% HGARs),"group"] <- 'HGARs'
  HR_df[which(rownames(HR_df) %in% PCARs),"group"] <- 'PCARs'
  HR_df[which(rownames(HR_df) %in% HLARs),"group"] <- 'HLARs'
  HR_df[which(rownames(HR_df) %in% CARs),"group"] <- 'CARs'
  
  HM_df$group <- 'others'
  HM_df[which(rownames(HM_df) %in% HGARs),"group"] <- 'HGARs'
  HM_df[which(rownames(HM_df) %in% PCARs),"group"] <- 'PCARs'
  HM_df[which(rownames(HM_df) %in% HLARs),"group"] <- 'HLARs'
  HM_df[which(rownames(HM_df) %in% CARs),"group"] <- 'CARs'
  
  RM_df$group <- 'others'
  RM_df[which(rownames(RM_df) %in% HGARs),"group"] <- 'HGARs'
  RM_df[which(rownames(RM_df) %in% PCARs),"group"] <- 'PCARs'
  RM_df[which(rownames(RM_df) %in% HLARs),"group"] <- 'HLARs'
  RM_df[which(rownames(RM_df) %in% CARs),"group"] <- 'CARs'
  
  wilcox_df <- rbind(HR_df,HM_df,RM_df)
  rownames(wilcox_df) <- NULL
  
  char_cmd <- paste0('./res/step_118_fig_230422/DESeq2_random/',cell_type_dot,'_df_DESeq2.rds')
  saveRDS(object = wilcox_df,file = char_cmd)
  gc()
  print(paste(cell_type,'done',sep = ' '))
}

# overlap with previous results -------------------------------------------
cell_type <- 'RG-1'

cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
old_list <- list.files('./res/step_74_fig_221104')
old_list <- old_list[grep(pattern = 'DESeq2',x = old_list,fixed = TRUE)]
old_list <- old_list[grep(pattern = cell_type_dot,x = old_list,fixed = TRUE)]
old_list <- paste('./res/step_74_fig_221104',old_list,sep = '/')
old_list <- readRDS(old_list)

new_list <- list.files('./res/step_118_fig_230422/DESeq2_random')
new_list <- new_list[grep(pattern = cell_type_dot,x = new_list,fixed = TRUE)]
new_list <- paste('./res/step_118_fig_230422/DESeq2_random',new_list,sep = '/')
new_list <- readRDS(new_list)

ggvenn(data = list(old = names(old_list)[old_list == 'human_specific'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'HGARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'human_macaque_conserved'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'PCARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'human_loss'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'HLARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'species_conserved'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'CARs']))))

# DESeq2 donor ------------------------------------------------------------------
#load data
consensus_peakset <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')
human_peak_matrix <- readRDS(file = './res/step_118_fig_230422/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_118_fig_230422/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_118_fig_230422/mouse_peak_matrix.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

cell_type_list <- base::Reduce(f = base::intersect,x = list(human = as.character(Greenleaf_ATAC_ArchR$cell_type),
                                                            macaque = as.character(macaque_multiome_ArchR$cell_type),
                                                            mouse = as.character(mouse_multiome_ArchR$Gex_macaque_cell_type)))

#for loop
for (cell_type in cell_type_list) {
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load cell type peaks
  human_peaks <- list.files(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls')
  human_peaks <- human_peaks[grep(pattern = cell_type_dot,x = human_peaks,fixed = TRUE)]
  human_peaks <- paste('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls',human_peaks,sep = '/')
  human_peaks <- readRDS(human_peaks)
  
  macaque_peaks <- list.files(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls')
  macaque_peaks <- macaque_peaks[grep(pattern = cell_type_dot,x = macaque_peaks,fixed = TRUE)]
  macaque_peaks <- paste('./processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls',macaque_peaks,sep = '/')
  macaque_peaks <- readRDS(macaque_peaks)
  
  mouse_peaks <- list.files(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls')
  mouse_peaks <- mouse_peaks[grep(pattern = cell_type_dot,x = mouse_peaks,fixed = TRUE)]
  mouse_peaks <- paste('./processed_data/221008_summary/mouse_multiome_ArchR_221009/PeakCalls',mouse_peaks,sep = '/')
  mouse_peaks <- readRDS(mouse_peaks)
  
  #get consensus peaks that overlap with original cell type peaks in at least one species
  human_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$human_coord,'GRanges'),subject = human_peaks))
  macaque_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$macaque_coord,'GRanges'),subject = macaque_peaks))
  mouse_overlap <- as.numeric(countOverlaps(query = as(consensus_peakset$human_coord,'GRanges'),subject = mouse_peaks))
  peak_list <- consensus_peakset$human_coord[which((human_overlap + macaque_overlap + mouse_overlap) > 0)]
  
  #get subsetted peak matrix
  subset_human_peak_matrix <- human_peak_matrix[,rownames(Greenleaf_ATAC_ArchR@cellColData)[Greenleaf_ATAC_ArchR$cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$human_coord
  rownames(subset_human_peak_matrix) <- temp[rownames(subset_human_peak_matrix)]$human_coord
  
  subset_macaque_peak_matrix <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$macaque_coord
  rownames(subset_macaque_peak_matrix) <- temp[rownames(subset_macaque_peak_matrix)]$human_coord
  
  subset_mouse_peak_matrix <- mouse_peak_matrix[,rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Gex_macaque_cell_type == cell_type]]
  temp <- consensus_peakset
  names(temp) <- temp$mouse_coord
  rownames(subset_mouse_peak_matrix) <- temp[rownames(subset_mouse_peak_matrix)]$human_coord
  
  subset_human_peak_matrix <- subset_human_peak_matrix[peak_list,]
  subset_macaque_peak_matrix <- subset_macaque_peak_matrix[peak_list,]
  subset_mouse_peak_matrix <- subset_mouse_peak_matrix[peak_list,]
  gc()
  
  #sampled by donors
  human_donor_list <- unique(Greenleaf_ATAC_ArchR$Sample)
  temp <- rownames(subset_human_peak_matrix)
  subset_human_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = human_donor_list,FUN = function(x){
    aggregated_peak_matrix <- rowSums(subset_human_peak_matrix[,which(Greenleaf_ATAC_ArchR@cellColData[colnames(subset_human_peak_matrix),"Sample"] == x)])
    return(aggregated_peak_matrix)
  }))
  colnames(subset_human_peak_matrix) <- human_donor_list
  rownames(subset_human_peak_matrix) <- temp
  subset_human_peak_matrix <- subset_human_peak_matrix[,colSums(subset_human_peak_matrix) > 0]
  
  macaque_donor_list <- unique(macaque_multiome_ArchR$Sample)
  temp <- rownames(subset_macaque_peak_matrix)
  subset_macaque_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = macaque_donor_list,FUN = function(x){
    aggregated_peak_matrix <- rowSums(subset_macaque_peak_matrix[,which(macaque_multiome_ArchR@cellColData[colnames(subset_macaque_peak_matrix),"Sample"] == x)])
    return(aggregated_peak_matrix)
  }))
  colnames(subset_macaque_peak_matrix) <- macaque_donor_list
  rownames(subset_macaque_peak_matrix) <- temp
  subset_macaque_peak_matrix <- subset_macaque_peak_matrix[,colSums(subset_macaque_peak_matrix) > 0]
  
  mouse_donor_list <- unique(mouse_multiome_ArchR$Sample)
  temp <- rownames(subset_mouse_peak_matrix)
  subset_mouse_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = mouse_donor_list,FUN = function(x){
    aggregated_peak_matrix <- rowSums(subset_mouse_peak_matrix[,which(mouse_multiome_ArchR@cellColData[colnames(subset_mouse_peak_matrix),"Sample"] == x)])
    return(aggregated_peak_matrix)
  }))
  colnames(subset_mouse_peak_matrix) <- mouse_donor_list
  rownames(subset_mouse_peak_matrix) <- temp
  subset_mouse_peak_matrix <- subset_mouse_peak_matrix[,colSums(subset_mouse_peak_matrix) > 0]
  
  #create DESeq2
  dds_HR <- cbind(subset_human_peak_matrix,subset_macaque_peak_matrix)
  meta_data <- data.frame(sample = colnames(dds_HR),species = c(rep('human',times = ncol(subset_human_peak_matrix)),
                                                                rep('macaque',times = ncol(subset_macaque_peak_matrix))))
  rownames(meta_data) <- meta_data$sample
  dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
  dds_HR <- DESeq(object = dds_HR)
  
  dds_HM <- cbind(subset_human_peak_matrix,subset_mouse_peak_matrix)
  meta_data <- data.frame(sample = colnames(dds_HM),species = c(rep('human',times = ncol(subset_human_peak_matrix)),
                                                                rep('mouse',times = ncol(subset_mouse_peak_matrix))))
  rownames(meta_data) <- meta_data$sample
  dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
  dds_HM <- DESeq(object = dds_HM)
  
  dds_RM <- cbind(subset_macaque_peak_matrix,subset_mouse_peak_matrix)
  meta_data <- data.frame(sample = colnames(dds_RM),species = c(rep('macaque',times = ncol(subset_macaque_peak_matrix)),
                                                                rep('mouse',times = ncol(subset_mouse_peak_matrix))))
  rownames(meta_data) <- meta_data$sample
  dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
  dds_RM <- DESeq(object = dds_RM)
  
  #DESeq2 df
  HR_df <- results(object = dds_HR,contrast = c('species','human','macaque'))
  HR_df <- na.omit(HR_df)
  colnames(HR_df) <- gsub(pattern = 'log2FoldChange',replacement = 'log2FC',x = colnames(HR_df),fixed = TRUE)
  colnames(HR_df) <- gsub(pattern = 'padj',replacement = 'fdr',x = colnames(HR_df),fixed = TRUE)
  
  HM_df <- results(object = dds_HM,contrast = c('species','human','mouse'))
  HM_df <- na.omit(HM_df)
  colnames(HM_df) <- gsub(pattern = 'log2FoldChange',replacement = 'log2FC',x = colnames(HM_df),fixed = TRUE)
  colnames(HM_df) <- gsub(pattern = 'padj',replacement = 'fdr',x = colnames(HM_df),fixed = TRUE)
  
  RM_df <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
  RM_df <- na.omit(RM_df)
  colnames(RM_df) <- gsub(pattern = 'log2FoldChange',replacement = 'log2FC',x = colnames(RM_df),fixed = TRUE)
  colnames(RM_df) <- gsub(pattern = 'padj',replacement = 'fdr',x = colnames(RM_df),fixed = TRUE)
  
  #human specific
  res_HR <- rownames(HR_df)[HR_df$log2FC > 1 & HR_df$fdr < 0.01]
  res_HM <- rownames(HM_df)[HM_df$log2FC > 1 & HM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  HGARs <- peak_list
  
  #primate conserved
  res_HR <- rownames(HR_df)[abs(HR_df$log2FC) <= 1]
  res_HM <- rownames(HM_df)[HM_df$log2FC > 1 & HM_df$fdr < 0.01]
  res_RM <- rownames(RM_df)[RM_df$log2FC > 1 & RM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM,res_RM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  PCARs <- peak_list
  
  #human lost
  res_HR <- rownames(HR_df)[HR_df$log2FC < -1 & HR_df$fdr < 0.01]
  res_HM <- rownames(HM_df)[HM_df$log2FC < -1 & HM_df$fdr < 0.01]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM))
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$mouse_coord,'GRanges'),subject = mouse_peaks) > 0]
  HLARs <- peak_list
  
  #species conserved
  res_HR <- rownames(HR_df)[abs(HR_df$log2FC) <= 1]
  res_HM <- rownames(HM_df)[abs(HM_df$log2FC) <= 1]
  res_RM <- rownames(RM_df)[abs(RM_df$log2FC) <= 1]
  peak_list <- base::Reduce(f = base::intersect,x = list(res_HR,res_HM,res_RM))
  peak_list <- peak_list[countOverlaps(query = as(peak_list,'GRanges'),subject = human_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$macaque_coord,'GRanges'),subject = macaque_peaks) > 0]
  peak_list <- peak_list[countOverlaps(query = as(consensus_peakset[peak_list]$mouse_coord,'GRanges'),subject = mouse_peaks) > 0]
  CARs <- peak_list
  
  #save list
  HR_df <- as.data.frame(HR_df)
  HM_df <- as.data.frame(HM_df)
  RM_df <- as.data.frame(RM_df)
  
  HR_df$peak <- rownames(HR_df)
  HR_df$compare <- 'HR'
  HM_df$peak <- rownames(HM_df)
  HM_df$compare <- 'HM'
  RM_df$peak <- rownames(RM_df)
  RM_df$compare <- 'RM'
  
  HR_df$group <- 'others'
  HR_df[which(rownames(HR_df) %in% HGARs),"group"] <- 'HGARs'
  HR_df[which(rownames(HR_df) %in% PCARs),"group"] <- 'PCARs'
  HR_df[which(rownames(HR_df) %in% HLARs),"group"] <- 'HLARs'
  HR_df[which(rownames(HR_df) %in% CARs),"group"] <- 'CARs'
  
  HM_df$group <- 'others'
  HM_df[which(rownames(HM_df) %in% HGARs),"group"] <- 'HGARs'
  HM_df[which(rownames(HM_df) %in% PCARs),"group"] <- 'PCARs'
  HM_df[which(rownames(HM_df) %in% HLARs),"group"] <- 'HLARs'
  HM_df[which(rownames(HM_df) %in% CARs),"group"] <- 'CARs'
  
  RM_df$group <- 'others'
  RM_df[which(rownames(RM_df) %in% HGARs),"group"] <- 'HGARs'
  RM_df[which(rownames(RM_df) %in% PCARs),"group"] <- 'PCARs'
  RM_df[which(rownames(RM_df) %in% HLARs),"group"] <- 'HLARs'
  RM_df[which(rownames(RM_df) %in% CARs),"group"] <- 'CARs'
  
  wilcox_df <- rbind(HR_df,HM_df,RM_df)
  rownames(wilcox_df) <- NULL
  
  char_cmd <- paste0('./res/step_118_fig_230422/DESeq2_donor/',cell_type_dot,'_df_DESeq2.rds')
  saveRDS(object = wilcox_df,file = char_cmd)
  gc()
  print(paste(cell_type,'done',sep = ' '))
}

# overlap with previous results -------------------------------------------
cell_type <- 'Ex-1'

cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
old_list <- list.files('./res/step_74_fig_221104')
old_list <- old_list[grep(pattern = 'DESeq2',x = old_list,fixed = TRUE)]
old_list <- old_list[grep(pattern = cell_type_dot,x = old_list,fixed = TRUE)]
old_list <- paste('./res/step_74_fig_221104',old_list,sep = '/')
old_list <- readRDS(old_list)

new_list <- list.files('./res/step_118_fig_230422/DESeq2_donor')
new_list <- new_list[grep(pattern = cell_type_dot,x = new_list,fixed = TRUE)]
new_list <- paste('./res/step_118_fig_230422/DESeq2_donor',new_list,sep = '/')
new_list <- readRDS(new_list)

ggvenn(data = list(old = names(old_list)[old_list == 'human_specific'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'HGARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'human_macaque_conserved'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'PCARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'human_loss'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'HLARs']))))
ggvenn(data = list(old = names(old_list)[old_list == 'species_conserved'],
                   new = unique(gsub(pattern = ':',replacement = '-',x = new_list$peak[new_list$group == 'CARs']))))

# overlap random with donor -----------------------------------------------
cell_type <- 'IP'

donor_list <- list.files('./res/step_118_fig_230422/DESeq2_donor')
donor_list <- donor_list[grep(pattern = cell_type_dot,x = donor_list,fixed = TRUE)]
donor_list <- paste('./res/step_118_fig_230422/DESeq2_donor',donor_list,sep = '/')
donor_list <- readRDS(donor_list)

random_list <- list.files('./res/step_118_fig_230422/DESeq2_random')
random_list <- random_list[grep(pattern = cell_type_dot,x = random_list,fixed = TRUE)]
random_list <- paste('./res/step_118_fig_230422/DESeq2_random',random_list,sep = '/')
random_list <- readRDS(random_list)

ggvenn(data = list(old = donor_list$peak[donor_list$group == 'HGARs'],
                   new = random_list$peak[random_list$group == 'HGARs']))
ggvenn(data = list(old = donor_list$peak[donor_list$group == 'PCARs'],
                   new = random_list$peak[random_list$group == 'PCARs']))
ggvenn(data = list(old = donor_list$peak[donor_list$group == 'HLARs'],
                   new = random_list$peak[random_list$group == 'HLARs']))
ggvenn(data = list(old = donor_list$peak[donor_list$group == 'CARs'],
                   new = random_list$peak[random_list$group == 'CARs']))

# overlap donor with Wilcox -----------------------------------------------
cell_type <- 'RG-1'

donor_list <- list.files('./res/step_118_fig_230422/DESeq2_donor')
donor_list <- donor_list[grep(pattern = cell_type_dot,x = donor_list,fixed = TRUE)]
donor_list <- paste('./res/step_118_fig_230422/DESeq2_donor',donor_list,sep = '/')
donor_list <- readRDS(donor_list)

random_list <- list.files('./res/step_118_fig_230422/Wilcox')
random_list <- random_list[grep(pattern = cell_type_dot,x = random_list,fixed = TRUE)]
random_list <- paste('./res/step_118_fig_230422/Wilcox',random_list,sep = '/')
random_list <- readRDS(random_list)

ggvenn(data = list(old = donor_list$peak[donor_list$group == 'HGARs'],
                   new = random_list$peak[random_list$group == 'HGARs']))
ggvenn(data = list(old = donor_list$peak[donor_list$group == 'PCARs'],
                   new = random_list$peak[random_list$group == 'PCARs']))
ggvenn(data = list(old = donor_list$peak[donor_list$group == 'HLARs'],
                   new = random_list$peak[random_list$group == 'HLARs']))
ggvenn(data = list(old = donor_list$peak[donor_list$group == 'CARs'],
                   new = random_list$peak[random_list$group == 'CARs']))

# Wilcox enrichheatmap ----------------------------------------------------
#load data
consensus_peakset <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

cell_type_list <- base::Reduce(f = base::intersect,x = list(human = as.character(Greenleaf_ATAC_ArchR$cell_type),
                                                            macaque = as.character(macaque_multiome_ArchR$cell_type),
                                                            mouse = as.character(mouse_multiome_ArchR$Gex_macaque_cell_type)))

for (cell_type in cell_type_list) {
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load DAP
  DAP_list <- list.files(path = './res/step_118_fig_230422/Wilcox')
  DAP_list <- DAP_list[grep(pattern = cell_type_dot,x = DAP_list,fixed = TRUE)]
  DAP_list <- paste('./res/step_118_fig_230422/Wilcox',DAP_list,sep = '/')
  DAP_list <- readRDS(DAP_list)
  DAP_list <- DAP_list[DAP_list$group != 'others',]
  
  #load coverage
  human_coverage <- list.files('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type')
  human_coverage <- human_coverage[grep(pattern = cell_type_dot,x = human_coverage,fixed = TRUE)]
  human_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type',human_coverage,sep = '/')
  human_coverage <- rtracklayer::import.bw(con = human_coverage)
  
  macaque_coverage <- list.files('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type')
  macaque_coverage <- macaque_coverage[grep(pattern = cell_type_dot,x = macaque_coverage,fixed = TRUE)]
  macaque_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type',macaque_coverage,sep = '/')
  macaque_coverage <- rtracklayer::import.bw(con = macaque_coverage)
  
  mouse_coverage <- list.files('./res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type')
  mouse_coverage <- mouse_coverage[grep(pattern = cell_type_dot,x = mouse_coverage,fixed = TRUE)]
  mouse_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type',mouse_coverage,sep = '/')
  mouse_coverage <- rtracklayer::import.bw(con = mouse_coverage)
  
  #generate group list
  group_list <- base::lapply(X = unique(DAP_list$group),FUN = function(x){
    peak_list <- unique(DAP_list[DAP_list$group == x,"peak"])
    temp <- rep(x,times = length(peak_list))
    names(temp) <- peak_list
    return(temp)
  })
  group_list <- base::Reduce(f = append,x = group_list)
  
  #human enrichheatmap
  coord_list <- consensus_peakset[names(group_list)]$human_coord
  coord_list <- base::lapply(X = coord_list,FUN = function(x){
    temp <- strsplit(x = x,split = ':')[[1]]
    temp_chrom <- temp[1]
    mid_point <- strsplit(x = temp[2],split = '-')[[1]]
    mid_point <- round((as.numeric(mid_point[1]) + as.numeric(mid_point[2])) / 2)
    return(paste(temp_chrom,mid_point,sep = ':'))
  })
  coord_list <- as(unlist(coord_list),'GRanges')
  names(coord_list) <- names(group_list)
  
  mat <- normalizeToMatrix(signal = human_coverage,target = coord_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p_human <- EnrichedHeatmap(mat = mat,row_split = factor(group_list[rownames(mat)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean)
  
  #macaque enrichheatmap
  coord_list <- consensus_peakset[names(group_list)]$macaque_coord
  coord_list <- base::lapply(X = coord_list,FUN = function(x){
    temp <- strsplit(x = x,split = ':')[[1]]
    temp_chrom <- temp[1]
    mid_point <- strsplit(x = temp[2],split = '-')[[1]]
    mid_point <- round((as.numeric(mid_point[1]) + as.numeric(mid_point[2])) / 2)
    return(paste(temp_chrom,mid_point,sep = ':'))
  })
  coord_list <- as(unlist(coord_list),'GRanges')
  names(coord_list) <- names(group_list)
  
  mat <- normalizeToMatrix(signal = macaque_coverage,target = coord_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p_macaque <- EnrichedHeatmap(mat = mat,row_split = factor(group_list[rownames(mat)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean)
  
  #mouse enrichheatmap
  coord_list <- consensus_peakset[names(group_list)]$mouse_coord
  coord_list <- base::lapply(X = coord_list,FUN = function(x){
    temp <- strsplit(x = x,split = ':')[[1]]
    temp_chrom <- temp[1]
    mid_point <- strsplit(x = temp[2],split = '-')[[1]]
    mid_point <- round((as.numeric(mid_point[1]) + as.numeric(mid_point[2])) / 2)
    return(paste(temp_chrom,mid_point,sep = ':'))
  })
  coord_list <- as(unlist(coord_list),'GRanges')
  names(coord_list) <- names(group_list)
  
  mat <- normalizeToMatrix(signal = mouse_coverage,target = coord_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p_mouse <- EnrichedHeatmap(mat = mat,row_split = factor(group_list[rownames(mat)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean)
  
  #plot
  col_fun <- colorRamp2(breaks = c(0,0.01,0.02),colors = c('#FFC30F','#C70039','#581845'))
  p <- EnrichedHeatmap(mat = p_human@matrix,row_split = factor(group_list[rownames(p_human@matrix)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean,col = col_fun,pos_line = FALSE,row_order = p_human@row_order,name = 'human',width = unit(1.5,'inches'),height = unit(9,'inches'),show_heatmap_legend = FALSE) + 
    EnrichedHeatmap(mat = p_macaque@matrix,row_split = factor(group_list[rownames(p_macaque@matrix)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean,col = col_fun,pos_line = FALSE,row_order = p_human@row_order,name = 'macaque',width = unit(1.5,'inches'),height = unit(9,'inches'),show_heatmap_legend = FALSE) + 
    EnrichedHeatmap(mat = p_mouse@matrix,row_split = factor(group_list[rownames(p_mouse@matrix)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean,col = col_fun,pos_line = FALSE,row_order = p_human@row_order,name = 'Normalized insertions',width = unit(1.5,'inches'),height = unit(9,'inches'),show_heatmap_legend = TRUE)
  
  temp_char <- paste0('./res/step_118_fig_230422/Wilcox_enrichHeatmap/',cell_type_dot,'_wilcox_enrichheatmap.pdf')
  pdf(file = temp_char,width = 8,height = 12)
  print(p)
  dev.off()
  
  #echo
  temp_char <- paste(cell_type_dot,'done!',sep = ' ')
  gc()
  print(temp_char)
}

# DESeq2 donor enrichheatmap ----------------------------------------------------
#load data
consensus_peakset <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

cell_type_list <- base::Reduce(f = base::intersect,x = list(human = as.character(Greenleaf_ATAC_ArchR$cell_type),
                                                            macaque = as.character(macaque_multiome_ArchR$cell_type),
                                                            mouse = as.character(mouse_multiome_ArchR$Gex_macaque_cell_type)))
cell_type_list <- cell_type_list[!(cell_type_list %in% c('InMGE'))]

for (cell_type in cell_type_list) {
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load DAP
  DAP_list <- list.files(path = './res/step_118_fig_230422/DESeq2_donor')
  DAP_list <- DAP_list[grep(pattern = cell_type_dot,x = DAP_list,fixed = TRUE)]
  DAP_list <- paste('./res/step_118_fig_230422/DESeq2_donor',DAP_list,sep = '/')
  DAP_list <- readRDS(DAP_list)
  DAP_list <- DAP_list[DAP_list$group != 'others',]
  
  #load coverage
  human_coverage <- list.files('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type')
  human_coverage <- human_coverage[grep(pattern = cell_type_dot,x = human_coverage,fixed = TRUE)]
  human_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type',human_coverage,sep = '/')
  human_coverage <- rtracklayer::import.bw(con = human_coverage)
  
  macaque_coverage <- list.files('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type')
  macaque_coverage <- macaque_coverage[grep(pattern = cell_type_dot,x = macaque_coverage,fixed = TRUE)]
  macaque_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type',macaque_coverage,sep = '/')
  macaque_coverage <- rtracklayer::import.bw(con = macaque_coverage)
  
  mouse_coverage <- list.files('./res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type')
  mouse_coverage <- mouse_coverage[grep(pattern = cell_type_dot,x = mouse_coverage,fixed = TRUE)]
  mouse_coverage <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type',mouse_coverage,sep = '/')
  mouse_coverage <- rtracklayer::import.bw(con = mouse_coverage)
  
  #generate group list
  group_list <- base::lapply(X = unique(DAP_list$group),FUN = function(x){
    peak_list <- unique(DAP_list[DAP_list$group == x,"peak"])
    temp <- rep(x,times = length(peak_list))
    names(temp) <- peak_list
    return(temp)
  })
  group_list <- base::Reduce(f = append,x = group_list)
  
  #human enrichheatmap
  coord_list <- consensus_peakset[names(group_list)]$human_coord
  coord_list <- base::lapply(X = coord_list,FUN = function(x){
    temp <- strsplit(x = x,split = ':')[[1]]
    temp_chrom <- temp[1]
    mid_point <- strsplit(x = temp[2],split = '-')[[1]]
    mid_point <- round((as.numeric(mid_point[1]) + as.numeric(mid_point[2])) / 2)
    return(paste(temp_chrom,mid_point,sep = ':'))
  })
  coord_list <- as(unlist(coord_list),'GRanges')
  names(coord_list) <- names(group_list)
  
  mat <- normalizeToMatrix(signal = human_coverage,target = coord_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p_human <- EnrichedHeatmap(mat = mat,row_split = factor(group_list[rownames(mat)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean)
  
  #macaque enrichheatmap
  coord_list <- consensus_peakset[names(group_list)]$macaque_coord
  coord_list <- base::lapply(X = coord_list,FUN = function(x){
    temp <- strsplit(x = x,split = ':')[[1]]
    temp_chrom <- temp[1]
    mid_point <- strsplit(x = temp[2],split = '-')[[1]]
    mid_point <- round((as.numeric(mid_point[1]) + as.numeric(mid_point[2])) / 2)
    return(paste(temp_chrom,mid_point,sep = ':'))
  })
  coord_list <- as(unlist(coord_list),'GRanges')
  names(coord_list) <- names(group_list)
  
  mat <- normalizeToMatrix(signal = macaque_coverage,target = coord_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p_macaque <- EnrichedHeatmap(mat = mat,row_split = factor(group_list[rownames(mat)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean)
  
  #mouse enrichheatmap
  coord_list <- consensus_peakset[names(group_list)]$mouse_coord
  coord_list <- base::lapply(X = coord_list,FUN = function(x){
    temp <- strsplit(x = x,split = ':')[[1]]
    temp_chrom <- temp[1]
    mid_point <- strsplit(x = temp[2],split = '-')[[1]]
    mid_point <- round((as.numeric(mid_point[1]) + as.numeric(mid_point[2])) / 2)
    return(paste(temp_chrom,mid_point,sep = ':'))
  })
  coord_list <- as(unlist(coord_list),'GRanges')
  names(coord_list) <- names(group_list)
  
  mat <- normalizeToMatrix(signal = mouse_coverage,target = coord_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  p_mouse <- EnrichedHeatmap(mat = mat,row_split = factor(group_list[rownames(mat)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean)
  
  #plot
  col_fun <- colorRamp2(breaks = c(0,0.01,0.02),colors = c('#FFC30F','#C70039','#581845'))
  p <- EnrichedHeatmap(mat = p_human@matrix,row_split = factor(group_list[rownames(p_human@matrix)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean,col = col_fun,pos_line = FALSE,row_order = p_human@row_order,name = 'human',width = unit(1.5,'inches'),height = unit(9.5,'inches'),show_heatmap_legend = FALSE) + 
    EnrichedHeatmap(mat = p_macaque@matrix,row_split = factor(group_list[rownames(p_macaque@matrix)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean,col = col_fun,pos_line = FALSE,row_order = p_human@row_order,name = 'macaque',width = unit(1.5,'inches'),height = unit(9.5,'inches'),show_heatmap_legend = FALSE) + 
    EnrichedHeatmap(mat = p_mouse@matrix,row_split = factor(group_list[rownames(p_mouse@matrix)],levels = c('HGARs','PCARs','HLARs','CARs')),use_raster = TRUE,raster_resize_mat = mean,col = col_fun,pos_line = FALSE,row_order = p_human@row_order,name = 'Normalized insertions',width = unit(1.5,'inches'),height = unit(9.5,'inches'),show_heatmap_legend = TRUE)
  
  temp_char <- paste0('./res/step_118_fig_230422/DESeq2_donor_enrichHeatmap/',cell_type_dot,'_DESeq2_donor_enrichheatmap.pdf')
  pdf(file = temp_char,width = 8,height = 12)
  print(p)
  dev.off()
  
  #echo
  temp_char <- paste(cell_type_dot,'done!',sep = ' ')
  gc()
  print(temp_char)
}
