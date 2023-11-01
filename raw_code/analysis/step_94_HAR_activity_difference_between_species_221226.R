#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: HAR activity difference between species                         ##
## Data: 2022.12.26                                                                ##
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

# redo DF analysis --------------------------------------------------------
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
  
  #human mouse conserved
  temp_1 <- rownames(DF_human_macaque)[DF_human_macaque$log2FC > lfc_cutoff & DF_human_macaque$fdr < fdr_cutoff]
  temp_2 <- rownames(DF_human_mouse)[abs(DF_human_mouse$log2FC) <= lfc_cutoff]
  temp_3 <- rownames(DF_macaque_mouse)[DF_macaque_mouse$log2FC < -1*lfc_cutoff & DF_macaque_mouse$fdr < fdr_cutoff]
  human_mouse_conserved <- dplyr::intersect(x = temp_1,y = temp_2)
  human_mouse_conserved <- dplyr::intersect(x = human_mouse_conserved,y = temp_3)
  
  temp <- human_HAR[human_mouse_conserved]
  temp <- countOverlaps(query = temp,subject = human_peakset)
  human_mouse_conserved <- names(temp)[temp > 0]
  
  names(mouse_HAR) <- mouse_HAR$ori_peak
  temp <- mouse_HAR[human_mouse_conserved]
  temp <- countOverlaps(query = temp,subject = mouse_peakset)
  human_mouse_conserved <- names(temp)[temp > 0]
  
  #macaque specific
  temp_1 <- rownames(DF_human_macaque)[DF_human_macaque$log2FC < -1*lfc_cutoff & DF_human_macaque$fdr < fdr_cutoff]
  temp_2 <- rownames(DF_macaque_mouse)[DF_macaque_mouse$log2FC > lfc_cutoff & DF_macaque_mouse$fdr < fdr_cutoff]
  macaque_specific <- dplyr::intersect(x = temp_1,y = temp_2)
  
  names(macaque_HAR) <- macaque_HAR$ori_peak
  temp <- macaque_HAR[macaque_specific]
  temp <- countOverlaps(query = temp,subject = macaque_peakset)
  macaque_specific <- names(temp)[temp > 0]
  
  #mouse specific
  temp_1 <- rownames(DF_human_mouse)[DF_human_mouse$log2FC < -1*lfc_cutoff & DF_human_mouse$fdr < fdr_cutoff]
  temp_2 <- rownames(DF_macaque_mouse)[DF_macaque_mouse$log2FC < -1*lfc_cutoff & DF_macaque_mouse$fdr < fdr_cutoff]
  mouse_specific <- dplyr::intersect(x = temp_1,y = temp_2)
  
  names(mouse_HAR) <- mouse_HAR$ori_peak
  temp <- mouse_HAR[mouse_specific]
  temp <- countOverlaps(query = temp,subject = mouse_peakset)
  mouse_specific <- names(temp)[temp > 0]
  
  #human loss
  temp_1 <- rownames(DF_human_macaque)[DF_human_macaque$log2FC < -1*lfc_cutoff & DF_human_macaque$fdr < fdr_cutoff]
  temp_2 <- rownames(DF_human_mouse)[DF_human_mouse$log2FC < -1*lfc_cutoff & DF_human_mouse$fdr < fdr_cutoff]
  temp_3 <- rownames(DF_macaque_mouse)[abs(DF_macaque_mouse$log2FC) <= lfc_cutoff]
  human_loss <- dplyr::intersect(x = temp_1,temp_2)
  human_loss <- dplyr::intersect(x = human_loss,y = temp_3)
  
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
                         human_mouse_conserved = human_mouse_conserved,
                         macaque_specific = macaque_specific,
                         mouse_specific = mouse_specific,
                         human_loss = human_loss,
                         species_conserved = species_conserved)
  char <- paste0('./res/step_94_fig_221226/DF_three_species/',cell_type_dot,'_DF_list.rds')
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
  DF_list <- list.files(path = './res/step_94_fig_221226/DF_three_species')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_94_fig_221226/DF_three_species',DF_list,sep = '/')
  DF_list <- readRDS(file = DF_list)
  
  temp <- c(rep('human_specific',times = length(DF_list$human_specific)),
            rep('primate_conserved',times = length(DF_list$primate_conserved)),
            rep('human_mouse_conserved',times = length(DF_list$human_mouse_conserved)),
            rep('macaque_specific',times = length(DF_list$macaque_specific)),
            rep('mouse_specific',times = length(DF_list$mouse_specific)),
            rep('human_loss',times = length(DF_list$human_loss)),
            rep('species_conserved',times = length(DF_list$species_conserved)))
  names(temp) <- c(DF_list$human_specific,DF_list$primate_conserved,DF_list$human_mouse_conserved,DF_list$macaque_specific,
                   DF_list$mouse_specific,DF_list$human_loss,DF_list$species_conserved)
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
  p1 <- EnrichedHeatmap(mat = mat,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_mouse_conserved','macaque_specific','mouse_specific','human_loss','species_conserved')),use_raster = TRUE)
  
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
  p2 <- EnrichedHeatmap(mat = mat,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_mouse_conserved','macaque_specific','mouse_specific','human_loss','species_conserved')),use_raster = TRUE)
  
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
  p3 <- EnrichedHeatmap(mat = mat,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_mouse_conserved','macaque_specific','mouse_specific','human_loss','species_conserved')),use_raster = TRUE)
  
  #plot
  col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
  char <- paste0('./res/step_94_fig_221226/Enrich_heatmap_plot/',cell_type_dot,'_enrich_heatmap.pdf')
  pdf(file = char,width = 3.6,height = 8.5)
  print(EnrichedHeatmap(mat = p1@matrix,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_mouse_conserved','macaque_specific','mouse_specific','human_loss','species_conserved')),
                        use_raster = TRUE,col = col_fun,name = 'insertion',pos_line = FALSE) + 
          EnrichedHeatmap(mat = p2@matrix,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_mouse_conserved','macaque_specific','mouse_specific','human_loss','species_conserved')),
                          use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order) + 
          EnrichedHeatmap(mat = p3@matrix,row_split = factor(DF_list,levels = c('human_specific','primate_conserved','human_mouse_conserved','macaque_specific','mouse_specific','human_loss','species_conserved')),
                          use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order))
  dev.off()
  
  #report
  char <- paste(cell_type_dot,'done!',sep = ' ')
  print(char)
}

# HAR group number --------------------------------------------------------
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

#lapply count HAR group number
count_matrix <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  #cell type
  cell_type <- x
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

  #load HAR_list
  temp_HAR <- list.files(path = './res/step_94_fig_221226/DF_three_species')
  temp_HAR <- temp_HAR[grep(pattern = cell_type_dot,x = temp_HAR,fixed = TRUE)]
  temp_HAR <- paste('./res/step_94_fig_221226/DF_three_species',temp_HAR,sep = '/')
  temp_HAR <- readRDS(file = temp_HAR)

  #count
  group_list <- names(temp_HAR)
  temp <- unlist(base::lapply(X = group_list,FUN = function(i){
    temp <- temp_HAR[[i]]
    return(length(temp))
  }))
  temp <- data.frame(group_list = group_list,number = temp,cell_type = cell_type)
  return(temp)
}))
count_matrix$group_list <- factor(count_matrix$group_list,levels = rev(c('human_specific','primate_conserved','human_mouse_conserved','macaque_specific','mouse_specific','human_loss','species_conserved')))
count_matrix$cell_type <- factor(count_matrix$cell_type,levels = cell_type_list)

#ggplot
pdf(file = './res/step_94_fig_221226/HAR_group_number_barplot.pdf',width = 6,height = 4)
ggplot(data = count_matrix,aes(x = group_list,y = number,fill = cell_type)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_bw() + theme(aspect.ratio = 0.8) + 
  scale_fill_manual(values = color_param$celltype[cell_type_list]) + 
  coord_flip() + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5),
        axis.text = element_text(size = 10)) + 
  xlab('Group') + ylab('number') + labs(title = 'HAR group number')
dev.off()

# HAR related gene expression ---------------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#load gene list
gene_list <- readRDS(file = './res/step_93_fig_221216/gene_list_affected_by_human_specific_HAR.rds')
gene_list <- unlist(gene_list)
names(gene_list) <- NULL
gene_list <- unique(gene_list)

#modify gene list
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

#gene plot
for(i in gene_list){
  #load annotation
  macaque_anno <- readRDS(file = '/home/sunym/temp/macaque_anno.rds')
  mouse_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCh38_to_GRCm39.csv')
  
  #gene_name
  human_gene <- i
  macaque_gene <- macaque_anno[macaque_anno[,2] == human_gene,4]
  mouse_gene <- mouse_anno[mouse_anno[,2] == human_gene,4]
  
  #plot
  p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = human_gene,pt.size = 0.1,slot = 'data',
                    cols = as.character(ArchRPalettes$whitePurple[3:9]),order = FALSE) + 
    theme(aspect.ratio = 1) + NoLegend() + labs(title = human_gene) + NoAxes()
  
  p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = macaque_gene,pt.size = 0.1,slot = 'data',
                    cols = as.character(ArchRPalettes$whitePurple[3:9]),order = FALSE) + 
    theme(aspect.ratio = 1) + NoLegend() + labs(title = macaque_gene) + NoAxes()
  
  p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = mouse_gene,pt.size = 0.1,slot = 'data',
                    cols = as.character(ArchRPalettes$whitePurple[3:9]),order = FALSE) + 
    theme(aspect.ratio = 1) + NoLegend() + labs(title = mouse_gene) + NoAxes()
  
  #save pic
  char <- paste0('./res/step_94_fig_221226/three_species_',human_gene,'_featureplot.pdf')
  pdf(file = char,width = 18,height = 6)
  print(p1+p2+p3+plot_layout(ncol = 3))
  dev.off()
}