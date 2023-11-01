#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: RG peak comparison across species                               ##
## Data: 2022.05.11                                                                ##
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
setwd('/data/User/sunym/project/Brain/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(parallel)
library(Seurat)
library(ArchR)
library(rtracklayer)
library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(ComplexHeatmap)
library(parallel)
library(patchwork)
library(ggpubr)
library(ggvenn)
library(circlize)
library(ChIPpeakAnno)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# target gene of DAP ------------------------------------------------------
# # the first step of analysis is done on 202.205.131.32, cause ArchR is fucking shit.
# Brain_ATAC_peakset <- readRDS(file = './res/step_27_fig_220510/Brain_ATAC_peakset.rds')
# macaque_multiome_ArchR
# Greenleaf_ATAC_ArchR
# mouse_ATAC_ArchR
# 
# #add peakset
# #human
# temp <- Greenleaf_ATAC_ArchR@cellColData
# ArrowFiles
# file.exists(ArrowFiles)
# Greenleaf_ATAC_ArchR <- ArchRProject(
#   ArrowFiles = ArrowFiles, 
#   outputDirectory = '/data/User/sunym/temp/Greenleaf_ATAC_ArchR',
#   copyArrows = TRUE,
#   geneAnnotation = getGeneAnnotation(ArchRProj = Greenleaf_ATAC_ArchR),
#   genomeAnnotation = getGenomeAnnotation(ArchRProj = Greenleaf_ATAC_ArchR)
# )
# Greenleaf_ATAC_ArchR <- Greenleaf_ATAC_ArchR[rownames(temp)]
# Greenleaf_ATAC_ArchR@cellColData <- temp
# 
# Greenleaf_ATAC_ArchR <- addIterativeLSI(ArchRProj = Greenleaf_ATAC_ArchR,
#                                         useMatrix = "TileMatrix", 
#                                         name = "IterativeLSI", 
#                                         iterations = 3, 
#                                         clusterParams = list(
#                                           resolution = c(0.6), 
#                                           sampleCells = 10000, 
#                                           n.start = 10,
#                                           maxClusters = 30
#                                         ), 
#                                         varFeatures = 20000, 
#                                         totalFeatures = 500000, 
#                                         dimsToUse = 1:30, 
#                                         corCutOff = 0.5, 
#                                         force = TRUE)
# 
# Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = Brain_ATAC_peakset$human,
#                                    genomeAnnotation = getGenomeAnnotation(Greenleaf_ATAC_ArchR),force = TRUE)
# Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR,verbose = TRUE,force = TRUE)
# 
# #macaque
# temp <- macaque_multiome_ArchR@cellColData
# ArrowFiles
# file.exists(ArrowFiles)
# macaque_multiome_ArchR <- ArchRProject(
#   ArrowFiles = ArrowFiles, 
#   outputDirectory = '/data/User/sunym/temp/macaque_multiome_ArchR',
#   copyArrows = TRUE,
#   geneAnnotation = getGeneAnnotation(ArchRProj = macaque_multiome_ArchR),
#   genomeAnnotation = getGenomeAnnotation(ArchRProj = macaque_multiome_ArchR)
# )
# macaque_multiome_ArchR <- macaque_multiome_ArchR[rownames(temp)]
# macaque_multiome_ArchR@cellColData <- temp
# 
# macaque_multiome_ArchR <- addIterativeLSI(ArchRProj = macaque_multiome_ArchR,
#                                           useMatrix = "TileMatrix", 
#                                           name = "IterativeLSI", 
#                                           iterations = 3, 
#                                           clusterParams = list(
#                                             resolution = c(0.6), 
#                                             sampleCells = 10000, 
#                                             n.start = 10,
#                                             maxClusters = 30
#                                           ), 
#                                           varFeatures = 20000, 
#                                           totalFeatures = 500000, 
#                                           dimsToUse = 1:30, 
#                                           corCutOff = 0.5, 
#                                           force = TRUE)
# 
# macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = Brain_ATAC_peakset$macaque,
#                                      genomeAnnotation = getGenomeAnnotation(macaque_multiome_ArchR),force = TRUE)
# macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)
# 
# #mouse
# temp <- mouse_ATAC_ArchR@cellColData
# ArrowFiles
# file.exists(ArrowFiles)
# mouse_ATAC_ArchR <- ArchRProject(
#   ArrowFiles = ArrowFiles, 
#   outputDirectory = '/data/User/sunym/temp/mouse_ATAC_ArchR',
#   copyArrows = TRUE,
#   geneAnnotation = getGeneAnnotation(ArchRProj = mouse_ATAC_ArchR),
#   genomeAnnotation = getGenomeAnnotation(ArchRProj = mouse_ATAC_ArchR)
# )
# mouse_ATAC_ArchR <- mouse_ATAC_ArchR[rownames(temp)]
# mouse_ATAC_ArchR@cellColData <- temp
# 
# mouse_ATAC_ArchR <- addIterativeLSI(ArchRProj = mouse_ATAC_ArchR,
#                                     useMatrix = "TileMatrix", 
#                                     name = "IterativeLSI", 
#                                     iterations = 3, 
#                                     clusterParams = list(
#                                       resolution = c(0.6), 
#                                       sampleCells = 10000, 
#                                       n.start = 10,
#                                       maxClusters = 30
#                                     ), 
#                                     varFeatures = 20000, 
#                                     totalFeatures = 500000, 
#                                     dimsToUse = 1:30, 
#                                     corCutOff = 0.5, 
#                                     force = TRUE)
# 
# mouse_ATAC_ArchR <- addPeakSet(ArchRProj = mouse_ATAC_ArchR,peakSet = Brain_ATAC_peakset$mouse,
#                                genomeAnnotation = getGenomeAnnotation(mouse_ATAC_ArchR),force = TRUE)
# mouse_ATAC_ArchR <- addPeakMatrix(ArchRProj = mouse_ATAC_ArchR,verbose = TRUE,force = TRUE)
# 
# #calculate P2G
# #human
# Greenleaf_ATAC_ArchR <- addPeak2GeneLinks(
#   ArchRProj = Greenleaf_ATAC_ArchR,
#   reducedDims = "IterativeLSI",
#   useMatrix = 'GeneIntegrationMatrix'
# )
# p2g <- getPeak2GeneLinks(
#   ArchRProj = Greenleaf_ATAC_ArchR,
#   corCutOff = 0.45,
#   resolution = 1,
#   returnLoops = FALSE
# )
# saveRDS(p2g,file = './res/step_28_fig_220511/human_p2g.rds')
# 
# #macaque
# macaque_multiome_ArchR <- addPeak2GeneLinks(
#   ArchRProj = macaque_multiome_ArchR,
#   reducedDims = "IterativeLSI",
#   useMatrix = 'GeneExpressionMatrix'
# )
# p2g <- getPeak2GeneLinks(
#   ArchRProj = macaque_multiome_ArchR,
#   corCutOff = 0.45,
#   resolution = 1,
#   returnLoops = FALSE
# )
# saveRDS(p2g,file = './res/step_28_fig_220511/macaque_p2g.rds')
# 
# #mouse
# mouse_ATAC_ArchR <- addPeak2GeneLinks(
#   ArchRProj = mouse_ATAC_ArchR,
#   reducedDims = "IterativeLSI",
#   useMatrix = 'GeneIntegrationMatrix'
# )
# p2g <- getPeak2GeneLinks(
#   ArchRProj = mouse_ATAC_ArchR,
#   corCutOff = 0.45,
#   resolution = 1,
#   returnLoops = FALSE
# )
# saveRDS(p2g,file = './res/step_28_fig_220511/mouse_p2g.rds')

# count target of the p2g -------------------------------------------------
#load data
human_p2g <- readRDS(file = './res/step_28_fig_220511/human_p2g.rds')
macaque_p2g <- readRDS(file = './res/step_28_fig_220511/macaque_p2g.rds')
mouse_p2g <- readRDS(file = './res/step_28_fig_220511/mouse_p2g.rds')

#load peakset
Brain_ATAC_peakset <- readRDS(file = './res/step_27_fig_220510/Brain_ATAC_peakset.rds')
group_list <- readRDS(file = './res/step_27_fig_220510/Brain_peak_matrix_RG_gene_list_kmeans_group.rds')

#create function
count_target <- function(peak_list,p2g){
  #select peak anchor
  p2g@metadata$peakSet$idx <- 1:length(p2g@metadata$peakSet)
  p2g@metadata$peakSet$name <- paste(p2g@metadata$peakSet@seqnames,as.character(p2g@metadata$peakSet@ranges),sep = '-')
  temp <- p2g@metadata$peakSet[p2g@metadata$peakSet$name %in% peak_list]
  #count gene
  gene_list <- p2g[p2g$idxATAC %in% temp$idx,]
  gene_list <- p2g@metadata$geneSet$name[gene_list$idxRNA]
  #create data.frame
  temp <- data.frame(table(gene_list))
  colnames(temp) <- c('gene','counts')
  temp$counts <- temp$counts
  rownames(temp) <- temp$gene
  temp$p <- NA
  for (i in rownames(temp)) {
    temp[i,'p'] <- dhyper(x = sum(gene_list == i),
                          m = sum(p2g@metadata$geneSet$name[p2g$idxRNA] == i),
                          n = sum(p2g@metadata$geneSet$name[p2g$idxRNA] != i),
                          k = length(gene_list))
  }
  temp$p <- -log10(temp$p)
  #return
  return(temp)
}
## km3 -- human specific ---------------------------------------------------
peak_list <- group_list[group_list == 3]
peak_list <- names(peak_list)
temp <- count_target(peak_list = peak_list,p2g = human_p2g)
temp$group <- 'km3'
target_matrix <- temp

peak_list <- group_list[group_list != 3]
peak_list <- sample(peak_list,size = sum(group_list == 3),replace = FALSE)
peak_list <- names(peak_list)
temp <- count_target(peak_list = peak_list,p2g = human_p2g)
temp$group <- 'background'
target_matrix <- rbind(target_matrix,temp)

target_matrix$label <- target_matrix$gene
target_matrix[target_matrix$group == 'background',"label"] <- NA
target_matrix[target_matrix$p < 3,"label"] <- NA

ggplot(data = target_matrix,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','km3' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')

## km4 -- macaque specific ---------------------------------------------------
peak_list <- group_list[group_list == 4]
peak_list <- names(peak_list)
peak_list <- Brain_ATAC_peakset$macaque[Brain_ATAC_peakset$macaque$name %in% peak_list]
peak_list <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')
temp <- count_target(peak_list = peak_list,p2g = macaque_p2g)
temp$group <- 'km4'
target_matrix <- temp

peak_list <- group_list[group_list != 4]
peak_list <- sample(peak_list,size = sum(group_list == 4),replace = FALSE)
peak_list <- names(peak_list)
peak_list <- Brain_ATAC_peakset$macaque[Brain_ATAC_peakset$macaque$name %in% peak_list]
peak_list <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')
temp <- count_target(peak_list = peak_list,p2g = macaque_p2g)
temp$group <- 'background'
target_matrix <- rbind(target_matrix,temp)

target_matrix$label <- target_matrix$gene
target_matrix[target_matrix$group == 'background',"label"] <- NA
target_matrix[target_matrix$p < 3,"label"] <- NA

ggplot(data = target_matrix,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','km4' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')

## km1 -- mouse specific ---------------------------------------------------
peak_list <- group_list[group_list == 1]
peak_list <- names(peak_list)
peak_list <- Brain_ATAC_peakset$mouse[Brain_ATAC_peakset$mouse$name %in% peak_list]
peak_list <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')
temp <- count_target(peak_list = peak_list,p2g = mouse_p2g)
temp$group <- 'km1'
target_matrix <- temp

peak_list <- group_list[group_list != 1]
peak_list <- sample(peak_list,size = sum(group_list == 1),replace = FALSE)
peak_list <- names(peak_list)
peak_list <- Brain_ATAC_peakset$mouse[Brain_ATAC_peakset$mouse$name %in% peak_list]
peak_list <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')
temp <- count_target(peak_list = peak_list,p2g = mouse_p2g)
temp$group <- 'background'
target_matrix <- rbind(target_matrix,temp)

target_matrix$label <- target_matrix$gene
target_matrix[target_matrix$group == 'background',"label"] <- NA
target_matrix[target_matrix$p < 3,"label"] <- NA

ggplot(data = target_matrix,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','km1' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')

# focus on km3 target -- human specific -----------------------------------
peak_list <- group_list[group_list == 3]
peak_list <- names(peak_list)
temp <- count_target(peak_list = peak_list,p2g = human_p2g)
temp$group <- 'km3'
target_matrix <- temp

peak_list <- group_list[group_list != 3]
peak_list <- sample(peak_list,size = sum(group_list == 3),replace = FALSE)
peak_list <- names(peak_list)
temp <- count_target(peak_list = peak_list,p2g = human_p2g)
temp$group <- 'background'
target_matrix <- rbind(target_matrix,temp)

target_matrix$label <- target_matrix$gene
target_matrix[target_matrix$group == 'background',"label"] <- NA
target_matrix[target_matrix$p < 3,"label"] <- NA

ggplot(data = target_matrix,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','km3' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')

# export cell type level coverage -----------------------------------------

# #human
# temp
# ArrowFiles
# file.exists(ArrowFiles)
# Greenleaf_ATAC_ArchR <- ArchRProject(
#   ArrowFiles = ArrowFiles,
#   outputDirectory = './Greenleaf_ATAC_ArchR',
#   copyArrows = TRUE,
#   geneAnnotation = getGeneAnnotation(ArchRProj = temp),
#   genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
# )
# Greenleaf_ATAC_ArchR <- Greenleaf_ATAC_ArchR[rownames(temp@cellColData)]
# Greenleaf_ATAC_ArchR@cellColData <- temp@cellColData
# 
# table(Greenleaf_ATAC_ArchR$cell_type)
# coverage_file <- getGroupBW(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',normMethod = 'nFrags',tileSize = 25,maxCells = 100000,verbose = TRUE)
# 
# #macaque
# temp
# ArrowFiles
# file.exists(ArrowFiles)
# macaque_multiome_ArchR <- ArchRProject(
#   ArrowFiles = ArrowFiles,
#   outputDirectory = './macaque_multiome_ArchR',
#   copyArrows = TRUE,
#   geneAnnotation = getGeneAnnotation(ArchRProj = temp),
#   genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
# )
# macaque_multiome_ArchR <- macaque_multiome_ArchR[rownames(temp@cellColData)]
# macaque_multiome_ArchR@cellColData <- temp@cellColData
# 
# table(macaque_multiome_ArchR$cell_type)
# coverage_file <- getGroupBW(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',normMethod = 'nFrags',tileSize = 25,maxCells = 100000,verbose = TRUE)
# 
# #mouse
# temp
# ArrowFiles
# file.exists(ArrowFiles)
# mouse_ATAC_ArchR <- ArchRProject(
#   ArrowFiles = ArrowFiles,
#   outputDirectory = './mouse_ATAC_ArchR',
#   copyArrows = TRUE,
#   geneAnnotation = getGeneAnnotation(ArchRProj = temp),
#   genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
# )
# mouse_ATAC_ArchR <- mouse_ATAC_ArchR[rownames(temp@cellColData)]
# mouse_ATAC_ArchR@cellColData <- temp@cellColData
# 
# table(mouse_ATAC_ArchR$cell_type)
# coverage_file <- getGroupBW(ArchRProj = mouse_ATAC_ArchR,groupBy = 'cell_type',normMethod = 'nFrags',tileSize = 25,maxCells = 100000,verbose = TRUE)
