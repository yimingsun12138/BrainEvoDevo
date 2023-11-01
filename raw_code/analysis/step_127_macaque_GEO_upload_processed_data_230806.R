#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque GEO upload processed data                               ##
## Data: 2023.08.06                                                                ##
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
library(openxlsx)
library(readxl)
library(monocle)
library(destiny)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0613')

# macaque snRNA_seq -------------------------------------------------------

#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_RNA_Seurat_220803.rds')
macaque_RNA_Seurat$sample <- NA
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor == 'A50A',"sample"] <- 'E80_rep_1'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor == 'A82B',"sample"] <- 'E90_rep_2'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor == 'A84B',"sample"] <- 'E84_rep_1'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor == 'A84C',"sample"] <- 'E84_rep_2'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor == 'A68A',"sample"] <- 'E92_rep_1'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor == 'A68B',"sample"] <- 'E92_rep_2'
table(macaque_RNA_Seurat$donor,macaque_RNA_Seurat$sample)

#rename
temp <- gsub(pattern = '^A..._',replacement = '',x = colnames(macaque_RNA_Seurat),fixed = FALSE)
temp <- paste(macaque_RNA_Seurat$sample,temp,sep = '_')
macaque_RNA_Seurat <- RenameCells(object = macaque_RNA_Seurat,new.names = temp)

#export
macaque_raw_counts <- macaque_RNA_Seurat@assays$RNA@counts
macaque_meta_data <- macaque_RNA_Seurat@meta.data[,c('nCount_RNA','nFeature_RNA','batch','donor','sample','cell_type')]

write.table(x = macaque_raw_counts,file = '/home/sunym/temp/macaque_raw_counts_table.txt',sep = '\t')
write.table(x = macaque_meta_data,file = '/home/sunym/temp/macaque_meta_data.txt',sep = '\t')

system('md5sum /home/sunym/temp/macaque_raw_counts_table.txt',intern = TRUE)
#"b8d4ffb5a3117ca79d12eaf100d62566  /home/sunym/temp/macaque_raw_counts_table.txt"
system('md5sum /home/sunym/temp/macaque_meta_data.txt',intern = TRUE)
#"1cfb1da995559ffb8a4ba1e2b8668036  /home/sunym/temp/macaque_meta_data.txt"

# macaque scATAC_seq ------------------------------------------------------
macaque_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_ATAC_ArchR_220912/')
macaque_peakset <- getPeakSet(ArchRProj = macaque_ATAC_ArchR)
mcols(macaque_peakset) <- NULL
names(macaque_peakset) <- as.character(macaque_peakset)

macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_ATAC_ArchR,useMatrix = 'PeakMatrix')
macaque_peakset <- macaque_peak_matrix@rowRanges
mcols(macaque_peakset) <- NULL
names(macaque_peakset) <- NULL

macaque_peak_matrix <- macaque_peak_matrix@assays@data$PeakMatrix
table(macaque_ATAC_ArchR$donor)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A50A#',replacement = 'E80_rep_1_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A82A#',replacement = 'E90_rep_1_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A82B#',replacement = 'E90_rep_2_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A84B#',replacement = 'E84_rep_1_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A84C#',replacement = 'E84_rep_2_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A68A#',replacement = 'E92_rep_1_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A68B#',replacement = 'E92_rep_2_',x = colnames(macaque_peak_matrix),fixed = FALSE)

macaque_meta_data <- as.data.frame(macaque_ATAC_ArchR@cellColData)
rownames(macaque_meta_data) <- gsub(pattern = '^A50A#',replacement = 'E80_rep_1_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A82A#',replacement = 'E90_rep_1_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A82B#',replacement = 'E90_rep_2_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A84B#',replacement = 'E84_rep_1_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A84C#',replacement = 'E84_rep_2_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A68A#',replacement = 'E92_rep_1_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A68B#',replacement = 'E92_rep_2_',x = rownames(macaque_meta_data),fixed = FALSE)
macaque_meta_data <- macaque_meta_data[colnames(macaque_peak_matrix),c("donor","batch","TSSEnrichment","nFrags","ReadsInPeaks","cell_type")]

#export
export.bed(object = macaque_peakset,con = '/home/sunym/temp/macaque_peakset.bed')
writeMM(obj = macaque_peak_matrix,file = '/home/sunym/temp/macaque_peak_matrix.mtx')
write.table(x = macaque_meta_data,file = '/home/sunym/temp/macaque_meta_data.txt',sep = '\t')

system('md5sum /home/sunym/temp/macaque_peakset.bed',intern = TRUE)
#"a077def960fa8d93361e7a26454f4ffb  /home/sunym/temp/macaque_peakset.bed"
system('md5sum /home/sunym/temp/macaque_peak_matrix.mtx',intern = TRUE)
#"5521b3978d83a8ba79b6cd2d17f2ed46  /home/sunym/temp/macaque_peak_matrix.mtx"
system('md5sum /home/sunym/temp/macaque_meta_data.txt',intern = TRUE)
#"045fa8a090525a1e1e03b99c9663f06c  /home/sunym/temp/macaque_meta_data.txt"

# macaque scMultiome ------------------------------------------------------

## RNA ---------------------------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
macaque_multiome_Seurat$sample <- NA
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$donor == 'A50A',"sample"] <- 'E80_rep_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$donor == 'A50B',"sample"] <- 'E80_rep_2'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$donor == 'A82A',"sample"] <- 'E90_rep_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$donor == 'A82B',"sample"] <- 'E90_rep_2'
table(macaque_multiome_Seurat$donor,macaque_multiome_Seurat$sample)

#rename
temp <- gsub(pattern = '^A...#',replacement = '',x = colnames(macaque_multiome_Seurat),fixed = FALSE)
temp <- paste(macaque_multiome_Seurat$sample,temp,sep = '_')
macaque_multiome_Seurat <- RenameCells(object = macaque_multiome_Seurat,new.names = temp)

#export expression matrix
macaque_raw_counts <- macaque_multiome_Seurat@assays$RNA@counts
write.table(x = macaque_raw_counts,file = '/home/sunym/temp/macaque_raw_counts_table.txt',sep = '\t')

system('md5sum /home/sunym/temp/macaque_raw_counts_table.txt',intern = TRUE)
#"63c9f5d0af95293a14c0230545973995  /home/sunym/temp/macaque_raw_counts_table.txt"

## ATAC --------------------------------------------------------------------
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
macaque_peakset <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
mcols(macaque_peakset) <- NULL
names(macaque_peakset) <- as.character(macaque_peakset)

macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix')
macaque_peakset <- macaque_peak_matrix@rowRanges
mcols(macaque_peakset) <- NULL
names(macaque_peakset) <- NULL

macaque_peak_matrix <- macaque_peak_matrix@assays@data$PeakMatrix
table(macaque_multiome_ArchR$donor)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A50A#',replacement = 'E80_rep_1_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A50B#',replacement = 'E80_rep_2_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A82A#',replacement = 'E90_rep_1_',x = colnames(macaque_peak_matrix),fixed = FALSE)
colnames(macaque_peak_matrix) <- gsub(pattern = '^A82B#',replacement = 'E90_rep_2_',x = colnames(macaque_peak_matrix),fixed = FALSE)

macaque_meta_data <- as.data.frame(macaque_multiome_ArchR@cellColData)
rownames(macaque_meta_data) <- gsub(pattern = '^A50A#',replacement = 'E80_rep_1_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A50B#',replacement = 'E80_rep_2_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A82A#',replacement = 'E90_rep_1_',x = rownames(macaque_meta_data),fixed = FALSE)
rownames(macaque_meta_data) <- gsub(pattern = '^A82B#',replacement = 'E90_rep_2_',x = rownames(macaque_meta_data),fixed = FALSE)
macaque_meta_data <- macaque_meta_data[colnames(macaque_peak_matrix),c("donor","TSSEnrichment","nFrags","ReadsInPeaks","cell_type")]

temp_meta_data <- macaque_multiome_Seurat@meta.data
macaque_meta_data <- cbind(temp_meta_data[rownames(macaque_meta_data),c("nCount_RNA","nFeature_RNA")],macaque_meta_data)

#export
export.bed(object = macaque_peakset,con = '/home/sunym/temp/macaque_peakset.bed')
writeMM(obj = macaque_peak_matrix,file = '/home/sunym/temp/macaque_peak_matrix.mtx')
write.table(x = macaque_meta_data,file = '/home/sunym/temp/macaque_meta_data.txt',sep = '\t')

system('md5sum /home/sunym/temp/macaque_peakset.bed',intern = TRUE)
#"0622303bcb55581e50e28dbf4cba3424  /home/sunym/temp/macaque_peakset.bed"
system('md5sum /home/sunym/temp/macaque_peak_matrix.mtx',intern = TRUE)
#"74c00ad6412cae2b436d5ba156d2ef0d  /home/sunym/temp/macaque_peak_matrix.mtx"
system('md5sum /home/sunym/temp/macaque_meta_data.txt',intern = TRUE)
#"b2f61e1b627c567c0d34199822dde991  /home/sunym/temp/macaque_meta_data.txt"

# mouse scMultiome --------------------------------------------------------

## RNA ---------------------------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
mouse_multiome_Seurat$sample <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor == 'E145_1',"sample"] <- 'E14_rep_1'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor == 'E145_2',"sample"] <- 'E14_rep_2'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor == 'E155_1',"sample"] <- 'E15_rep_1'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor == 'E155_2',"sample"] <- 'E15_rep_2'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor == 'E155_3',"sample"] <- 'E15_rep_3'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor == 'E165_2',"sample"] <- 'E16_rep_1'
table(mouse_multiome_Seurat$sample)

#rename
temp <- gsub(pattern = '^E.....#',replacement = '',x = colnames(mouse_multiome_Seurat),fixed = FALSE)
temp <- paste(mouse_multiome_Seurat$sample,temp,sep = '_')
mouse_multiome_Seurat <- RenameCells(object = mouse_multiome_Seurat,new.names = temp)

#export expression matrix
mouse_raw_counts <- mouse_multiome_Seurat@assays$RNA@counts
write.table(x = mouse_raw_counts,file = '/home/sunym/temp/mouse_raw_counts_table.txt',sep = '\t')

system('md5sum /home/sunym/temp/mouse_raw_counts_table.txt',intern = TRUE)
# "45371e292dae45c3a2d30525a01cca3d  /home/sunym/temp/mouse_raw_counts_table.txt"

## ATAC --------------------------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009')
mouse_multiome_ArchR$cell_type <- mouse_multiome_ArchR$Gex_cell_type
mouse_peakset <- getPeakSet(ArchRProj = mouse_multiome_ArchR)
mcols(mouse_peakset) <- NULL
names(mouse_peakset) <- as.character(mouse_peakset)

mouse_peak_matrix <- getMatrixFromProject(ArchRProj = mouse_multiome_ArchR,useMatrix = 'PeakMatrix')
mouse_peakset <- mouse_peak_matrix@rowRanges
mcols(mouse_peakset) <- NULL
names(mouse_peakset) <- NULL

mouse_peak_matrix <- mouse_peak_matrix@assays@data$PeakMatrix
table(mouse_multiome_ArchR$donor)
colnames(mouse_peak_matrix) <- gsub(pattern = '^E145_1#',replacement = 'E14_rep_1_',x = colnames(mouse_peak_matrix),fixed = FALSE)
colnames(mouse_peak_matrix) <- gsub(pattern = '^E145_2#',replacement = 'E14_rep_2_',x = colnames(mouse_peak_matrix),fixed = FALSE)
colnames(mouse_peak_matrix) <- gsub(pattern = '^E155_1#',replacement = 'E15_rep_1_',x = colnames(mouse_peak_matrix),fixed = FALSE)
colnames(mouse_peak_matrix) <- gsub(pattern = '^E155_2#',replacement = 'E15_rep_2_',x = colnames(mouse_peak_matrix),fixed = FALSE)
colnames(mouse_peak_matrix) <- gsub(pattern = '^E155_3#',replacement = 'E15_rep_3_',x = colnames(mouse_peak_matrix),fixed = FALSE)
colnames(mouse_peak_matrix) <- gsub(pattern = '^E165_2#',replacement = 'E16_rep_1_',x = colnames(mouse_peak_matrix),fixed = FALSE)

mouse_meta_data <- as.data.frame(mouse_multiome_ArchR@cellColData)
rownames(mouse_meta_data) <- gsub(pattern = '^E145_1#',replacement = 'E14_rep_1_',x = rownames(mouse_meta_data),fixed = FALSE)
rownames(mouse_meta_data) <- gsub(pattern = '^E145_2#',replacement = 'E14_rep_2_',x = rownames(mouse_meta_data),fixed = FALSE)
rownames(mouse_meta_data) <- gsub(pattern = '^E155_1#',replacement = 'E15_rep_1_',x = rownames(mouse_meta_data),fixed = FALSE)
rownames(mouse_meta_data) <- gsub(pattern = '^E155_2#',replacement = 'E15_rep_2_',x = rownames(mouse_meta_data),fixed = FALSE)
rownames(mouse_meta_data) <- gsub(pattern = '^E155_3#',replacement = 'E15_rep_3_',x = rownames(mouse_meta_data),fixed = FALSE)
rownames(mouse_meta_data) <- gsub(pattern = '^E165_2#',replacement = 'E16_rep_1_',x = rownames(mouse_meta_data),fixed = FALSE)
mouse_meta_data <- mouse_meta_data[colnames(mouse_peak_matrix),c("donor","TSSEnrichment","nFrags","ReadsInPeaks","cell_type")]

temp_meta_data <- mouse_multiome_Seurat@meta.data
mouse_meta_data <- cbind(temp_meta_data[rownames(mouse_meta_data),c("nCount_RNA","nFeature_RNA")],mouse_meta_data)

#export
export.bed(object = mouse_peakset,con = '/home/sunym/temp/mouse_peakset.bed')
writeMM(obj = mouse_peak_matrix,file = '/home/sunym/temp/mouse_peak_matrix.mtx')
write.table(x = mouse_meta_data,file = '/home/sunym/temp/mouse_meta_data.txt',sep = '\t')

system('md5sum /home/sunym/temp/mouse_peakset.bed',intern = TRUE)
# "395f63b0ae6ca8de6d8df045a06eb23c  /home/sunym/temp/mouse_peakset.bed"
system('md5sum /home/sunym/temp/mouse_peak_matrix.mtx',intern = TRUE)
# "8733def8049c14d553b0b399b51cb248  /home/sunym/temp/mouse_peak_matrix.mtx"
system('md5sum /home/sunym/temp/mouse_meta_data.txt',intern = TRUE)
# "675a4f10a77f1a9c080e558f3929dbb8  /home/sunym/temp/mouse_meta_data.txt"