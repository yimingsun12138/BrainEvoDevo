#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: create species peak matrix normalized by cell type ReadsInPeaks ##
## Data: 2022.06.06                                                                ##
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
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnrichedHeatmap)
library(circlize)
library(scales)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(ggpointdensity)
library(networkD3)
library(htmlwidgets)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# create species peak matrix ----------------------------------------------
#all run on 202.205.131.32, ArchR is fucking shit.

## human peak matrix -------------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

#create new ArchR object
temp <- readRDS(file = '/data/User/sunym/project/Brain/data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
ArrowFiles <- list.files('/data/User/sunym/project/Brain/data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/arrow_file/')
ArrowFiles <- paste('/data/User/sunym/project/Brain/data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/arrow_file',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

Greenleaf_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/temp/Greenleaf_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

Greenleaf_ATAC_ArchR <- Greenleaf_ATAC_ArchR[rownames(Greenleaf_ATAC_ArchR@cellColData) %in% rownames(temp@cellColData)]

#add meta data
Greenleaf_ATAC_ArchR$cell_type <- temp@cellColData[rownames(Greenleaf_ATAC_ArchR@cellColData),"cell_type"]
Greenleaf_ATAC_ArchR$Reads_In_Cell_Type_Peaks <- NA
cell_type_list <- unique(Greenleaf_ATAC_ArchR$cell_type)
for (i in cell_type_list) {
  #get cell type peakset
  char <- list.files('/data/User/sunym/project/Brain/data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
  ii <- sub(pattern = '-',replacement = '.',x = i,fixed = TRUE)
  ii <- sub(pattern = '/',replacement = '.',x = ii,fixed = TRUE)
  char <- char[grep(pattern = paste0('^',ii,'-'),x = char,fixed = FALSE)]
  char <- paste('/data/User/sunym/project/Brain/data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',char,sep = '/')
  temp_peakset <- readRDS(file = char)
  #add peakset
  Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = temp_peakset,
                                     genomeAnnotation = getGenomeAnnotation(Greenleaf_ATAC_ArchR),force = TRUE)
  Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR,verbose = TRUE,force = TRUE)
  cell_list <- rownames(Greenleaf_ATAC_ArchR@cellColData)[Greenleaf_ATAC_ArchR$cell_type == i]
  Greenleaf_ATAC_ArchR@cellColData[cell_list,'Reads_In_Cell_Type_Peaks'] <- Greenleaf_ATAC_ArchR@cellColData[cell_list,'ReadsInPeaks']
  print(paste(i,'done!',sep = ' '))
}

#add species orthologous peakset
Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = Brain_ATAC_peakset$human,
                                   genomeAnnotation = getGenomeAnnotation(Greenleaf_ATAC_ArchR),force = TRUE)
Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR,verbose = TRUE,force = TRUE)
human_peak_matrix <- getMatrixFromProject(ArchRProj = Greenleaf_ATAC_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

#create human peak matrix
meta_data <- as.data.frame(human_peak_matrix@colData)
temp <- human_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- paste(human_peak_matrix@rowRanges@seqnames,as.character(human_peak_matrix@rowRanges@ranges),sep = '-')
human_peak_matrix <- CreateSeuratObject(counts = temp,project = 'human',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

human_peak_matrix$pseudo_bulk <- NA
for (i in unique(human_peak_matrix$cell_type)) {
  j <- sum(human_peak_matrix$cell_type == i)
  pseudo_bulk <- sample(rep(c('rep1','rep2','rep3'),times = c(round(j/3),round(j/3),(j-2*round(j/3)))))
  pseudo_bulk <- paste(i,pseudo_bulk,sep = '#')
  human_peak_matrix@meta.data[human_peak_matrix$cell_type == i,'pseudo_bulk'] <- pseudo_bulk
}

#save data
saveRDS(human_peak_matrix,file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')

## macaque peak matrix -----------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

#create new ArchR object
temp <- readRDS(file = '/data/User/sunym/project/Brain/ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
ArrowFiles <- list.files('/data/User/sunym/project/Brain/ArchR/arrow_file/macaque_multiome_data/')
ArrowFiles <- paste('/data/User/sunym/project/Brain/ArchR/arrow_file/macaque_multiome_data/',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/temp/macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

macaque_multiome_ArchR <- macaque_multiome_ArchR[rownames(macaque_multiome_ArchR@cellColData) %in% rownames(temp@cellColData)]

#add meta data
macaque_multiome_ArchR$cell_type <- temp@cellColData[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]
macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks <- NA
cell_type_list <- unique(macaque_multiome_ArchR$cell_type)
for (i in cell_type_list) {
  #get cell type peakset
  char <- list.files('/data/User/sunym/project/Brain/ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/')
  ii <- sub(pattern = '-',replacement = '.',x = i,fixed = TRUE)
  ii <- sub(pattern = '/',replacement = '.',x = ii,fixed = TRUE)
  char <- char[grep(pattern = paste0('^',ii,'-'),x = char,fixed = FALSE)]
  char <- paste('/data/User/sunym/project/Brain/ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls',char,sep = '/')
  temp_peakset <- readRDS(file = char)
  #add peakset
  macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = temp_peakset,
                                       genomeAnnotation = getGenomeAnnotation(macaque_multiome_ArchR),force = TRUE)
  macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)
  cell_list <- rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$cell_type == i]
  macaque_multiome_ArchR@cellColData[cell_list,'Reads_In_Cell_Type_Peaks'] <- macaque_multiome_ArchR@cellColData[cell_list,'ReadsInPeaks']
  print(paste(i,'done!',sep = ' '))
}

#add species orthologous peakset
macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = Brain_ATAC_peakset$macaque,
                                     genomeAnnotation = getGenomeAnnotation(macaque_multiome_ArchR),force = TRUE)
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)
macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

#create macaque peak matrix
meta_data <- as.data.frame(macaque_peak_matrix@colData)
temp <- macaque_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- paste(macaque_peak_matrix@rowRanges@seqnames,as.character(macaque_peak_matrix@rowRanges@ranges),sep = '-')
macaque_peak_matrix <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

macaque_peak_matrix$pseudo_bulk <- NA
for (i in unique(macaque_peak_matrix$cell_type)) {
  j <- sum(macaque_peak_matrix$cell_type == i)
  pseudo_bulk <- sample(rep(c('rep1','rep2','rep3'),times = c(round(j/3),round(j/3),(j-2*round(j/3)))))
  pseudo_bulk <- paste(i,pseudo_bulk,sep = '#')
  macaque_peak_matrix@meta.data[macaque_peak_matrix$cell_type == i,'pseudo_bulk'] <- pseudo_bulk
}

#save data
saveRDS(macaque_peak_matrix,file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')

## mouse peak matrix -------------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

#create new ArchR object
temp <- readRDS(file = '/data/User/sunym/project/Brain/data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')
ArrowFiles <- list.files('/data/User/sunym/project/Brain/data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/arrow_file/')
ArrowFiles <- paste('/data/User/sunym/project/Brain/data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/arrow_file',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

mouse_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/temp/mouse_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

mouse_ATAC_ArchR <- mouse_ATAC_ArchR[rownames(mouse_ATAC_ArchR@cellColData) %in% rownames(temp@cellColData)]

#add meta data
mouse_ATAC_ArchR$cell_type <- temp@cellColData[rownames(mouse_ATAC_ArchR@cellColData),"cell_type"]
mouse_ATAC_ArchR$Reads_In_Cell_Type_Peaks <- NA
cell_type_list <- unique(mouse_ATAC_ArchR$cell_type)
for (i in cell_type_list) {
  #get cell type peakset
  char <- list.files('/data/User/sunym/project/Brain/data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/')
  ii <- sub(pattern = '-',replacement = '.',x = i,fixed = TRUE)
  ii <- sub(pattern = '/',replacement = '.',x = ii,fixed = TRUE)
  char <- char[grep(pattern = paste0('^',ii,'-'),x = char,fixed = FALSE)]
  char <- paste('/data/User/sunym/project/Brain/data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls',char,sep = '/')
  temp_peakset <- readRDS(file = char)
  #add peakset
  mouse_ATAC_ArchR <- addPeakSet(ArchRProj = mouse_ATAC_ArchR,peakSet = temp_peakset,
                                 genomeAnnotation = getGenomeAnnotation(mouse_ATAC_ArchR),force = TRUE)
  mouse_ATAC_ArchR <- addPeakMatrix(ArchRProj = mouse_ATAC_ArchR,verbose = TRUE,force = TRUE)
  cell_list <- rownames(mouse_ATAC_ArchR@cellColData)[mouse_ATAC_ArchR$cell_type == i]
  mouse_ATAC_ArchR@cellColData[cell_list,'Reads_In_Cell_Type_Peaks'] <- mouse_ATAC_ArchR@cellColData[cell_list,'ReadsInPeaks']
  print(paste(i,'done!',sep = ' '))
}

#add species orthologous peakset
mouse_ATAC_ArchR <- addPeakSet(ArchRProj = mouse_ATAC_ArchR,peakSet = Brain_ATAC_peakset$mouse,
                               genomeAnnotation = getGenomeAnnotation(mouse_ATAC_ArchR),force = TRUE)
mouse_ATAC_ArchR <- addPeakMatrix(ArchRProj = mouse_ATAC_ArchR,verbose = TRUE,force = TRUE)
mouse_peak_matrix <- getMatrixFromProject(ArchRProj = mouse_ATAC_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

#create mouse peak matrix
meta_data <- as.data.frame(mouse_peak_matrix@colData)
temp <- mouse_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- paste(mouse_peak_matrix@rowRanges@seqnames,as.character(mouse_peak_matrix@rowRanges@ranges),sep = '-')
mouse_peak_matrix <- CreateSeuratObject(counts = temp,project = 'mouse',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

mouse_peak_matrix$pseudo_bulk <- NA
for (i in unique(mouse_peak_matrix$cell_type)) {
  j <- sum(mouse_peak_matrix$cell_type == i)
  pseudo_bulk <- sample(rep(c('rep1','rep2','rep3'),times = c(round(j/3),round(j/3),(j-2*round(j/3)))))
  pseudo_bulk <- paste(i,pseudo_bulk,sep = '#')
  mouse_peak_matrix@meta.data[mouse_peak_matrix$cell_type == i,'pseudo_bulk'] <- pseudo_bulk
}

#save data
saveRDS(mouse_peak_matrix,file = './res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

# investigate reads in cell type peaks and reads in peaks -----------------
#load data
human_peak_matrix <- readRDS(file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')
macaque_peak_matrix <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
mouse_peak_matrix <- readRDS(file = './res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

#human
pdf(file = './res/step_35_fig_220606/human_reads_in_peaks_dotplot.pdf',width = 6,height = 6)
ggplot(data = human_peak_matrix@meta.data,aes(x = Reads_In_Cell_Type_Peaks,y = ReadsInPeaks)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'blue',size = 0.5) + 
  geom_abline(intercept = 0,slope = 1,color = 'red',size = 0.5,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = 'Human Reads in Peaks')
dev.off()

#macaque
pdf(file = './res/step_35_fig_220606/macaque_reads_in_peaks_dotplot.pdf',width = 6,height = 6)
ggplot(data = macaque_peak_matrix@meta.data,aes(x = Reads_In_Cell_Type_Peaks,y = ReadsInPeaks)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'blue',size = 0.5) + 
  geom_abline(intercept = 0,slope = 1,color = 'red',size = 0.5,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = 'Macaque Reads in Peaks')
dev.off()

#mouse
pdf(file = './res/step_35_fig_220606/mouse_reads_in_peaks_dotplot.pdf',width = 6,height = 6)
ggplot(data = macaque_peak_matrix@meta.data,aes(x = Reads_In_Cell_Type_Peaks,y = ReadsInPeaks)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'blue',size = 0.5) + 
  geom_abline(intercept = 0,slope = 1,color = 'red',size = 0.5,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = 'Mouse Reads in Peaks')
dev.off()

#I decide to use cell type specific peakset calculated ReadsInPeaks to generate bw file.

# create Brain peak matrix ------------------------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')
human_peak_Seurat <- readRDS(file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')
macaque_peak_Seurat <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
mouse_peak_Seurat <- readRDS(file = './res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

#species peak matrix
human_peak_matrix <- Seurat::AggregateExpression(object = human_peak_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'pseudo_bulk',slot = 'counts',verbose = TRUE)
human_peak_matrix <- human_peak_matrix$RNA
colnames(human_peak_matrix) <- paste('human',colnames(human_peak_matrix),sep = '#')

macaque_peak_matrix <- Seurat::AggregateExpression(object = macaque_peak_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'pseudo_bulk',slot = 'counts',verbose = TRUE)
macaque_peak_matrix <- macaque_peak_matrix$RNA
colnames(macaque_peak_matrix) <- paste('macaque',colnames(macaque_peak_matrix),sep = '#')

mouse_peak_matrix <- Seurat::AggregateExpression(object = mouse_peak_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'pseudo_bulk',slot = 'counts',verbose = TRUE)
mouse_peak_matrix <- mouse_peak_matrix$RNA
colnames(mouse_peak_matrix) <- paste('mouse',colnames(mouse_peak_matrix),sep = '#')

#rename features
names(Brain_ATAC_peakset$macaque) <- paste(Brain_ATAC_peakset$macaque@seqnames,as.character(Brain_ATAC_peakset$macaque@ranges),sep = '-')
names(Brain_ATAC_peakset$mouse) <- paste(Brain_ATAC_peakset$mouse@seqnames,as.character(Brain_ATAC_peakset$mouse@ranges),sep = '-')
rownames(macaque_peak_matrix) <- Brain_ATAC_peakset$macaque[rownames(macaque_peak_matrix)]$name
rownames(mouse_peak_matrix) <- Brain_ATAC_peakset$mouse[rownames(mouse_peak_matrix)]$name

#create Brain peak matrix
Brain_peak_matrix <- cbind(human_peak_matrix,macaque_peak_matrix[rownames(human_peak_matrix),],mouse_peak_matrix[rownames(human_peak_matrix),])
meta_data <- data.frame(sample=colnames(Brain_peak_matrix))
meta_data$species <- unlist(base::lapply(X = meta_data$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][1])
}))
meta_data$cell_type <- unlist(base::lapply(X = meta_data$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))
rownames(meta_data) <- meta_data$sample

meta_data$nFrags <- unlist(base::lapply(X = rownames(meta_data),FUN = function(x){
  temp <- strsplit(x = x,split = '#')[[1]]
  if(temp[1] == 'human'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(human_peak_Seurat)[human_peak_Seurat$pseudo_bulk == temp]
    num <- sum(human_peak_Seurat@meta.data[temp,'nFrags'])
    return(num)
  }else if(temp[1] == 'macaque'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(macaque_peak_Seurat)[macaque_peak_Seurat$pseudo_bulk == temp]
    num <- sum(macaque_peak_Seurat@meta.data[temp,'nFrags'])
    return(num)
  }else if(temp[1] == 'mouse'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(mouse_peak_Seurat)[mouse_peak_Seurat$pseudo_bulk == temp]
    num <- sum(mouse_peak_Seurat@meta.data[temp,'nFrags'])
    return(num)
  }else{
    stop('error!')
  }
}))

meta_data$ReadsInTSS <- unlist(base::lapply(X = rownames(meta_data),FUN = function(x){
  temp <- strsplit(x = x,split = '#')[[1]]
  if(temp[1] == 'human'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(human_peak_Seurat)[human_peak_Seurat$pseudo_bulk == temp]
    num <- sum(human_peak_Seurat@meta.data[temp,'ReadsInTSS'])
    return(num)
  }else if(temp[1] == 'macaque'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(macaque_peak_Seurat)[macaque_peak_Seurat$pseudo_bulk == temp]
    num <- sum(macaque_peak_Seurat@meta.data[temp,'ReadsInTSS'])
    return(num)
  }else if(temp[1] == 'mouse'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(mouse_peak_Seurat)[mouse_peak_Seurat$pseudo_bulk == temp]
    num <- sum(mouse_peak_Seurat@meta.data[temp,'ReadsInTSS'])
    return(num)
  }else{
    stop('error!')
  }
}))

meta_data$ReadsInPeaks <- unlist(base::lapply(X = rownames(meta_data),FUN = function(x){
  temp <- strsplit(x = x,split = '#')[[1]]
  if(temp[1] == 'human'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(human_peak_Seurat)[human_peak_Seurat$pseudo_bulk == temp]
    num <- sum(human_peak_Seurat@meta.data[temp,'ReadsInPeaks'])
    return(num)
  }else if(temp[1] == 'macaque'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(macaque_peak_Seurat)[macaque_peak_Seurat$pseudo_bulk == temp]
    num <- sum(macaque_peak_Seurat@meta.data[temp,'ReadsInPeaks'])
    return(num)
  }else if(temp[1] == 'mouse'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(mouse_peak_Seurat)[mouse_peak_Seurat$pseudo_bulk == temp]
    num <- sum(mouse_peak_Seurat@meta.data[temp,'ReadsInPeaks'])
    return(num)
  }else{
    stop('error!')
  }
}))

meta_data$Reads_In_Cell_Type_Peaks <- unlist(base::lapply(X = rownames(meta_data),FUN = function(x){
  temp <- strsplit(x = x,split = '#')[[1]]
  if(temp[1] == 'human'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(human_peak_Seurat)[human_peak_Seurat$pseudo_bulk == temp]
    num <- sum(human_peak_Seurat@meta.data[temp,'Reads_In_Cell_Type_Peaks'])
    return(num)
  }else if(temp[1] == 'macaque'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(macaque_peak_Seurat)[macaque_peak_Seurat$pseudo_bulk == temp]
    num <- sum(macaque_peak_Seurat@meta.data[temp,'Reads_In_Cell_Type_Peaks'])
    return(num)
  }else if(temp[1] == 'mouse'){
    temp <- paste(temp[2],temp[3],sep = '#')
    temp <- colnames(mouse_peak_Seurat)[mouse_peak_Seurat$pseudo_bulk == temp]
    num <- sum(mouse_peak_Seurat@meta.data[temp,'Reads_In_Cell_Type_Peaks'])
    return(num)
  }else{
    stop('error!')
  }
}))

Brain_peak_matrix <- CreateSeuratObject(counts = Brain_peak_matrix,project = 'ATAC',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
saveRDS(Brain_peak_matrix,file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')

# sample distance analysis ------------------------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')

#subset dataset
Brain_peak_matrix <- Brain_peak_matrix[,!(Brain_peak_matrix$cell_type %in% c('Cycling','End','End/Per','Ependymal','Mic','OPC','Per'))]

#normalize by reads in cell type peaks
temp <- Brain_peak_matrix@assays$RNA@counts
temp <- do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  return((temp[,x]/Brain_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"])*median(Brain_peak_matrix$Reads_In_Cell_Type_Peaks))
}))
colnames(temp) <- colnames(Brain_peak_matrix@assays$RNA@counts)
rownames(temp) <- rownames(Brain_peak_matrix@assays$RNA@counts)

temp <- log1p(temp)
Brain_peak_matrix@assays$RNA@data <- temp

Brain_peak_matrix <- FindVariableFeatures(object = Brain_peak_matrix,assay = 'RNA',selection.method = 'disp',nfeatures = 10000,verbose = TRUE)
Brain_peak_matrix <- ScaleData(object = Brain_peak_matrix,features = VariableFeatures(Brain_peak_matrix),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)

#PCA
Brain_peak_matrix <- RunPCA(object = Brain_peak_matrix,assay = 'RNA',features = VariableFeatures(Brain_peak_matrix),npcs = 50,verbose = TRUE)
temp <- Brain_peak_matrix@meta.data
temp <- cbind(temp,Brain_peak_matrix@reductions$pca@cell.embeddings)
pdf(file = './res/step_35_fig_220606/Brain_peak_matrix_PCA_plot_with_10000_features.pdf',width = 6,height = 5)
ggplot(temp,aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'PCA plot')
dev.off()

# export cell type coverage file ------------------------------------------

## human -------------------------------------------------------------------

#load peak matrix
human_peak_matrix <- readRDS(file = '/data/User/sunym/project/Brain/res/step_35_fig_220606/human_peak_matrix_Seurat.rds')

#create ArchR project
temp
ArrowFiles
file.exists(ArrowFiles)
Greenleaf_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = './Greenleaf_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)
Greenleaf_ATAC_ArchR <- Greenleaf_ATAC_ArchR[colnames(human_peak_matrix)]

#add meta data
Greenleaf_ATAC_ArchR$nFrags_matrix <- human_peak_matrix$nFrags
Greenleaf_ATAC_ArchR$ReadsInTSS_matrix <- human_peak_matrix$ReadsInTSS
Greenleaf_ATAC_ArchR$cell_type <- human_peak_matrix$cell_type
Greenleaf_ATAC_ArchR$Reads_In_Cell_Type_Peaks <- human_peak_matrix$Reads_In_Cell_Type_Peaks
Greenleaf_ATAC_ArchR$ReadsInPeaks <- human_peak_matrix$ReadsInPeaks

#export
table(Greenleaf_ATAC_ArchR$cell_type)
table(Greenleaf_ATAC_ArchR$nFrags_matrix == Greenleaf_ATAC_ArchR$nFrags)
table(Greenleaf_ATAC_ArchR$ReadsInTSS_matrix == Greenleaf_ATAC_ArchR$ReadsInTSS)

Greenleaf_ATAC_ArchR$ReadsInTSS <- Greenleaf_ATAC_ArchR$Reads_In_Cell_Type_Peaks
coverage_file <- getGroupBW(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)

Greenleaf_ATAC_ArchR$ReadsInTSS <- Greenleaf_ATAC_ArchR$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)

## macaque -----------------------------------------------------------------

#load peak matrix
macaque_peak_matrix <- readRDS(file = '/data/User/sunym/project/Brain/res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')

#create ArchR project
temp
ArrowFiles
file.exists(ArrowFiles)
macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = './macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)
macaque_multiome_ArchR <- macaque_multiome_ArchR[colnames(macaque_peak_matrix)]

#add meta data
macaque_multiome_ArchR$nFrags_matrix <- macaque_peak_matrix$nFrags
macaque_multiome_ArchR$ReadsInTSS_matrix <- macaque_peak_matrix$ReadsInTSS
macaque_multiome_ArchR$cell_type <- macaque_peak_matrix$cell_type
macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks <- macaque_peak_matrix$Reads_In_Cell_Type_Peaks
macaque_multiome_ArchR$ReadsInPeaks <- macaque_peak_matrix$ReadsInPeaks

#export
table(macaque_multiome_ArchR$cell_type)
table(macaque_multiome_ArchR$nFrags_matrix == macaque_multiome_ArchR$nFrags)
table(macaque_multiome_ArchR$ReadsInTSS_matrix == macaque_multiome_ArchR$ReadsInTSS)

macaque_multiome_ArchR$ReadsInTSS <- macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
coverage_file <- getGroupBW(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)

macaque_multiome_ArchR$ReadsInTSS <- macaque_multiome_ArchR$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)

## mouse -------------------------------------------------------------------

#load peak matrix
mouse_peak_matrix <- readRDS(file = '/data/User/sunym/project/Brain/res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

#create ArchR project
temp
ArrowFiles
file.exists(ArrowFiles)
mouse_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = './mouse_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)
mouse_ATAC_ArchR <- mouse_ATAC_ArchR[colnames(mouse_peak_matrix)]

#add meta data
mouse_ATAC_ArchR$nFrags_matrix <- mouse_peak_matrix$nFrags
mouse_ATAC_ArchR$ReadsInTSS_matrix <- mouse_peak_matrix$ReadsInTSS
mouse_ATAC_ArchR$cell_type <- mouse_peak_matrix$cell_type
mouse_ATAC_ArchR$Reads_In_Cell_Type_Peaks <- mouse_peak_matrix$Reads_In_Cell_Type_Peaks
mouse_ATAC_ArchR$ReadsInPeaks <- mouse_peak_matrix$ReadsInPeaks

#export
table(mouse_ATAC_ArchR$cell_type)
table(mouse_ATAC_ArchR$nFrags_matrix == mouse_ATAC_ArchR$nFrags)
table(mouse_ATAC_ArchR$ReadsInTSS_matrix == mouse_ATAC_ArchR$ReadsInTSS)

mouse_ATAC_ArchR$ReadsInTSS <- mouse_ATAC_ArchR$Reads_In_Cell_Type_Peaks
coverage_file <- getGroupBW(ArchRProj = mouse_ATAC_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)

mouse_ATAC_ArchR$ReadsInTSS <- mouse_ATAC_ArchR$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = mouse_ATAC_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)
