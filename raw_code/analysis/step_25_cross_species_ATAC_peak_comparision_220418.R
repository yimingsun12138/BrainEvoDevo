#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: cross species ATAC peak comparision                             ##
## Data: 2022.04.18                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')

# macaque liftover to human -----------------------------------------------
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

#get peakset
macaque_peakset <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
names(macaque_peakset) <- paste(as.character(macaque_peakset@seqnames),as.character(macaque_peakset@ranges),sep = '-')
human_peakset <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
names(human_peakset) <- paste(as.character(human_peakset@seqnames),as.character(human_peakset@ranges),sep = '-')

#liftover to human
macaque_peakset_lifted <- my_unique_peakset_liftover(ori_GRanges = macaque_peakset,
                                                     UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                     chain_file = './data/reference/rheMac10ToHg38.over.chain',
                                                     liftOver_mismatch = 0.1,length_mismatch = 0.1,
                                                     mapped_chr = unique(as.character(human_peakset@seqnames)),
                                                     tmp_path = './res/step_25_fig_220418/tmp_220420_1243')

#merge with human peakset
macaque_human_merged_peakset <- my_bedtools_merge(peakset_x = macaque_peakset_lifted$mapped,
                                                  peakset_y = human_peakset,
                                                  bedtools_path = '/data/User/sunym/software/bedtools/bedtools',
                                                  d = 0,tmp_path = './res/step_25_fig_220418/tmp_220420_1314')

#re-liftover to macaque
macaque_human_merged_peakset <- my_unique_peakset_liftover(ori_GRanges = macaque_human_merged_peakset,
                                                           UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                           chain_file = './data/reference/hg38ToRheMac10.over.chain',
                                                           liftOver_mismatch = 0.1,length_mismatch = 0.1,
                                                           mapped_chr = unique(as.character(macaque_peakset@seqnames)),
                                                           tmp_path = './res/step_25_fig_220418/tmp_220420_1317')

# mouse liftover to human -------------------------------------------------
p1 <- plotEmbedding(ArchRProj = mouse_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

#get peakset
mouse_peakset <- getPeakSet(ArchRProj = mouse_ATAC_ArchR)
names(mouse_peakset) <- paste(as.character(mouse_peakset@seqnames),as.character(mouse_peakset@ranges),sep = '-')
mouse_peakset$peakname <- names(mouse_peakset)

#liftover to human
# mouse_peakset_lifted <- my_unique_peakset_liftover(ori_GRanges = mouse_peakset,
#                                                    UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
#                                                    chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/mm10ToHg38.over.chain',
#                                                    liftOver_mismatch = 0.1,length_mismatch = 0.1,
#                                                    mapped_chr = unique(as.character(human_peakset@seqnames)),
#                                                    tmp_path = './res/step_25_fig_220418/tmp_220420_1556')

mouse_peakset_lifted <- my_unique_peakset_liftover(ori_GRanges = mouse_peakset,
                                                   UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                   chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/mm10ToHg38.over.chain',
                                                   liftOver_mismatch = 0.8,length_mismatch = 0.8,
                                                   mapped_chr = unique(as.character(human_peakset@seqnames)),
                                                   tmp_path = './res/step_25_fig_220418/tmp_220420_1556')

#merge with macaque_human_merged_peakset
macaque_human_mouse_merged_peakset <- my_bedtools_merge(peakset_x = mouse_peakset_lifted$mapped,
                                                        peakset_y = macaque_human_merged_peakset$ori,
                                                        bedtools_path = '/data/User/sunym/software/bedtools/bedtools',
                                                        d = 0,tmp_path = './res/step_25_fig_220418/tmp_220420_1634')

#relift to mouse
Brain_ATAC_peakset_mouse <- my_unique_peakset_liftover(ori_GRanges = macaque_human_mouse_merged_peakset,
                                                       UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                       chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/hg38ToMm10.over.chain',
                                                       liftOver_mismatch = 0.8,length_mismatch = 0.8,
                                                       mapped_chr = unique(as.character(mouse_peakset@seqnames)),
                                                       tmp_path = './res/step_25_fig_220418/tmp_220420_1641')

#relift to macaque
Brain_ATAC_peakset_macaque <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset_mouse$ori,
                                                         UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                         chain_file = './data/reference/hg38ToRheMac10.over.chain',
                                                         liftOver_mismatch = 0.1,length_mismatch = 0.1,
                                                         mapped_chr = unique(as.character(macaque_peakset@seqnames)),
                                                         tmp_path = './res/step_25_fig_220418/tmp_220420_1756')
#save data
Brain_ATAC_peakset <- SimpleList(
  human = Brain_ATAC_peakset_macaque$ori,
  macaque = Brain_ATAC_peakset_macaque$mapped,
  mouse = Brain_ATAC_peakset_mouse$mapped[Brain_ATAC_peakset_mouse$mapped$name %in% names(Brain_ATAC_peakset_macaque$ori)]
)
saveRDS(Brain_ATAC_peakset,file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

# create human peak matrix ------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

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
Greenleaf_ATAC_ArchR$cell_type <- temp@cellColData[rownames(Greenleaf_ATAC_ArchR@cellColData),"cell_type"]
Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = Brain_ATAC_peakset$human,
                                   genomeAnnotation = getGenomeAnnotation(Greenleaf_ATAC_ArchR),force = TRUE)
Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR,verbose = TRUE,force = TRUE)
human_peak_matrix <- getMatrixFromProject(ArchRProj = Greenleaf_ATAC_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

#create human peak matrix
meta_data <- as.data.frame(human_peak_matrix@colData)
temp <- human_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- paste(human_peak_matrix@rowRanges@seqnames,as.character(human_peak_matrix@rowRanges@ranges),sep = '-')
human_peak_matrix <- CreateSeuratObject(counts = temp,project = 'human',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

temp <- readRDS(file = '/data/User/sunym/project/Brain/data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
rep_list <- base::lapply(names(temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups),FUN = function(x){
  temp <- names(temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[x]])
  return(paste(x,temp,sep = '#'))
})
rep_list <- unlist(rep_list)

temp <- do.call(what = cbind,args = base::lapply(X = rep_list,FUN = function(x){
  arg <- strsplit(x = x,split = '#')[[1]]
  cell_list <- temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[arg[1]]][[arg[2]]]
  temp <- rowSums(human_peak_matrix@assays$RNA@counts[,cell_list])
  temp <- data.frame(temp)
  colnames(temp) <- paste('human',x,sep = '#')
  rownames(temp) <- rownames(human_peak_matrix@assays$RNA@counts)
  return(temp)
}))

#save data
human_peak_matrix <- temp
saveRDS(human_peak_matrix,file = './res/step_25_fig_220418/human_peak_matrix.rds')

# create macaque peak matrix ------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

#create new ArchR object
temp <- readRDS(file = '/data/User/sunym/project/Brain/ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
ArrowFiles <- list.files('/data/User/sunym/project/Brain/ArchR/arrow_file/macaque_multiome_data/')
ArrowFiles <- paste('/data/User/sunym/project/Brain/ArchR/arrow_file/macaque_multiome_data',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/temp/macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

macaque_multiome_ArchR <- macaque_multiome_ArchR[rownames(macaque_multiome_ArchR@cellColData) %in% rownames(temp@cellColData)]
macaque_multiome_ArchR$cell_type <- temp@cellColData[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]
macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = Brain_ATAC_peakset$macaque,
                                   genomeAnnotation = getGenomeAnnotation(macaque_multiome_ArchR),force = TRUE)
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)
macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

#create human peak matrix
meta_data <- as.data.frame(macaque_peak_matrix@colData)
temp <- macaque_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- paste(macaque_peak_matrix@rowRanges@seqnames,as.character(macaque_peak_matrix@rowRanges@ranges),sep = '-')
macaque_peak_matrix <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

temp <- readRDS(file = '/data/User/sunym/project/Brain/ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
rep_list <- base::lapply(names(temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups),FUN = function(x){
  temp <- names(temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[x]])
  return(paste(x,temp,sep = '#'))
})
rep_list <- unlist(rep_list)

temp <- do.call(what = cbind,args = base::lapply(X = rep_list,FUN = function(x){
  arg <- strsplit(x = x,split = '#')[[1]]
  cell_list <- temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[arg[1]]][[arg[2]]]
  temp <- rowSums(macaque_peak_matrix@assays$RNA@counts[,cell_list])
  temp <- data.frame(temp)
  colnames(temp) <- paste('macaque',x,sep = '#')
  rownames(temp) <- rownames(macaque_peak_matrix@assays$RNA@counts)
  return(temp)
}))

#save data
macaque_peak_matrix <- temp
saveRDS(macaque_peak_matrix,file = './res/step_25_fig_220418/macaque_peak_matrix.rds')

# create mouse peak matrix ------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

#create new ArchR object
temp <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')
ArrowFiles <- list.files('./data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/arrow_file/')
ArrowFiles <- paste('./data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/arrow_file',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

mouse_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/temp/mouse_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

mouse_ATAC_ArchR <- mouse_ATAC_ArchR[rownames(mouse_ATAC_ArchR@cellColData) %in% rownames(temp@cellColData)]
mouse_ATAC_ArchR$cell_type <- temp@cellColData[rownames(mouse_ATAC_ArchR@cellColData),"cell_type"]
mouse_ATAC_ArchR <- addPeakSet(ArchRProj = mouse_ATAC_ArchR,peakSet = Brain_ATAC_peakset$mouse,
                               genomeAnnotation = getGenomeAnnotation(mouse_ATAC_ArchR),force = TRUE)
mouse_ATAC_ArchR <- addPeakMatrix(ArchRProj = mouse_ATAC_ArchR,verbose = TRUE,force = TRUE)
mouse_peak_matrix <- getMatrixFromProject(ArchRProj = mouse_ATAC_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

#create human peak matrix
meta_data <- as.data.frame(mouse_peak_matrix@colData)
temp <- mouse_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- paste(mouse_peak_matrix@rowRanges@seqnames,as.character(mouse_peak_matrix@rowRanges@ranges),sep = '-')
mouse_peak_matrix <- CreateSeuratObject(counts = temp,project = 'human',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

temp <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')
rep_list <- base::lapply(names(temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups),FUN = function(x){
  temp <- names(temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[x]])
  return(paste(x,temp,sep = '#'))
})
rep_list <- unlist(rep_list)

temp <- do.call(what = cbind,args = base::lapply(X = rep_list,FUN = function(x){
  arg <- strsplit(x = x,split = '#')[[1]]
  cell_list <- temp@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[arg[1]]][[arg[2]]]
  temp <- rowSums(mouse_peak_matrix@assays$RNA@counts[,cell_list])
  temp <- data.frame(temp)
  colnames(temp) <- paste('mouse',x,sep = '#')
  rownames(temp) <- rownames(mouse_peak_matrix@assays$RNA@counts)
  return(temp)
}))

#save data
mouse_peak_matrix <- temp
saveRDS(mouse_peak_matrix,file = './res/step_25_fig_220418/mouse_peak_matrix.rds')

# create Brain peak matrix ------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')
human_peak_matrix <- readRDS(file = './res/step_25_fig_220418/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_25_fig_220418/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_25_fig_220418/mouse_peak_matrix.rds')

#rename peak matrix
names(Brain_ATAC_peakset$macaque) <- paste(Brain_ATAC_peakset$macaque@seqnames,as.character(Brain_ATAC_peakset$macaque@ranges),sep = '-')
rownames(macaque_peak_matrix) <- Brain_ATAC_peakset$macaque[rownames(macaque_peak_matrix)]$name
names(Brain_ATAC_peakset$mouse) <- paste(Brain_ATAC_peakset$mouse@seqnames,as.character(Brain_ATAC_peakset$mouse@ranges),sep = '-')
rownames(mouse_peak_matrix) <- Brain_ATAC_peakset$mouse[rownames(mouse_peak_matrix)]$name

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

Brain_peak_matrix <- CreateSeuratObject(counts = Brain_peak_matrix,project = 'ATAC',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
saveRDS(Brain_peak_matrix,file = './res/step_25_fig_220418/Brain_peak_matrix.rds')

# sample distance analysis ------------------------------------------------
#load data
Brain_peak_matrix <- readRDS(file = './res/step_25_fig_220418/Brain_peak_matrix.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

Brain_peak_matrix <- Brain_peak_matrix[,!(Brain_peak_matrix$cell_type %in% c('Cycling','End','End/Per','Ependymal','Mic','OPC','Per'))]
Brain_peak_matrix <- NormalizeData(object = Brain_peak_matrix,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
Brain_peak_matrix <- FindVariableFeatures(object = Brain_peak_matrix,assay = 'RNA',selection.method = 'dispersion',nfeatures = 50000,verbose = TRUE)

#PCA
Brain_peak_matrix <- ScaleData(object = Brain_peak_matrix,features = VariableFeatures(Brain_peak_matrix),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE)
Brain_peak_matrix <- RunPCA(object = Brain_peak_matrix,assay = 'RNA',features = VariableFeatures(Brain_peak_matrix),npcs = 50,verbose = TRUE)

temp <- Brain_peak_matrix@meta.data
temp <- cbind(temp,Brain_peak_matrix@reductions$pca@cell.embeddings)
ggplot(temp,aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + theme(aspect.ratio = 1)

#UMAP
Brain_peak_matrix <- RunUMAP(object = Brain_peak_matrix,dims = 1:8,reduction = 'pca')
DimPlot(object = Brain_peak_matrix,group.by = 'cell_type',label = TRUE,repel = TRUE)

temp <- Brain_peak_matrix@meta.data
temp <- cbind(temp,Brain_peak_matrix@reductions$umap@cell.embeddings)

pdf(file = './res/step_25_fig_220418/Brain_peak_matrix_umap_dimplot.pdf',width = 8,height = 7)
ggplot(temp,aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'UMAP plot')
dev.off()

#cluster analysis
mat <- Brain_peak_matrix[,Brain_peak_matrix$cell_type %in% c('Ex-3','Ex-2','Ex-1','IP','RG','InCGE','InMGE')]
group <- kmeans(mat@assays$RNA@scale.data,centers = 10)$cluster
col_fun = colorRamp2(c(-2, 0, 4), c("green", "white", "red"))
col_anno <- HeatmapAnnotation(cell_type = mat$cell_type,which = 'column')

pdf(file = './res/step_25_fig_220418/Brain_peak_matrix_variable_peak_pseudobulk_heatmap.pdf',width = 14,height = 8)
Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
        show_column_names = TRUE,show_row_names = FALSE,
        cluster_rows = cluster_within_group(mat = t(as.matrix(mat@assays$RNA@scale.data)),factor = group),
        cluster_columns = cluster_within_group(mat = as.matrix(mat@assays$RNA@scale.data),factor = mat$cell_type),
        col = col_fun,height = unit(4,'inches'),width = unit(12,'inches'),name = 'z-score',top_annotation = col_anno)
dev.off()

# RG analysis -------------------------------------------------------------
Brain_peak_matrix <- readRDS(file = './res/step_25_fig_220418/Brain_peak_matrix.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

Brain_peak_matrix <- Brain_peak_matrix[,!(Brain_peak_matrix$cell_type %in% c('Cycling','End','End/Per','Ependymal','Mic','OPC','Per'))]
gene_list <- FindMarkers(object = Brain_peak_matrix,ident.1 = 'RG',group.by = 'cell_type',assay = 'RNA',slot = 'counts',logfc.threshold = 1,test.use = 'DESeq2',verbose = TRUE,only.pos = TRUE)
RG_peak_matrix <- Brain_peak_matrix[rownames(gene_list),Brain_peak_matrix$cell_type == 'RG']
RG_peak_matrix <- ScaleData(object = RG_peak_matrix,features = rownames(RG_peak_matrix),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE)

mat <- RG_peak_matrix
group <- kmeans(mat@assays$RNA@scale.data,centers = 8)$cluster
col_anno <- HeatmapAnnotation(cell_type = mat$species,which = 'column')
col_fun = colorRamp2(c(-2,0,2), c("green", "white", "red"))

pdf(file = './res/step_25_fig_220418/RG_peak_matrix_heatmap.pdf',width = 8,height = 16)
Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
        show_column_names = TRUE,show_row_names = FALSE,
        cluster_rows = FALSE,
        row_split = factor(group,levels = c(5,1,4,3,2,8,6,7)),
        cluster_columns = TRUE,clustering_method_columns = 'ward.D',
        height = unit(12,'inches'),width = unit(6,'inches'),
        name = 'z-score',top_annotation = col_anno,col = col_fun)
dev.off()

saveRDS(group,file = './res/step_25_fig_220418/RG_marker_kmeans_group.rds')

# IP analysis -------------------------------------------------------------
Brain_peak_matrix <- readRDS(file = './res/step_25_fig_220418/Brain_peak_matrix.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

Brain_peak_matrix <- Brain_peak_matrix[,!(Brain_peak_matrix$cell_type %in% c('Cycling','End','End/Per','Ependymal','Mic','OPC','Per'))]
gene_list <- FindMarkers(object = Brain_peak_matrix,ident.1 = 'IP',group.by = 'cell_type',assay = 'RNA',slot = 'counts',logfc.threshold = 1,test.use = 'DESeq2',verbose = TRUE,only.pos = TRUE)
IP_peak_matrix <- Brain_peak_matrix[rownames(gene_list),Brain_peak_matrix$cell_type == 'IP']
IP_peak_matrix <- ScaleData(object = IP_peak_matrix,features = rownames(IP_peak_matrix),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE)

mat <- IP_peak_matrix
group <- kmeans(mat@assays$RNA@scale.data,centers = 10)$cluster
col_anno <- HeatmapAnnotation(cell_type = mat$species,which = 'column')
col_fun = colorRamp2(c(-2,0,2), c("green", "white", "red"))

pdf(file = './res/step_25_fig_220418/IP_peak_matrix_heatmap.pdf',width = 8,height = 16)
Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
        show_column_names = TRUE,show_row_names = FALSE,
        cluster_rows = FALSE,
        row_split = factor(group,levels = 1:10),
        cluster_columns = TRUE,clustering_method_columns = 'ward.D',
        height = unit(12,'inches'),width = unit(6,'inches'),
        name = 'z-score',top_annotation = col_anno,col = col_fun)
dev.off()

# Ex_3 analysis -------------------------------------------------------------
Brain_peak_matrix <- readRDS(file = './res/step_25_fig_220418/Brain_peak_matrix.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_25_fig_220418/Brain_ATAC_Peakset.rds')

Brain_peak_matrix <- Brain_peak_matrix[,!(Brain_peak_matrix$cell_type %in% c('Cycling','End','End/Per','Ependymal','Mic','OPC','Per'))]
gene_list <- FindMarkers(object = Brain_peak_matrix,ident.1 = 'Ex-3',group.by = 'cell_type',assay = 'RNA',slot = 'counts',logfc.threshold = 1,test.use = 'DESeq2',verbose = TRUE,only.pos = TRUE)
Ex_3_peak_matrix <- Brain_peak_matrix[rownames(gene_list),Brain_peak_matrix$cell_type == 'Ex-3']
Ex_3_peak_matrix <- ScaleData(object = Ex_3_peak_matrix,features = rownames(Ex_3_peak_matrix),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE)

mat <- Ex_3_peak_matrix
group <- kmeans(mat@assays$RNA@scale.data,centers = 10)$cluster
col_anno <- HeatmapAnnotation(cell_type = mat$species,which = 'column')
col_fun = colorRamp2(c(-2,0,2), c("green", "white", "red"))

pdf(file = './res/step_25_fig_220418/Ex_3_peak_matrix_heatmap.pdf',width = 8,height = 16)
Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
        show_column_names = TRUE,show_row_names = FALSE,
        cluster_rows = FALSE,
        row_split = factor(group,levels = 1:10),
        cluster_columns = TRUE,clustering_method_columns = 'ward.D',
        height = unit(12,'inches'),width = unit(6,'inches'),
        name = 'z-score',top_annotation = col_anno,col = col_fun)
dev.off()