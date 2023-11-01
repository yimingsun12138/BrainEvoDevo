#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: validate species DAP method                                     ##
## Data: 2022.05.25                                                                ##
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

# try macaque 0.1/0.1 mouse 0.4/0.4 ---------------------------------------
#load data
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')

human_peakset <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
macaque_peakset <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
mouse_peakset <- getPeakSet(ArchRProj = mouse_ATAC_ArchR)

#chack data
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)
#all checked!

## macaque lift to human ---------------------------------------------------
lifted_macaque_peakset <- my_unique_peakset_liftover(ori_GRanges = macaque_peakset,
                                                     UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                     chain_file = './data/reference/rheMac10ToHg38.over.chain',
                                                     liftOver_mismatch = 0.1,length_filter = TRUE,length_mismatch = 0.1,
                                                     chr_filter = TRUE,mapped_chr = unique(as.character(human_peakset@seqnames)),
                                                     overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220525/')

#merge with human peakset
Brain_ATAC_peakset <- my_bedtools_merge(peakset_x = lifted_macaque_peakset$mapped,peakset_y = human_peakset,
                                        bedtools_path = '/data/User/sunym/software/bedtools/bedtools',
                                        d = 0,bedtools_param = NULL,tmp_path = '/data/User/sunym/temp/tmp_220525')

#re-lift to macaque
Brain_ATAC_peakset <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset,
                                                 UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                 chain_file = './data/reference/hg38ToRheMac10.over.chain',
                                                 liftOver_mismatch = 0.1,length_filter = TRUE,length_mismatch = 0.1,
                                                 chr_filter = TRUE,mapped_chr = unique(as.character(macaque_peakset@seqnames)),
                                                 overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220525')

## mouse lift to human -----------------------------------------------------
lifted_mouse_peakset <- my_unique_peakset_liftover(ori_GRanges = mouse_peakset,
                                                   UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                   chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/mm10ToHg38.over.chain',
                                                   liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                   chr_filter = TRUE,mapped_chr = unique(as.character(human_peakset@seqnames)),
                                                   overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220525')

#merge with Brain_ATAC_peakset
Brain_ATAC_peakset <- my_bedtools_merge(peakset_x = Brain_ATAC_peakset$ori,peakset_y = lifted_mouse_peakset$mapped,
                                        bedtools_path = '/data/User/sunym/software/bedtools/bedtools',
                                        d = 0,bedtools_param = NULL,tmp_path = '/data/User/sunym/temp/tmp_220525')


## create unique Brain_ATAC_peakset ----------------------------------------
#re-lift to mouse
Brain_ATAC_peakset <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset,
                                                 UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                 chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/hg38ToMm10.over.chain',
                                                 liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                 chr_filter = TRUE,mapped_chr = unique(as.character(mouse_peakset@seqnames)),
                                                 overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220525')

#re-lift to macaque
Brain_ATAC_peakset_macaque <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset$ori,
                                                         UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                         chain_file = './data/reference/hg38ToRheMac10.over.chain',
                                                         liftOver_mismatch = 0.1,length_filter = TRUE,length_mismatch = 0.1,
                                                         chr_filter = TRUE,mapped_chr = unique(as.character(macaque_peakset@seqnames)),
                                                         overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220525')

#re-lift to mouse
Brain_ATAC_peakset_mouse <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset_macaque$ori,
                                                       UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                       chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/hg38ToMm10.over.chain',
                                                       liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                       chr_filter = TRUE,mapped_chr = unique(as.character(mouse_peakset@seqnames)),
                                                       overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220525')

#save data
table(names(Brain_ATAC_peakset_macaque$ori) == names(Brain_ATAC_peakset_mouse$ori))
Brain_ATAC_peakset <- SimpleList(
  human = Brain_ATAC_peakset_macaque$ori,
  macaque = Brain_ATAC_peakset_macaque$mapped,
  mouse = Brain_ATAC_peakset_mouse$mapped
)
saveRDS(Brain_ATAC_peakset,file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')

# create peak matrix ------------------------------------------------------
#all run on 202.205.131.32, ArchR is fucking shit.
## create human peak matrix ------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')

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

human_peak_matrix$pseudo_bulk <- NA
for (i in unique(human_peak_matrix$cell_type)) {
  j <- sum(human_peak_matrix$cell_type == i)
  pseudo_bulk <- sample(rep(c('rep1','rep2','rep3'),times = c(round(j/3),round(j/3),(j-2*round(j/3)))))
  pseudo_bulk <- paste(i,pseudo_bulk,sep = '#')
  human_peak_matrix@meta.data[human_peak_matrix$cell_type == i,'pseudo_bulk'] <- pseudo_bulk
}
human_peak_matrix <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'RNA',return.seurat = FALSE,
                                                 group.by = 'pseudo_bulk',slot = 'counts',verbose = TRUE)
human_peak_matrix <- human_peak_matrix$RNA

#save data
saveRDS(human_peak_matrix,file = './res/step_31_fig_220525/human_peak_matrix.rds')

## create macaque peak matrix ------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')

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
macaque_multiome_ArchR$cell_type <- temp@cellColData[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]
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
macaque_peak_matrix <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'RNA',return.seurat = FALSE,
                                                   group.by = 'pseudo_bulk',slot = 'counts',verbose = TRUE)
macaque_peak_matrix <- macaque_peak_matrix$RNA

#save data
saveRDS(macaque_peak_matrix,file = './res/step_31_fig_220525/macaque_peak_matrix.rds')

## create mouse peak matrix ------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')

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
mouse_ATAC_ArchR$cell_type <- temp@cellColData[rownames(mouse_ATAC_ArchR@cellColData),"cell_type"]
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
mouse_peak_matrix <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'RNA',return.seurat = FALSE,
                                                 group.by = 'pseudo_bulk',slot = 'counts',verbose = TRUE)
mouse_peak_matrix <- mouse_peak_matrix$RNA

#save data
saveRDS(mouse_peak_matrix,file = './res/step_31_fig_220525/mouse_peak_matrix.rds')

# create Brain_peak_matrix ------------------------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')
human_peak_matrix <- readRDS(file = './res/step_31_fig_220525/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_31_fig_220525/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_31_fig_220525/mouse_peak_matrix.rds')

#rename peak matrix
names(Brain_ATAC_peakset$macaque) <- paste(Brain_ATAC_peakset$macaque@seqnames,as.character(Brain_ATAC_peakset$macaque@ranges),sep = '-')
rownames(macaque_peak_matrix) <- Brain_ATAC_peakset$macaque[rownames(macaque_peak_matrix)]$name
names(Brain_ATAC_peakset$mouse) <- paste(Brain_ATAC_peakset$mouse@seqnames,as.character(Brain_ATAC_peakset$mouse@ranges),sep = '-')
rownames(mouse_peak_matrix) <- Brain_ATAC_peakset$mouse[rownames(mouse_peak_matrix)]$name

colnames(human_peak_matrix) <- paste('human',colnames(human_peak_matrix),sep = '#')
colnames(macaque_peak_matrix) <- paste('macaque',colnames(macaque_peak_matrix),sep = '#')
colnames(mouse_peak_matrix) <- paste('mouse',colnames(mouse_peak_matrix),sep = '#')
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

#save data
Brain_peak_matrix <- CreateSeuratObject(counts = Brain_peak_matrix,project = 'ATAC',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
saveRDS(Brain_peak_matrix,file = './res/step_31_fig_220525/Brain_peak_matrix_Seurat.rds')

# sample distance analysis ------------------------------------------------
#load data
Brain_peak_matrix <- readRDS(file = './res/step_31_fig_220525/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')

Brain_peak_matrix <- Brain_peak_matrix[,!(Brain_peak_matrix$cell_type %in% c('Cycling','End','End/Per','Ependymal','Mic','OPC','Per'))]
Brain_peak_matrix <- NormalizeData(object = Brain_peak_matrix,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
Brain_peak_matrix <- FindVariableFeatures(object = Brain_peak_matrix,assay = 'RNA',selection.method = 'dispersion',nfeatures = 3000,verbose = TRUE)

#PCA
Brain_peak_matrix <- ScaleData(object = Brain_peak_matrix,features = VariableFeatures(Brain_peak_matrix),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE)
Brain_peak_matrix <- RunPCA(object = Brain_peak_matrix,assay = 'RNA',features = VariableFeatures(Brain_peak_matrix),npcs = 50,verbose = TRUE)

temp <- Brain_peak_matrix@meta.data
temp <- cbind(temp,Brain_peak_matrix@reductions$pca@cell.embeddings)
ggplot(temp,aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + theme(aspect.ratio = 1)

#UMAP
Brain_peak_matrix <- RunUMAP(object = Brain_peak_matrix,dims = 1:22,reduction = 'pca')
DimPlot(object = Brain_peak_matrix,group.by = 'cell_type',label = TRUE,repel = TRUE)

temp <- Brain_peak_matrix@meta.data
temp <- cbind(temp,Brain_peak_matrix@reductions$umap@cell.embeddings)

pdf(file = './res/step_31_fig_220525/Brain_peak_matrix_umap_dimplot.pdf',width = 6,height = 5)
ggplot(temp,aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'UMAP plot')
dev.off()

#cluster analysis
hc <- hclust(d = dist(t(Brain_peak_matrix@assays$RNA@data[VariableFeatures(Brain_peak_matrix),])),method = 'ward.D')
pdf(file = './res/step_31_fig_220525/Brain_peak_matrix_hclust.pdf',width = 12,height = 6)
plot(hc)
dev.off()

# try RG analysis ---------------------------------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')
Brain_peak_matrix <- readRDS(file = './res/step_31_fig_220525/Brain_peak_matrix_Seurat.rds')

macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')

RG_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == 'RG']

#get peak_list
human_RG_peak <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/RG-reproduciblePeaks.gr.rds')
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_RG_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

macaque_RG_peak <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/RG-reproduciblePeaks.gr.rds')
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- countOverlaps(query = Brain_ATAC_peakset$macaque,subject = macaque_RG_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

mouse_RG_peak <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/RG-reproduciblePeaks.gr.rds')
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- countOverlaps(query = Brain_ATAC_peakset$mouse,subject = mouse_RG_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(RG_peak_matrix))

#macaque lfc
macaque_lfc <- FindMarkers(object = RG_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#mouse lfc
mouse_lfc <- FindMarkers(object = RG_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#peak lfc matrix
peak_lfc_matrix <- data.frame(peak = peak_list,macaque_lfc = 0,macaque_sig = 0,mouse_lfc = 0,mouse_sig = 0)
rownames(peak_lfc_matrix) <- peak_lfc_matrix$peak
peak_lfc_matrix[rownames(macaque_lfc),"macaque_lfc"] <- macaque_lfc$avg_log2FC
peak_lfc_matrix[rownames(macaque_lfc),"macaque_sig"] <- macaque_lfc$sig
peak_lfc_matrix[rownames(mouse_lfc),"mouse_lfc"] <- mouse_lfc$avg_log2FC
peak_lfc_matrix[rownames(mouse_lfc),"mouse_sig"] <- mouse_lfc$sig

#try lfc=1/0.7, sig=2
peak_lfc_matrix$group <- 'mixed'

temp <- rownames(macaque_lfc)[abs(macaque_lfc$avg_log2FC) > 1 & macaque_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'
temp <- rownames(mouse_lfc)[abs(mouse_lfc$avg_log2FC) > 1 & mouse_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'

peak_lfc_matrix[abs(peak_lfc_matrix$macaque_lfc) < 0.7 & 
                  abs(peak_lfc_matrix$mouse_lfc) < 0.7 & 
                  countOverlaps(query = Brain_ATAC_peakset$macaque[peak_lfc_matrix$peak],subject = macaque_RG_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$mouse[peak_lfc_matrix$peak],subject = mouse_RG_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$human[peak_lfc_matrix$peak],subject = human_RG_peak) > 0,
                "group"] <- 'conserved'

# #kmeans
# mat <- RG_peak_matrix
# mat <- NormalizeData(object = mat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
# VariableFeatures(mat) <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'differential']
# mat <- ScaleData(object = mat,features = VariableFeatures(mat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)

# group <- kmeans(mat@assays$RNA@scale.data[VariableFeatures(mat),],centers = 6)$cluster
# group <- readRDS(file = './res/step_29_fig_220513/RG_peak_kmeans_group.rds')
# ii <- group$peak_name
# group <- group$cluster
# names(group) <- ii
# col_anno <- HeatmapAnnotation(species = mat$species,which = 'column')
# 
# p1 <- Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
#               show_column_names = FALSE,show_row_names = FALSE,
#               cluster_rows = FALSE,row_split = factor(group[rownames(mat@assays$RNA@scale.data)],
#                                                       levels = unique(group)),
#               cluster_columns = FALSE,
#               column_split = factor(mat$species,levels = c('human','macaque','mouse')),
#               height = unit(9,'inches'),width = unit(6,'inches'),
#               name = 'z-score',top_annotation = col_anno,border = TRUE)
# 
# temp <- peak_lfc_matrix[rownames(p1@matrix)[p1@row_order],c("macaque_lfc","mouse_lfc")]
# p2 <- Heatmap(matrix = as.matrix(temp),
#               show_column_names = TRUE,show_row_names = FALSE,
#               cluster_rows = FALSE,row_split = factor(group[rownames(temp)],levels = unique(group)),
#               cluster_columns = FALSE,
#               height = unit(9,'inches'),width = unit(1,'inches'),
#               name = 'lfc',border = TRUE)
# 
# p1+p2
# 
# #modify groups
# temp <- peak_lfc_matrix[peak_lfc_matrix$group == 'conserved',]
# ii <- rep('conserved',times = nrow(temp))
# names(ii) <- temp$peak
# group <- append(group,ii)
# group <- data.frame(peak_name = names(group),raw_cluster = as.character(group))
# 
# #cluster 3 macaque_human_conserved
# #cluster 5 macaque_specific
# #cluster 6 mouse_macaque_conserved
# #cluster 2 mouse_specific
# #cluster 4 mouse_human_conserved
# #cluster 1 human_specific
# 
# group$cluster <- NA
# group[group$raw_cluster == '3',"cluster"] <- 'macaque_human_conserved'
# group[group$raw_cluster == '5',"cluster"] <- 'macaque_specific'
# group[group$raw_cluster == '6',"cluster"] <- 'mouse_macaque_conserved'
# group[group$raw_cluster == '2',"cluster"] <- 'mouse_specific'
# group[group$raw_cluster == '4',"cluster"] <- 'mouse_human_conserved'
# group[group$raw_cluster == '1',"cluster"] <- 'human_specific'
# group[group$raw_cluster == 'conserved',"cluster"] <- 'species_conserved'
# 
# saveRDS(group,file = './res/step_31_fig_220525/RG_peak_kmeans_group.rds')

#kmeans
mat <- RG_peak_matrix
mat <- NormalizeData(object = mat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
VariableFeatures(mat) <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'differential']
mat <- ScaleData(object = mat,features = VariableFeatures(mat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)

group <- readRDS(file = './res/step_31_fig_220525/RG_peak_kmeans_group.rds')
ii <- group$peak_name
group <- group$cluster
names(group) <- ii
col_anno <- HeatmapAnnotation(species = mat$species,which = 'column')

p1 <- Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
              show_column_names = FALSE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(mat@assays$RNA@scale.data)],
                                                      levels = c('human_specific','macaque_human_conserved',
                                                                 'mouse_human_conserved','macaque_specific',
                                                                 'mouse_macaque_conserved','mouse_specific',
                                                                 'species_conserved')),
              cluster_columns = FALSE,
              column_split = factor(mat$species,levels = c('human','macaque','mouse')),
              height = unit(10,'inches'),width = unit(6,'inches'),
              name = 'z-score',top_annotation = col_anno,border = TRUE,
              row_title = 'merged peakset')

temp <- peak_lfc_matrix[rownames(p1@matrix)[p1@row_order],c("macaque_lfc","mouse_lfc")]
p2 <- Heatmap(matrix = as.matrix(temp),
              show_column_names = TRUE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(temp)],levels = c('human_specific','macaque_human_conserved',
                                                                                       'mouse_human_conserved','macaque_specific',
                                                                                       'mouse_macaque_conserved','mouse_specific',
                                                                                       'species_conserved')),
              cluster_columns = FALSE,
              height = unit(10,'inches'),width = unit(1,'inches'),
              name = 'lfc',border = TRUE,row_title = 'merged peakset')

pdf(file = './res/step_31_fig_220525/RG_differential_peakset_z_score_heatmap.pdf',width = 10,height = 12)
p1+p2
dev.off()

#boxplot show lfc
peak_lfc_matrix <- peak_lfc_matrix[peak_lfc_matrix$group == 'differential',]
peak_lfc_matrix$cluster <- group[rownames(peak_lfc_matrix)]
temp <- data.frame(lfc = peak_lfc_matrix$macaque_lfc,species = 'macaque',cluster = peak_lfc_matrix$cluster)
temp <- rbind(temp,data.frame(lfc = peak_lfc_matrix$mouse_lfc,species = 'mouse',cluster = peak_lfc_matrix$cluster))

pdf(file = './res/step_31_fig_220525/RG_differential_peakset_lfc_boxplot.pdf',width = 7,height = 5)
ggplot(data = temp,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  facet_wrap(~cluster,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1))
dev.off()

# enrichheatmap
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/human_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- as.character(Brain_ATAC_peakset$macaque[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$macaque[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/macaque_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- as.character(Brain_ATAC_peakset$mouse[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$mouse[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/mouse_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.004,0.008),colors = c('#FFC30F','#C70039','#581845'))
pdf(file = './res/step_31_fig_220525/RG_peakset_CPM_center_besed_heatmap.pdf',width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                         'mouse_human_conserved','macaque_specific',
                                                                                         'mouse_macaque_conserved','mouse_specific',
                                                                                         'species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE)
dev.off()

#enrichment distribution
my_insertion_plot <- function(group=group,p1 = p1,p2 = p2,p3 = p3,group.by=group.by){
  temp <- p1@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'mouse'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5','mouse' = '#FDB31A')) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group.by) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_specific')
p_macaque_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_specific')
p_mouse_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_specific')
p_macaque_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_human_conserved')
p_mouse_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_human_conserved')
p_mouse_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

pdf(file = './res/step_31_fig_220525/RG_peakset_mean_insertion_dot_line_plot.pdf',width = 18,height = 7)
p_human_specific+p_macaque_specific+p_mouse_specific+p_macaque_human_conserved+
  p_mouse_human_conserved+p_mouse_macaque_conserved+p_species_conserved + plot_layout(ncol = 4)
dev.off()

#it worked!

# species contribution to final RG peakset --------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_31_fig_220525/Brain_ATAC_peak_macaque_0.1_mouse_0.4.rds')

#get RG peakset
human_RG_peak <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/RG-reproduciblePeaks.gr.rds')
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_RG_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

macaque_RG_peak <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/RG-reproduciblePeaks.gr.rds')
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- countOverlaps(query = Brain_ATAC_peakset$macaque,subject = macaque_RG_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

mouse_RG_peak <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/RG-reproduciblePeaks.gr.rds')
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- countOverlaps(query = Brain_ATAC_peakset$mouse,subject = mouse_RG_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% names(Brain_ATAC_peakset$human))

#show the contribution of these peak
#human
human_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$human[peak_list],subject = human_RG_peak)
human_contributed_peak <- human_contributed_peak[human_contributed_peak > 0]
human_contributed_peak <- names(human_contributed_peak)

#macaque
macaque_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$macaque[peak_list],subject = macaque_RG_peak)
macaque_contributed_peak <- macaque_contributed_peak[macaque_contributed_peak > 0]
macaque_contributed_peak <- names(macaque_contributed_peak)

#mouse
mouse_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$mouse[peak_list],subject = mouse_RG_peak)
mouse_contributed_peak <- mouse_contributed_peak[mouse_contributed_peak > 0]
mouse_contributed_peak <- names(mouse_contributed_peak)

#ggplot
temp <- list(human = human_contributed_peak,macaque = macaque_contributed_peak,mouse = mouse_contributed_peak)

pdf(file = './res/step_31_fig_220525/species_contribution_to_final_RG_peakset.pdf',width = 6.5,height = 6.5)
ggvenn(data = temp,c('human','macaque','mouse'))
dev.off()

# validate species conserved peakset parameter ----------------------------
RG_peak_matrix <- readRDS(file = './res/step_31_fig_220525/Brain_peak_matrix_Seurat.rds')
RG_peak_matrix <- RG_peak_matrix[,RG_peak_matrix$cell_type == 'RG']

table(peak_list %in% rownames(RG_peak_matrix))

#macaque lfc
macaque_lfc <- FindMarkers(object = RG_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#mouse lfc
mouse_lfc <- FindMarkers(object = RG_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#peak lfc matrix
peak_lfc_matrix <- data.frame(peak = peak_list,macaque_lfc = 0,macaque_sig = 0,mouse_lfc = 0,mouse_sig = 0)
rownames(peak_lfc_matrix) <- peak_lfc_matrix$peak
peak_lfc_matrix[rownames(macaque_lfc),"macaque_lfc"] <- macaque_lfc$avg_log2FC
peak_lfc_matrix[rownames(macaque_lfc),"macaque_sig"] <- macaque_lfc$sig
peak_lfc_matrix[rownames(mouse_lfc),"mouse_lfc"] <- mouse_lfc$avg_log2FC
peak_lfc_matrix[rownames(mouse_lfc),"mouse_sig"] <- mouse_lfc$sig

#try lfc=1, sig=2
peak_lfc_matrix$group <- 'mixed'

temp <- rownames(macaque_lfc)[abs(macaque_lfc$avg_log2FC) > 1 & macaque_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'
temp <- rownames(mouse_lfc)[abs(mouse_lfc$avg_log2FC) > 1 & mouse_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'

peak_lfc_matrix[abs(peak_lfc_matrix$macaque_lfc) <= 1 & 
                  abs(peak_lfc_matrix$mouse_lfc) <= 1 & 
                  countOverlaps(query = Brain_ATAC_peakset$macaque[peak_lfc_matrix$peak],subject = macaque_RG_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$mouse[peak_lfc_matrix$peak],subject = mouse_RG_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$human[peak_lfc_matrix$peak],subject = human_RG_peak) > 0,
                "group"] <- 'conserved'

#try enrich heatmap
group <- readRDS(file = './res/step_31_fig_220525/RG_peak_kmeans_group.rds')
group <- group[group$raw_cluster != 'conserved',]
ii <- group$peak_name
group <- group$cluster
names(group) <- ii

ii <- peak_lfc_matrix$peak[peak_lfc_matrix$group == 'conserved']
temp <- rep('species_conserved',times = length(ii))
names(temp) <- ii
group <- append(group,temp)
table(group)
table(duplicated(names(group)))

# enrichheatmap
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/human_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- as.character(Brain_ATAC_peakset$macaque[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$macaque[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/macaque_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- as.character(Brain_ATAC_peakset$mouse[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$mouse[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/mouse_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.004,0.008),colors = c('#FFC30F','#C70039','#581845'))
pdf(file = './res/step_31_fig_220525/RG_peakset_CPM_center_besed_heatmap_conserved_peak_modified.pdf',width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                         'mouse_human_conserved','macaque_specific',
                                                                                         'mouse_macaque_conserved','mouse_specific',
                                                                                         'species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE)
dev.off()

#enrichment distribution
my_insertion_plot <- function(group=group,p1 = p1,p2 = p2,p3 = p3,group.by=group.by){
  temp <- p1@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'mouse'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5','mouse' = '#FDB31A')) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group.by) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_specific')
p_macaque_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_specific')
p_mouse_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_specific')
p_macaque_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_human_conserved')
p_mouse_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_human_conserved')
p_mouse_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

pdf(file = './res/step_31_fig_220525/RG_peakset_mean_insertion_dot_line_plot_conserved_peak_modified.pdf',width = 18,height = 7)
p_human_specific+p_macaque_specific+p_mouse_specific+p_macaque_human_conserved+
  p_mouse_human_conserved+p_mouse_macaque_conserved+p_species_conserved + plot_layout(ncol = 4)
dev.off()

#sankey plot show relationship between species contribution analysis and kmeans analysis
peak_lfc_matrix$kmeans_cluster <- 'dropped'
peak_lfc_matrix[names(group),"kmeans_cluster"] <- as.character(group)

temp <- base::lapply(X = peak_lfc_matrix$peak,FUN = function(x){
  contribution_score <- 1*sum(x %in% mouse_contributed_peak) + 10*sum(x %in% macaque_contributed_peak) + 100*sum(x %in% human_contributed_peak)
  if(contribution_score == 1){
    return('mouse_specific')
  }else if(contribution_score == 10){
    return('macaque_specific')
  }else if(contribution_score == 100){
    return('human_specific')
  }else if(contribution_score == 11){
    return('mouse_macaque_conserved')
  }else if(contribution_score == 101){
    return('mouse_human_conserved')
  }else if(contribution_score == 110){
    return('macaque_human_conserved')
  }else if(contribution_score == 111){
    return('species_conserved')
  }else{
    return('error!')
  }
})
temp <- unlist(temp)
peak_lfc_matrix$peak_cluster <- temp

#scibet
pdf(file = './res/step_31_fig_220525/kmeans_cluster_and_peak_only_cluster_confusion_heatmap.pdf',width = 6,height = 6)
scibet::Confusion_heatmap(ori = peak_lfc_matrix$peak_cluster,prd = peak_lfc_matrix$kmeans_cluster) + 
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 16,face = 'bold')) + 
  ylab('kmeans cluster') + xlab('peak only cluster')
dev.off()

#sankey plot
confusion_matrix <- my_confusion_matrix(ori = paste(peak_lfc_matrix$peak_cluster,'peak_only',sep = '_'),prd = paste(peak_lfc_matrix$kmeans_cluster,'kmeans',sep = '_'))
nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

saveWidget(p1,file = './res/step_31_fig_220525/kmeans_cluster_and_peak_only_cluster_sankey_plot.html')
