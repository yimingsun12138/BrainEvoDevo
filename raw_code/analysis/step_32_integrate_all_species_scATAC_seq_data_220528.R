#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: integrate all species scATAC_seq data                           ##
## Data: 2022.05.28                                                                ##
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

# use macaque parameter 0.1 mouse parameter 0.4 ---------------------------
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

## macaque liftover to human -----------------------------------------------
lifted_macaque_peakset <- my_unique_peakset_liftover(ori_GRanges = macaque_peakset,
                                                     UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                     chain_file = './data/reference/rheMac10ToHg38.over.chain',
                                                     liftOver_mismatch = 0.1,length_filter = TRUE,length_mismatch = 0.1,
                                                     chr_filter = TRUE,mapped_chr = unique(as.character(human_peakset@seqnames)),
                                                     overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220528')

my_send_sms('liftover done!')

#merge with human peakset
macaque_human_merged_peakset <- my_bedtools_merge(peakset_x = lifted_macaque_peakset$mapped,
                                                  peakset_y = human_peakset,
                                                  bedtools_path = '/data/User/sunym/software/bedtools/bedtools',
                                                  d = 0,bedtools_param = NULL,
                                                  tmp_path = '/data/User/sunym/temp/tmp_220528')

#re-liftover to macaque
macaque_human_merged_peakset <- my_unique_peakset_liftover(ori_GRanges = macaque_human_merged_peakset,
                                                           UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                           chain_file = './data/reference/hg38ToRheMac10.over.chain',
                                                           liftOver_mismatch = 0.1,length_filter = TRUE,length_mismatch = 0.1,
                                                           chr_filter = TRUE,mapped_chr = unique(as.character(macaque_peakset@seqnames)),
                                                           overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220528')

my_send_sms('liftover done!')

## mouse liftover to human -------------------------------------------------
lifted_mouse_peakset <- my_unique_peakset_liftover(ori_GRanges = mouse_peakset,
                                                   UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                   chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/mm10ToHg38.over.chain',
                                                   liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                   chr_filter = TRUE,mapped_chr = unique(as.character(human_peakset@seqnames)),
                                                   overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220528')

my_send_sms('liftover done!')

#merge with macaque_human_merged_peakset
Brain_ATAC_peakset <- my_bedtools_merge(peakset_x = lifted_mouse_peakset$mapped,
                                        peakset_y = macaque_human_merged_peakset$ori,
                                        bedtools_path = '/data/User/sunym/software/bedtools/bedtools',
                                        d = 0,bedtools_param = NULL,
                                        tmp_path = '/data/User/sunym/temp/tmp_220528')

## create species unique peakset -------------------------------------------

#re-liftover to mouse
Brain_ATAC_peakset <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset,
                                                 UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                 chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/hg38ToMm10.over.chain',
                                                 liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                 chr_filter = TRUE,mapped_chr = unique(as.character(mouse_peakset@seqnames)),
                                                 overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220528')

my_send_sms('liftover done!')

#re-liftover to macaque
Brain_ATAC_peakset_macaque <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset$ori,
                                                         UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                         chain_file = './data/reference/hg38ToRheMac10.over.chain',
                                                         liftOver_mismatch = 0.1,length_filter = TRUE,length_mismatch = 0.1,
                                                         chr_filter = TRUE,mapped_chr = unique(as.character(macaque_peakset@seqnames)),
                                                         overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220528')

my_send_sms('liftover done!')

#re-liftover to mouse
Brain_ATAC_peakset_mouse <- my_unique_peakset_liftover(ori_GRanges = Brain_ATAC_peakset_macaque$ori,
                                                       UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                       chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/hg38ToMm10.over.chain',
                                                       liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                       chr_filter = TRUE,mapped_chr = unique(as.character(mouse_peakset@seqnames)),
                                                       overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220528')

my_send_sms('liftover done!')

#validate
table(names(Brain_ATAC_peakset_macaque$ori) == names(Brain_ATAC_peakset_mouse$ori))
Brain_ATAC_peakset <- SimpleList(
  human = Brain_ATAC_peakset_macaque$ori,
  macaque = Brain_ATAC_peakset_macaque$mapped,
  mouse = Brain_ATAC_peakset_mouse$mapped
)
saveRDS(Brain_ATAC_peakset,file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')

# create species peak matrix ----------------------------------------------
#all run on 202.205.131.32, ArchR is fucking shit.

## human peak matrix -------------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')

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
saveRDS(human_peak_matrix,file = './res/step_32_fig_220528/human_peak_matrix.rds')

## macaque peak matrix -----------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')

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
saveRDS(macaque_peak_matrix,file = './res/step_32_fig_220528/macaque_peak_matrix.rds')

## mouse peak matrix -------------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')

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
saveRDS(mouse_peak_matrix,file = './res/step_32_fig_220528/mouse_peak_matrix.rds')

# create Brain peak matrix ------------------------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')
human_peak_matrix <- readRDS(file = './res/step_32_fig_220528/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_32_fig_220528/macaque_peak_matrix.rds')
mouse_peak_matrix <- readRDS(file = './res/step_32_fig_220528/mouse_peak_matrix.rds')

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
saveRDS(Brain_peak_matrix,file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')

# sample distance analysis ------------------------------------------------
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')

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
Brain_peak_matrix <- RunUMAP(object = Brain_peak_matrix,dims = 1:10,reduction = 'pca')
DimPlot(object = Brain_peak_matrix,group.by = 'cell_type',label = TRUE,repel = TRUE)

temp <- Brain_peak_matrix@meta.data
temp <- cbind(temp,Brain_peak_matrix@reductions$umap@cell.embeddings)

pdf(file = './res/step_32_fig_220528/Brain_peak_matrix_umap_dimplot.pdf',width = 6,height = 5)
ggplot(temp,aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'UMAP plot')
dev.off()

#cluster analysis
hc <- hclust(d = dist(t(Brain_peak_matrix@assays$RNA@scale.data[VariableFeatures(Brain_peak_matrix),])),method = 'ward.D')

pdf(file = './res/step_32_fig_220528/Brain_peak_matrix_hclust.pdf',width = 12,height = 6)
plot(hc)
dev.off()

#Brain peak matrix create done!

#notice that the pseudo bulk lose the meta data of cell composition, re-generate it in step_32_220603