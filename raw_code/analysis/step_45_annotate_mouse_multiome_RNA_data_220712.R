#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: annotate mouse multiome RNA data                                ##
## Data: 2022.07.12                                                                ##
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
library(clustree)
library(ggvenn)
library(harmony)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_meta_data.rds')

#re create seurat object
mouse_multiome_Seurat <- mouse_multiome_Seurat@assays$RNA@counts[,rownames(meta_data)]
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat$donor <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"donor"]
mouse_multiome_Seurat$Age <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"Age"]

# re-process mouse multiome Seurat ----------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = c('nCount_RNA','donor','Age'),npcs = 50,preprocess = TRUE)
ndims <- 16
mouse_multiome_Seurat <- FindNeighbors(object = mouse_multiome_Seurat,reduction = 'pca',dims = 1:ndims,assay = 'RNA')
mouse_multiome_Seurat <- FindClusters(object = mouse_multiome_Seurat,resolution = c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5))
mouse_multiome_Seurat <- RunUMAP(object = mouse_multiome_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE,reduction = 'umap') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))

pdf(file = './res/step_45_fig_220712/macaque_multiome_Seurat_cluster_by_different_resolution.pdf',width = 18,height = 8)
clustree(x = mouse_multiome_Seurat@meta.data,prefix = 'RNA_snn_res.')
dev.off()

#the cluster seems not that robust

# label transfer from mouse RNA Seurat ------------------------------------
mouse_RNA_Seurat <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/mouse_RNA_seurat.rds')

gene_list <- dplyr::intersect(x = VariableFeatures(mouse_RNA_Seurat),y = rownames(mouse_multiome_Seurat@assays$RNA@counts))
mouse_RNA_Seurat <- my_process_seurat(object = mouse_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','nFeature_RNA','percent_mito','CC_Difference'),npcs = 50,preprocess = TRUE)
mouse_RNA_Seurat <- my_process_seurat(object = mouse_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = 1,group.by = 'New_cellType',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$RNA@scale.data,object = mouse_RNA_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')

anchors <- my_FindTransferAnchors(reference = mouse_RNA_Seurat,query = mouse_multiome_Seurat,ref_reduction = 'pca',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'RNA',l2.norm = FALSE,dims = 1:35,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = mouse_RNA_Seurat$New_cellType,l2.norm = FALSE,dims = 1:35,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$mouse_cell_type <- mouse_multiome_Seurat$predicted.id

p1 <- DimPlot(object = mouse_RNA_Seurat,reduction = 'umap',group.by = 'New_cellType',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'projected_UMAP',group.by = 'mouse_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
scibet::Confusion_heatmap(ori = mouse_multiome_Seurat$RNA_snn_res.0.5,prd = mouse_multiome_Seurat$mouse_cell_type)

#save meta_data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_45_fig_220712/mouse_multiome_Seurat_label_transfered_from_mouse_RNA_Seurat_meta_data.rds')

# label transfer from macaque multiome data -------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_multiome_Seurat_220712.rds')

temp <- mouse_multiome_Seurat@assays$RNA@counts
mouse_to_macaque_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_Mmul10.csv')
mouse_to_macaque_anno <- mouse_to_macaque_anno[,c(2,1,4,3)]
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_macaque_anno,filter_anno = TRUE,workers = 6,future.globals.maxSize = 20*(1034^3))
mouse_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

macaque_multiome_Seurat[['converted']] <- CreateAssayObject(counts = macaque_multiome_Seurat@assays$RNA@counts,min.cells = 0,min.features = 0)
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$cell_type %in% c('InCGE','OPC','Mic'))]
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(x = rownames(mouse_multiome_Seurat@assays$converted@counts),y = VariableFeatures(macaque_multiome_Seurat))

macaque_multiome_Seurat$species <- 'macaque'
mouse_multiome_Seurat$species <- 'mouse'
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque = macaque_multiome_Seurat,mouse = mouse_multiome_Seurat),assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque = c('nCount_converted','donor'),mouse = NULL),npcs = 50,reference_loading = 'macaque',integration_var = 'species',
                                           harmony_input_dim = 30,max.iter.harmony = 50,UMAP_dim = 30,resolution = 1)

temp_macaque <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'macaque']
temp_macaque$cell_type <- macaque_multiome_Seurat@meta.data[colnames(temp_macaque),"cell_type"]
temp_mouse <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'mouse']
anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_mouse,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:30,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:30,verbose = TRUE)
temp_mouse <- AddMetaData(object = temp_mouse,metadata = predictions)
mouse_multiome_Seurat$macaque_cell_type <- temp_mouse@meta.data[colnames(mouse_multiome_Seurat),"predicted.id"]
mouse_multiome_Seurat$macaque_cell_type_score <- temp_mouse@meta.data[colnames(mouse_multiome_Seurat),"prediction.score.max"]

DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(ArchRPalettes$stallion))


# label transfer from mouse RNA seurat ------------------------------------
mouse_RNA_Seurat <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scRNA_seq/processed_data/mouse_RNA_Seurat_pcaproject_label_transfer_Seuart_220321.rds')
length(VariableFeatures(mouse_RNA_Seurat))
table(VariableFeatures(mouse_RNA_Seurat) %in% rownames(mouse_multiome_Seurat@assays$RNA@counts))
gene_list <- dplyr::intersect(x = VariableFeatures(mouse_RNA_Seurat),y = rownames(mouse_multiome_Seurat@assays$RNA@counts))

mouse_multiome_Seurat$dataset <- 'own'
mouse_RNA_Seurat$dataset <- 'public'
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(own = mouse_multiome_Seurat[,mouse_multiome_Seurat$macaque_cell_type != 'InMGE'],public = mouse_RNA_Seurat),assay = 'RNA',variable_feature = gene_list,
                                           var_to_regress_list = list(own = NULL,public = c('nCount_RNA','batch')),npcs = 50,reference_loading = 'public',integration_var = 'dataset',
                                           harmony_input_dim = 20,max.iter.harmony = 50,UMAP_dim = 20,resolution = 1)

temp_own <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'own']
temp_public <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'public']
temp_public$cell_type <- mouse_RNA_Seurat@meta.data[colnames(temp_public),"cell_type"]

anchors <- my_FindTransferAnchors(reference = temp_public,query = temp_own,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:20,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_public$cell_type,l2.norm = FALSE,dims = 1:20,verbose = TRUE)
temp_own <- AddMetaData(object = temp_own,metadata = predictions)

mouse_multiome_Seurat$mouse_cell_type <- NA
mouse_multiome_Seurat$mouse_cell_type_score <- NA
mouse_multiome_Seurat@meta.data[colnames(temp_own),"mouse_cell_type"] <- temp_own$predicted.id
mouse_multiome_Seurat@meta.data[colnames(temp_own),"mouse_cell_type_score"] <- temp_own$prediction.score.max

p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'mouse_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(ArchRPalettes$stallion)) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(ArchRPalettes$stallion)) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

saveRDS(object = mouse_multiome_Seurat,file = '/data/User/sunym/trash/mouse_multiome_Seurat.rds')

#other plot
Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(temp_own),"cell_type"] <- temp_own$predicted.id
Brain_RNA_Seurat@meta.data[colnames(temp_public),"cell_type"] <- temp_public$cell_type
p1 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'dataset',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

p1 <- DimPlot(object = temp_macaque,group.by = 'cell_type',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = temp_mouse,group.by = 'predicted.id',label = TRUE,repel = TRUE)