#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: process mouse multiome RNAseq data                              ##
## Data: 2022.07.07                                                                ##
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

#library
library(Rmisc)
library(Seurat)
library(SingleCellExperiment)
library(sctransform)
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
library(patchwork)
library(networkD3)
library(htmlwidgets)
library(circlize)
library(harmony)
library(ArchR)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# load data ---------------------------------------------------------------
E145_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix/')
E155_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix/')
E155_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix/')

E145_1 <- E145_1$`Gene Expression`
E155_1 <- E155_1$`Gene Expression`
E155_2 <- E155_2$`Gene Expression`

#rename express matrix
colnames(E145_1) <- paste('E145_1#',colnames(E145_1),sep = '')
colnames(E155_1) <- paste('E155_1#',colnames(E155_1),sep = '')
colnames(E155_2) <- paste('E155_2#',colnames(E155_2),sep = '')

E145_1 <- CreateSeuratObject(counts = E145_1,project = 'E145_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_1 <- CreateSeuratObject(counts = E155_1,project = 'E155_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_2 <- CreateSeuratObject(counts = E155_2,project = 'E155_2',assay = 'RNA',min.cells = 0,min.features = 0)

# QC ----------------------------------------------------------------------
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
E145_1[['percent.mt']] <- PercentageFeatureSet(object = E145_1,features = MT_gene_list)
E155_1[['percent.mt']] <- PercentageFeatureSet(object = E155_1,features = MT_gene_list)
E155_2[['percent.mt']] <- PercentageFeatureSet(object = E155_2,features = MT_gene_list)

#E145_1
VlnPlot(E145_1,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E145_1,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E145_1 <- subset(E145_1,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200)

#E155_1
VlnPlot(E155_1,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_1,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_1 <- subset(E155_1,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200)

#E155_2
VlnPlot(E155_2,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_2 <- subset(E155_2,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200)

# merge -------------------------------------------------------------------
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

VlnPlot(mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(mouse_multiome_Seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")

#add meta data
temp <- unlist(base::lapply(X = colnames(mouse_multiome_Seurat),FUN = function(x){
  temp <- strsplit(x = x,split = '#',fixed = TRUE)
  temp <- temp[[1]][1]
  return(temp)
}))
mouse_multiome_Seurat$donor <- temp

temp <- unlist(base::lapply(X = mouse_multiome_Seurat$donor,FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  temp <- temp[[1]][1]
  return(temp)
}))
mouse_multiome_Seurat$Age <- temp

# process mouse multiome data ---------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
ndims <- 16
mouse_multiome_Seurat <- FindNeighbors(object = mouse_multiome_Seurat,reduction = 'pca',dims = 1:ndims,assay = 'RNA',verbose = TRUE)
mouse_multiome_Seurat <- FindClusters(object = mouse_multiome_Seurat,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),verbose = TRUE)
mouse_multiome_Seurat <- RunUMAP(object = mouse_multiome_Seurat,dims = 1:ndims,reduction = 'pca')

p1 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'donor',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

saveRDS(mouse_multiome_Seurat,file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')

# human region analysis ---------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
human_Whole_Brain <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_whole_brain_RNA_Seurat_nFeature_3000_dim_50_220709.rds')

#convert gene id
mouse_to_human <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
temp <- mouse_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)
mouse_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

gene_list <- dplyr::intersect(x = VariableFeatures(human_Whole_Brain),y = rownames(mouse_multiome_Seurat@assays$converted@counts))
human_Whole_Brain <- my_process_seurat(object = human_Whole_Brain,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA'),npcs = 50,preprocess = TRUE)
human_Whole_Brain <- my_process_seurat(object = human_Whole_Brain,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 50,resolution = 1,group.by = 'cell.type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#project
temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$converted@scale.data,object = human_Whole_Brain,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

#label transfer
anchors <- my_FindTransferAnchors(reference = human_Whole_Brain,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'converted',l2.norm = TRUE,dims = 1:50,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = human_Whole_Brain$area,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$human_area <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$human_area_score <- mouse_multiome_Seurat$prediction.score.max
predictions <- TransferData(anchorset = anchors,refdata = human_Whole_Brain$cell.type,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$human_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$human_cell_type_score <- mouse_multiome_Seurat$prediction.score.max
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_42_fig_220707/mouse_multiome_data_human_whole_brain_area_prediction_220710.rds')

# mouse region analysis ---------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
mouse_nerve_system <- readRDS(file = './data/public/Molecular_Architecture_of_the_Mouse_Nervous_System/mouse_nerve_system_RNA_Seurat_nFeature_3000_dim_50_220709.rds')
mouse_nerve_system <- mouse_nerve_system[,!(mouse_nerve_system$Developmental_compartment %in% c('Diencephalon,Mesencephalon,Rhombencephalon,Spinal cord','Neural crest','Spinal cord'))]

mouse_nerve_system <- FindVariableFeatures(object = mouse_nerve_system,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- dplyr::intersect(x = VariableFeatures(mouse_nerve_system),y = rownames(mouse_multiome_Seurat@assays$RNA@counts))
mouse_nerve_system <- my_process_seurat(object = mouse_nerve_system,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA'),npcs = 50,preprocess = TRUE)
mouse_nerve_system <- my_process_seurat(object = mouse_nerve_system,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 50,resolution = 1,group.by = 'Taxonomy_group',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#project
temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$RNA@scale.data,object = mouse_nerve_system,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')

#label transfer
anchors <- my_FindTransferAnchors(reference = mouse_nerve_system,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'RNA',l2.norm = TRUE,dims = 1:50,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = mouse_nerve_system$Region,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$mouse_area <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$mouse_area_score <- mouse_multiome_Seurat$prediction.score.max
predictions <- TransferData(anchorset = anchors,refdata = mouse_nerve_system$Taxonomy_group,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$mouse_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$mouse_cell_type_score <- mouse_multiome_Seurat$prediction.score.max
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_42_fig_220707/mouse_multiome_data_mouse_nerve_system_area_prediction_220710.rds')

# label transfer ----------------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_Seurat_reannotate_round_1_meta_data.rds')
macaque_multiome_Seurat$new_cell_type <- meta_data[colnames(macaque_multiome_Seurat),"cell_type"]

temp <- macaque_multiome_Seurat@assays$RNA@counts
macaque_to_mouse_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCm39.csv')
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_mouse_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)
macaque_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(mouse_multiome_Seurat@assays$RNA@counts))
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_converted','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'new_cell_type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$RNA@scale.data,object = macaque_multiome_Seurat,assayUsed = 'converted',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'converted',query_assay = 'RNA',l2.norm = TRUE,dims = 1:30,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$new_cell_type,l2.norm = TRUE,dims = 1:30,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$cell_type <- mouse_multiome_Seurat$predicted.id

p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'new_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_42_fig_220707/mouse_multiome_data_macaque_multiome_data_label_transfer_220710.rds')

# mouse data label transfer -----------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
mouse_RNA_Seurat <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/mouse_RNA_seurat.rds')
DimPlot(mouse_RNA_Seurat,group.by = 'New_cellType',label = TRUE,repel = TRUE,reduction = 'umap')

gene_list <- dplyr::intersect(x = VariableFeatures(mouse_RNA_Seurat),y = rownames(mouse_multiome_Seurat@assays$RNA@counts))
mouse_RNA_Seurat <- my_process_seurat(object = mouse_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','nFeature_RNA','percent_mito','CC_Difference'),npcs = 50,preprocess = TRUE)
mouse_RNA_Seurat <- my_process_seurat(object = mouse_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = 1,group.by = 'New_cellType',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$RNA@scale.data,object = mouse_RNA_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')

anchors <- my_FindTransferAnchors(reference = mouse_RNA_Seurat,query = mouse_multiome_Seurat,ref_reduction = 'pca',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'RNA',l2.norm = TRUE,dims = 1:35,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = mouse_RNA_Seurat$New_cellType,l2.norm = TRUE,dims = 1:35,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$mouse_cell_type <- mouse_multiome_Seurat$predicted.id

p1 <- DimPlot(object = mouse_RNA_Seurat,reduction = 'umap',group.by = 'New_cellType',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'projected_UMAP',group.by = 'mouse_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_42_fig_220707/mouse_multiome_data_mouse_RNA_data_label_transfer_220710.rds')

# doublet filter ----------------------------------------------------------

#load annotation
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_human_whole_brain_area_prediction_220710.rds')
mouse_multiome_Seurat$human_region <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"human_area"]
mouse_multiome_Seurat$human_region_score <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"human_area_score"]
mouse_multiome_Seurat$human_cell_type <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"human_cell_type"]
mouse_multiome_Seurat$human_cell_type_score <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"human_cell_type_score"]
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_macaque_multiome_data_label_transfer_220710.rds')
mouse_multiome_Seurat$macaque_cell_type <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"cell_type"]
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_mouse_nerve_system_area_prediction_220710.rds')
mouse_multiome_Seurat$mouse_region <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"mouse_area"]
mouse_multiome_Seurat$mouse_region_score <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"mouse_area_score"]
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_mouse_RNA_data_label_transfer_220710.rds')
mouse_multiome_Seurat$mouse_cell_type <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"mouse_cell_type"]
mouse_multiome_Seurat$mouse_cell_type_score <- meta_data[rownames(mouse_multiome_Seurat@meta.data),"prediction.score.max"]

#load doublet score
E145_1_doublet <- read.csv(file = './res/step_42_fig_220707/E145_1_scrublet.csv')
E155_1_doublet <- read.csv(file = './res/step_42_fig_220707/E155_1_scrublet.csv')
E155_2_doublet <- read.csv(file = './res/step_42_fig_220707/E155_2_scrublet.csv')

E145_1_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E145_1_doublet$barcodes,fixed = TRUE)
E155_1_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_1_doublet$barcodes,fixed = TRUE)
E155_2_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_2_doublet$barcodes,fixed = TRUE)

E145_1_doublet$barcodes <- paste('E145_1#',E145_1_doublet$barcodes,sep = '')
E155_1_doublet$barcodes <- paste('E155_1#',E155_1_doublet$barcodes,sep = '')
E155_2_doublet$barcodes <- paste('E155_2#',E155_2_doublet$barcodes,sep = '')

mouse_multiome_Seurat$doublet_score <- NA
mouse_multiome_Seurat$doublet <- NA

mouse_multiome_Seurat@meta.data[E145_1_doublet$barcodes,"doublet_score"] <- E145_1_doublet$score
mouse_multiome_Seurat@meta.data[E155_1_doublet$barcodes,"doublet_score"] <- E155_1_doublet$score
mouse_multiome_Seurat@meta.data[E155_2_doublet$barcodes,"doublet_score"] <- E155_2_doublet$score

mouse_multiome_Seurat@meta.data[E145_1_doublet$barcodes,"doublet"] <- E145_1_doublet$prediction
mouse_multiome_Seurat@meta.data[E155_1_doublet$barcodes,"doublet"] <- E155_1_doublet$prediction
mouse_multiome_Seurat@meta.data[E155_2_doublet$barcodes,"doublet"] <- E155_2_doublet$prediction

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('doublet_score')) + theme(aspect.ratio = 1)

pdf(file = './res/step_42_fig_220707/mouse_multiome_data_cluster_13_seems_doublet.pdf',width = 18,height = 6)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

pdf(file = './res/step_42_fig_220707/mouse_multiome_data_QC_vlnplot.pdf',width = 8,height = 4)
VlnPlot(object = mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),pt.size = 0.1,group.by = 'donor')
dev.off()

#cluster 13 seems doublet
pdf(file = './res/step_42_fig_220707/mouse_multiome_data_cluster_13_seems_doublet_vlnplot.pdf',width = 12,height = 6)
VlnPlot(object = mouse_multiome_Seurat,features = c('doublet_score'),pt.size = 0,group.by = 'RNA_snn_res.0.5')
dev.off()

mouse_multiome_Seurat <- mouse_multiome_Seurat[,which(!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('13')))]

#nFeatures nCounts filter
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('nCount_RNA')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('nFeature_RNA')) + theme(aspect.ratio = 1)
p1+p2+p3+p4+plot_layout(ncol = 4)

# filter weird clusters ---------------------------------------------------
#filter Meis2 expressed In
table(mouse_multiome_Seurat[,mouse_multiome_Seurat$RNA_snn_res.0.5 == '1']$macaque_cell_type)

pdf(file = './res/step_42_fig_220707/mouse_multiome_data_Meis2_expression_featureplot.pdf',width = 16,height = 10)
FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Adarb2','Meis2','Neurod6'),pt.size = 0.1)
dev.off()

p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'human_region',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
pdf(file = './res/step_42_fig_220707/mouse_cluster_1_12_region_predict_dimplot.pdf',width = 14,height = 7)
p1+p2+plot_layout(ncol = 2)
dev.off()

#filter cluster 1 and cluster 12 and cluster 10
mouse_multiome_Seurat <- mouse_multiome_Seurat[,!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('1','12','10'))]

#redo process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
ndims <- 16
mouse_multiome_Seurat <- FindNeighbors(object = mouse_multiome_Seurat,reduction = 'pca',dims = 1:ndims,assay = 'RNA',verbose = TRUE)
mouse_multiome_Seurat <- FindClusters(object = mouse_multiome_Seurat,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),verbose = TRUE)
mouse_multiome_Seurat <- RunUMAP(object = mouse_multiome_Seurat,dims = 1:ndims,reduction = 'pca')

p1 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
pdf(file = './res/step_42_fig_220707/mouse_multiome_data_redo_umap_round_1.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#InCGE and InMGE
pdf(file = './res/step_42_fig_220707/mouse_multiome_data_InMGE_marker_expression_featureplot.pdf',width = 6,height = 8)
FeaturePlot(object = mouse_multiome_Seurat,features = c('Sst','Npy','Lhx6','Nxph1','Nxph2'))
dev.off()
#cluster 2 is InMGE

pdf(file = './res/step_42_fig_220707/mouse_multiome_data_InCGE_marker_expression_featureplot.pdf',width = 6,height = 6)
FeaturePlot(object = mouse_multiome_Seurat,features = c('Pax6','Sp8','Cxcl14','Htr3a'))
dev.off()
#InCGE do not exist yet.

#cluster 7 and 9 marker
temp <- as.character(unique(mouse_multiome_Seurat$RNA_snn_res.0.5))
temp <- temp[temp != 7]
cluster_7_marker <- base::lapply(X = temp,FUN = function(x){
  temp <- FindMarkers(object = mouse_multiome_Seurat,ident.1 = '7',ident.2 = x,group.by = 'RNA_snn_res.0.5',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
  return(rownames(temp))
})
cluster_7_marker <- Reduce(f = dplyr::intersect,x = cluster_7_marker)

VlnPlot(object = mouse_multiome_Seurat,features = cluster_7_marker[109:125],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')
#Hmgb2 Ptov1 Sf3b2 Ran

pdf(file = './res/step_42_fig_220707/mouse_multiome_data_cluster_7_marker_featureplot.pdf',width = 6,height = 6)
FeaturePlot(object = mouse_multiome_Seurat,features = c('Hmgb2','Ptov1','Sf3b2','Ran'))
dev.off()

#cluster 7 can be filter
temp <- as.character(unique(mouse_multiome_Seurat$RNA_snn_res.0.5))
temp <- temp[temp != 9]
cluster_9_marker <- base::lapply(X = temp,FUN = function(x){
  temp <- FindMarkers(object = mouse_multiome_Seurat,ident.1 = '9',ident.2 = x,group.by = 'RNA_snn_res.0.5',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
  return(rownames(temp))
})
cluster_9_marker <- Reduce(f = dplyr::intersect,x = cluster_9_marker)

VlnPlot(object = mouse_multiome_Seurat,features = cluster_9_marker[85:98],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')
#Pbx3 Zfhx3 Gm44593 Meis1 lsl1 C79798 Sorcs3 Gm27032 Syt6 Plcxd3 Gm6994 Esrrg Grm7

FeaturePlot(object = mouse_multiome_Seurat,features = c('Pbx3','Zfhx3','Gm44593','Meis1','Isl1','C79798','Sorcs3','Gm27032','Syt6','Plcxd3','Gm6994','Esrrg','Grm7'))
#Zfhx3 Meis1 Syt6 checked
pdf(file = './res/step_42_fig_220707/mouse_multiome_data_cluster_9_marker_featureplot.pdf',width = 12,height = 6)
geneSetAveragePlot(genes = c('Zfhx3','Meis1','Syt6'),object = mouse_multiome_Seurat,object.class = 'Seurat',assay = 'RNA',embedding = 'umap',
                   plot.title = 'cluster 9 marker',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],
                   aspectratio = 1,point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)
dev.off()

pdf(file = './res/step_42_fig_220707/mouse_cluster_9_possible_region.pdf',width = 8,height = 8)
DimPlot(object = mouse_multiome_Seurat,group.by = 'human_region',label = TRUE,repel = TRUE)
dev.off()

#filter cluster 7 and 9
mouse_multiome_Seurat <- mouse_multiome_Seurat[,!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('7','9'))]
DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE)

#redo process round 2
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
ndims <- 12
mouse_multiome_Seurat <- FindNeighbors(object = mouse_multiome_Seurat,reduction = 'pca',dims = 1:ndims,assay = 'RNA',verbose = TRUE)
mouse_multiome_Seurat <- FindClusters(object = mouse_multiome_Seurat,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),verbose = TRUE)
mouse_multiome_Seurat <- RunUMAP(object = mouse_multiome_Seurat,dims = 1:ndims,reduction = 'pca')

p1 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#annotate
#Cycling marker
FeaturePlot(object = mouse_multiome_Seurat,features = c('Top2a','Mki67','Clspn','Aurka'))
#cluster 8 Cycling

#RG marker
FeaturePlot(object = mouse_multiome_Seurat,features = c('Sox2','Pax6','Hes5'))
#not that specific
#cluster 5 10 RG

#IP marker
FeaturePlot(object = mouse_multiome_Seurat,features = c('Eomes','Ppp1r17'))
#cluster 6 IP

#Per marker
FeaturePlot(object = mouse_multiome_Seurat,features = c('Foxc2','Pdgfrb'))
#cluster 12 Per

#End marker
FeaturePlot(object = mouse_multiome_Seurat,features = c('Cldn5','Pecam1'))
#cluster 11 End

#In marker
FeaturePlot(object = mouse_multiome_Seurat,features = c('Gad1','Gad2','Lhx6','Sst','Sp8','Nr2f2'))
#cluster 2 InMGE

#cluster 3 4 0 1 Ex
mouse_multiome_Seurat$cell_type <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('8'),"cell_type"] <- 'Cycling'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('5','10'),"cell_type"] <- 'RG'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('6'),"cell_type"] <- 'IP'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('12'),"cell_type"] <- 'Per'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('11'),"cell_type"] <- 'End'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('2'),"cell_type"] <- 'InMGE'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('3','4','0','1'),"cell_type"] <- 'Ex'

p1 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
pdf(file = './res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_dimplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(meta_data,file = './res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_meta_data.rds')

