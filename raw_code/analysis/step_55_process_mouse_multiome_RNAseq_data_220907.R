#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: process mouse multiome RNAseq data                              ##
## Data: 2022.09.07                                                                ##
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
E145_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E145_2/outs/filtered_feature_bc_matrix/')
E155_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix/')
E155_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix/')
E155_3 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E155_3/outs/filtered_feature_bc_matrix/')
E165_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E165_2/outs/filtered_feature_bc_matrix/')

E145_1 <- E145_1$`Gene Expression`
E145_2 <- E145_2$`Gene Expression`
E155_1 <- E155_1$`Gene Expression`
E155_2 <- E155_2$`Gene Expression`
E155_3 <- E155_3$`Gene Expression`
E165_2 <- E165_2$`Gene Expression`

#rename express matrix
colnames(E145_1) <- paste('E145_1#',colnames(E145_1),sep = '')
colnames(E145_2) <- paste('E145_2#',colnames(E145_2),sep = '')
colnames(E155_1) <- paste('E155_1#',colnames(E155_1),sep = '')
colnames(E155_2) <- paste('E155_2#',colnames(E155_2),sep = '')
colnames(E155_3) <- paste('E155_3#',colnames(E155_3),sep = '')
colnames(E165_2) <- paste('E165_2#',colnames(E165_2),sep = '')

E145_1 <- CreateSeuratObject(counts = E145_1,project = 'E145_1',assay = 'RNA',min.cells = 0,min.features = 0)
E145_2 <- CreateSeuratObject(counts = E145_2,project = 'E145_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_1 <- CreateSeuratObject(counts = E155_1,project = 'E155_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_2 <- CreateSeuratObject(counts = E155_2,project = 'E155_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_3 <- CreateSeuratObject(counts = E155_3,project = 'E155_3',assay = 'RNA',min.cells = 0,min.features = 0)
E165_2 <- CreateSeuratObject(counts = E165_2,project = 'E165_2',assay = 'RNA',min.cells = 0,min.features = 0)

# QC ----------------------------------------------------------------------
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
E145_1[['percent.mt']] <- PercentageFeatureSet(object = E145_1,features = MT_gene_list)
E145_2[['percent.mt']] <- PercentageFeatureSet(object = E145_2,features = MT_gene_list)
E155_1[['percent.mt']] <- PercentageFeatureSet(object = E155_1,features = MT_gene_list)
E155_2[['percent.mt']] <- PercentageFeatureSet(object = E155_2,features = MT_gene_list)
E155_3[['percent.mt']] <- PercentageFeatureSet(object = E155_3,features = MT_gene_list)
E165_2[['percent.mt']] <- PercentageFeatureSet(object = E165_2,features = MT_gene_list)

#E145_1
VlnPlot(E145_1,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E145_1,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E145_1 <- subset(E145_1,subset = nFeature_RNA < 7500 & nCount_RNA < 20000 & nFeature_RNA > 200)

#E145_2
VlnPlot(E145_2,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E145_2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E145_2 <- subset(E145_2,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200 & percent.mt < 10)

#E155_1
VlnPlot(E155_1,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_1,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_1 <- subset(E155_1,subset = nFeature_RNA < 7500 & nCount_RNA < 20000 & nFeature_RNA > 200)

#E155_2
VlnPlot(E155_2,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_2 <- subset(E155_2,subset = nFeature_RNA < 7500 & nCount_RNA < 20000 & nFeature_RNA > 200)

#E155_3
VlnPlot(E155_3,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_3,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_3 <- subset(E155_3,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200 & percent.mt < 10)

#E165_2
VlnPlot(E165_2,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E165_2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E165_2 <- subset(E165_2,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200 & percent.mt < 10)

# merge -------------------------------------------------------------------
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
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

mouse_multiome_Seurat$batch <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_1','E155_1','E155_2'),"batch"] <- '220707'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_2','E155_3','E165_2'),"batch"] <- '220907'

# process mouse multiome data ---------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 20
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.7',label = TRUE)
p1 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'donor',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#save temp data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/raw_mouse_multiome_Seurat_meta_data_220907.rds')

# human region analysis ---------------------------------------------------
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

#save data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/human_whole_brain_predicted_mouse_multiome_Seurat_meta_data_220908.rds')

#display
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE)

# label transfer from macaque data ----------------------------------------
#load data
E145_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix/')
E145_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E145_2/outs/filtered_feature_bc_matrix/')
E155_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix/')
E155_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix/')
E155_3 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E155_3/outs/filtered_feature_bc_matrix/')
E165_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E165_2/outs/filtered_feature_bc_matrix/')

E145_1 <- E145_1$`Gene Expression`
E145_2 <- E145_2$`Gene Expression`
E155_1 <- E155_1$`Gene Expression`
E155_2 <- E155_2$`Gene Expression`
E155_3 <- E155_3$`Gene Expression`
E165_2 <- E165_2$`Gene Expression`

colnames(E145_1) <- paste('E145_1#',colnames(E145_1),sep = '')
colnames(E145_2) <- paste('E145_2#',colnames(E145_2),sep = '')
colnames(E155_1) <- paste('E155_1#',colnames(E155_1),sep = '')
colnames(E155_2) <- paste('E155_2#',colnames(E155_2),sep = '')
colnames(E155_3) <- paste('E155_3#',colnames(E155_3),sep = '')
colnames(E165_2) <- paste('E165_2#',colnames(E165_2),sep = '')

E145_1 <- CreateSeuratObject(counts = E145_1,project = 'E145_1',assay = 'RNA',min.cells = 0,min.features = 0)
E145_2 <- CreateSeuratObject(counts = E145_2,project = 'E145_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_1 <- CreateSeuratObject(counts = E155_1,project = 'E155_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_2 <- CreateSeuratObject(counts = E155_2,project = 'E155_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_3 <- CreateSeuratObject(counts = E155_3,project = 'E155_3',assay = 'RNA',min.cells = 0,min.features = 0)
E165_2 <- CreateSeuratObject(counts = E165_2,project = 'E165_2',assay = 'RNA',min.cells = 0,min.features = 0)

MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
E145_1[['percent.mt']] <- PercentageFeatureSet(object = E145_1,features = MT_gene_list)
E145_2[['percent.mt']] <- PercentageFeatureSet(object = E145_2,features = MT_gene_list)
E155_1[['percent.mt']] <- PercentageFeatureSet(object = E155_1,features = MT_gene_list)
E155_2[['percent.mt']] <- PercentageFeatureSet(object = E155_2,features = MT_gene_list)
E155_3[['percent.mt']] <- PercentageFeatureSet(object = E155_3,features = MT_gene_list)
E165_2[['percent.mt']] <- PercentageFeatureSet(object = E165_2,features = MT_gene_list)

VlnPlot(E145_1,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E145_1,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E145_1 <- subset(E145_1,subset = nFeature_RNA < 7500 & nCount_RNA < 20000 & nFeature_RNA > 200)

VlnPlot(E145_2,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E145_2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E145_2 <- subset(E145_2,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200 & percent.mt < 10)

VlnPlot(E155_1,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_1,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_1 <- subset(E155_1,subset = nFeature_RNA < 7500 & nCount_RNA < 20000 & nFeature_RNA > 200)

VlnPlot(E155_2,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_2 <- subset(E155_2,subset = nFeature_RNA < 7500 & nCount_RNA < 20000 & nFeature_RNA > 200)

VlnPlot(E155_3,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E155_3,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E155_3 <- subset(E155_3,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200 & percent.mt < 10)

VlnPlot(E165_2,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(E165_2,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
E165_2 <- subset(E165_2,subset = nFeature_RNA < 7500 & nCount_RNA < 30000 & nFeature_RNA > 200 & percent.mt < 10)

mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

VlnPlot(mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(mouse_multiome_Seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")

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

mouse_multiome_Seurat$batch <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_1','E155_1','E155_2'),"batch"] <- '220707'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_2','E155_3','E165_2'),"batch"] <- '220907'

#pre-process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 20
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.7',label = TRUE)
p1 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,pt.size = 0.1,group.by = 'donor',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#add human region information
meta_data <- readRDS(file = './res/step_55_fig_220907/human_whole_brain_predicted_mouse_multiome_Seurat_meta_data_220908.rds')
mouse_multiome_Seurat$human_area <- meta_data[colnames(mouse_multiome_Seurat),"human_area"]
mouse_multiome_Seurat$human_area_score <- meta_data[colnames(mouse_multiome_Seurat),"human_area_score"]
mouse_multiome_Seurat$human_cell_type <- meta_data[colnames(mouse_multiome_Seurat),"human_cell_type"]
mouse_multiome_Seurat$human_cell_type_score <- meta_data[colnames(mouse_multiome_Seurat),"human_cell_type_score"]

#label transfer from macaque
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')

#convert gene id
mouse_to_human <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
temp <- mouse_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)
mouse_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(mouse_multiome_Seurat@assays$converted@counts))
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 32,resolution = 1,group.by = 'cell_type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#project
temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$converted@scale.data,object = macaque_multiome_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'converted',l2.norm = TRUE,dims = 1:32,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$cell_type,l2.norm = TRUE,dims = 1:32,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$macaque_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$macaque_cell_type_score <- mouse_multiome_Seurat$prediction.score.max

#dimplot
DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE,reduction = 'umap')

#save data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/macaque_cell_type_predicted_mouse_multiome_Seurat_meta_data_220909.rds')

# doublet filter ----------------------------------------------------------
#load annotation
E145_1_doublet <- read.csv(file = './res/step_42_fig_220707/E145_1_scrublet.csv')
E145_2_doublet <- read.csv(file = './res/step_55_fig_220907/E145_2_scrublet.csv')
E155_1_doublet <- read.csv(file = './res/step_42_fig_220707/E155_1_scrublet.csv')
E155_2_doublet <- read.csv(file = './res/step_42_fig_220707/E155_2_scrublet.csv')
E155_3_doublet <- read.csv(file = './res/step_55_fig_220907/E155_3_scrublet.csv')
E165_2_doublet <- read.csv(file = './res/step_55_fig_220907/E165_2_scrublet.csv')

E145_1_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E145_1_doublet$barcodes,fixed = TRUE)
E145_2_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E145_2_doublet$barcodes,fixed = TRUE)
E155_1_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_1_doublet$barcodes,fixed = TRUE)
E155_2_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_2_doublet$barcodes,fixed = TRUE)
E155_3_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_3_doublet$barcodes,fixed = TRUE)
E165_2_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E165_2_doublet$barcodes,fixed = TRUE)

E145_1_doublet$barcodes <- paste('E145_1#',E145_1_doublet$barcodes,sep = '')
E145_2_doublet$barcodes <- paste('E145_2#',E145_2_doublet$barcodes,sep = '')
E155_1_doublet$barcodes <- paste('E155_1#',E155_1_doublet$barcodes,sep = '')
E155_2_doublet$barcodes <- paste('E155_2#',E155_2_doublet$barcodes,sep = '')
E155_3_doublet$barcodes <- paste('E155_3#',E155_3_doublet$barcodes,sep = '')
E165_2_doublet$barcodes <- paste('E165_2#',E165_2_doublet$barcodes,sep = '')

mouse_multiome_Seurat$doublet_score <- NA
mouse_multiome_Seurat$doublet <- NA

mouse_multiome_Seurat@meta.data[E145_1_doublet$barcodes,"doublet_score"] <- E145_1_doublet$score
mouse_multiome_Seurat@meta.data[E145_2_doublet$barcodes,"doublet_score"] <- E145_2_doublet$score
mouse_multiome_Seurat@meta.data[E155_1_doublet$barcodes,"doublet_score"] <- E155_1_doublet$score
mouse_multiome_Seurat@meta.data[E155_2_doublet$barcodes,"doublet_score"] <- E155_2_doublet$score
mouse_multiome_Seurat@meta.data[E155_3_doublet$barcodes,"doublet_score"] <- E155_3_doublet$score
mouse_multiome_Seurat@meta.data[E165_2_doublet$barcodes,"doublet_score"] <- E165_2_doublet$score

mouse_multiome_Seurat@meta.data[E145_1_doublet$barcodes,"doublet"] <- E145_1_doublet$prediction
mouse_multiome_Seurat@meta.data[E145_2_doublet$barcodes,"doublet"] <- E145_2_doublet$prediction
mouse_multiome_Seurat@meta.data[E155_1_doublet$barcodes,"doublet"] <- E155_1_doublet$prediction
mouse_multiome_Seurat@meta.data[E155_2_doublet$barcodes,"doublet"] <- E155_2_doublet$prediction
mouse_multiome_Seurat@meta.data[E155_3_doublet$barcodes,"doublet"] <- E155_3_doublet$prediction
mouse_multiome_Seurat@meta.data[E165_2_doublet$barcodes,"doublet"] <- E165_2_doublet$prediction

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('doublet_score'),reduction = 'umap') + theme(aspect.ratio = 1)

pdf(file = './res/step_55_fig_220907/mouse_multiome_data_cluster_1_seems_doublet.pdf',width = 18,height = 6)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

pdf(file = './res/step_55_fig_220907/mouse_multiome_data_cluster_1_seems_doublet_vlnplot.pdf',width = 12,height = 6)
VlnPlot(object = mouse_multiome_Seurat,features = 'doublet_score',group.by = 'RNA_snn_res.0.5')
dev.off()

#filter cluster 1
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[mouse_multiome_Seurat$RNA_snn_res.0.5 != '1']]
table(mouse_multiome_Seurat$RNA_snn_res.0.5)
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('doublet_score'),reduction = 'umap') + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

# nCount nFeature filter --------------------------------------------------
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('nCount_RNA'),reduction = 'umap') + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('nFeature_RNA'),reduction = 'umap') + theme(aspect.ratio = 1)
p1+p2+p3+p4+plot_layout(ncol = 4)

#not very specific

# filter Meis2 expressed In -----------------------------------------------
DefaultAssay(object = mouse_multiome_Seurat) <- 'RNA'
pdf(file = './res/step_55_fig_220907/mouse_multiome_data_Meis2_expression_featureplot.pdf',width = 16,height = 10)
FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Adarb2','Meis2','Neurod6'),pt.size = 0.1,reduction = 'umap')
dev.off()

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#InCGE ?
FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Meis2','Sp8','Cxcl14','Htr3a'),pt.size = 0.1,reduction = 'umap',order = TRUE)

#filter cluster 14 6
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('14','6'))]]
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

# redo cluster and umap ---------------------------------------------------
ndims <- 15
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.5',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

# filter cells from different areas ---------------------------------------
#InMGE
FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Meis2','Sp8','Cxcl14','Htr3a','Neurod6'),pt.size = 0.1,reduction = 'umap',order = TRUE)
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 0 seems very weird
table(mouse_multiome_Seurat[,mouse_multiome_Seurat$RNA_snn_res.0.5 == '0']$human_area)
VlnPlot(object = mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','doublet_score'),group.by = 'RNA_snn_res.0.5',pt.size = 0)
#mostly from thalamus, not ex, filter it
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('0'))]]

#redo cluster and umap
ndims <- 15
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.5',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#save metadata
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/mouse_multiome_meta_data_220909_1638_temp_data.rds')

# redo human region analysis and macaque cell type prediction -------------
#load data
E145_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix/')
E145_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E145_2/outs/filtered_feature_bc_matrix/')
E155_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix/')
E155_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix/')
E155_3 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E155_3/outs/filtered_feature_bc_matrix/')
E165_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E165_2/outs/filtered_feature_bc_matrix/')

E145_1 <- E145_1$`Gene Expression`
E145_2 <- E145_2$`Gene Expression`
E155_1 <- E155_1$`Gene Expression`
E155_2 <- E155_2$`Gene Expression`
E155_3 <- E155_3$`Gene Expression`
E165_2 <- E165_2$`Gene Expression`

colnames(E145_1) <- paste('E145_1#',colnames(E145_1),sep = '')
colnames(E145_2) <- paste('E145_2#',colnames(E145_2),sep = '')
colnames(E155_1) <- paste('E155_1#',colnames(E155_1),sep = '')
colnames(E155_2) <- paste('E155_2#',colnames(E155_2),sep = '')
colnames(E155_3) <- paste('E155_3#',colnames(E155_3),sep = '')
colnames(E165_2) <- paste('E165_2#',colnames(E165_2),sep = '')

E145_1 <- CreateSeuratObject(counts = E145_1,project = 'E145_1',assay = 'RNA',min.cells = 0,min.features = 0)
E145_2 <- CreateSeuratObject(counts = E145_2,project = 'E145_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_1 <- CreateSeuratObject(counts = E155_1,project = 'E155_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_2 <- CreateSeuratObject(counts = E155_2,project = 'E155_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_3 <- CreateSeuratObject(counts = E155_3,project = 'E155_3',assay = 'RNA',min.cells = 0,min.features = 0)
E165_2 <- CreateSeuratObject(counts = E165_2,project = 'E165_2',assay = 'RNA',min.cells = 0,min.features = 0)

MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
E145_1[['percent.mt']] <- PercentageFeatureSet(object = E145_1,features = MT_gene_list)
E145_2[['percent.mt']] <- PercentageFeatureSet(object = E145_2,features = MT_gene_list)
E155_1[['percent.mt']] <- PercentageFeatureSet(object = E155_1,features = MT_gene_list)
E155_2[['percent.mt']] <- PercentageFeatureSet(object = E155_2,features = MT_gene_list)
E155_3[['percent.mt']] <- PercentageFeatureSet(object = E155_3,features = MT_gene_list)
E165_2[['percent.mt']] <- PercentageFeatureSet(object = E165_2,features = MT_gene_list)

mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

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

mouse_multiome_Seurat$batch <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_1','E155_1','E155_2'),"batch"] <- '220707'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_2','E155_3','E165_2'),"batch"] <- '220907'

#load meta data
meta_data <- readRDS(file = './res/step_55_fig_220907/mouse_multiome_meta_data_220909_1638_temp_data.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]

#save memory
remove(list = c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2'))
gc()

#pre-process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 20
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.5',label = TRUE)

#convert gene id
mouse_to_human <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
temp <- mouse_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)
mouse_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

#human region analysis
human_Whole_Brain <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_whole_brain_RNA_Seurat_nFeature_3000_dim_50_220709.rds')
gene_list <- dplyr::intersect(x = VariableFeatures(human_Whole_Brain),y = rownames(mouse_multiome_Seurat@assays$converted@counts))
human_Whole_Brain <- my_process_seurat(object = human_Whole_Brain,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA'),npcs = 50,preprocess = TRUE)
human_Whole_Brain <- my_process_seurat(object = human_Whole_Brain,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 50,resolution = 1,group.by = 'cell.type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$converted@scale.data,object = human_Whole_Brain,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

anchors <- my_FindTransferAnchors(reference = human_Whole_Brain,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'converted',l2.norm = TRUE,dims = 1:50,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = human_Whole_Brain$area,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$human_area <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$human_area_score <- mouse_multiome_Seurat$prediction.score.max
predictions <- TransferData(anchorset = anchors,refdata = human_Whole_Brain$cell.type,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$human_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$human_cell_type_score <- mouse_multiome_Seurat$prediction.score.max

#save meta data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/human_whole_brain_predicted_mouse_multiome_Seurat_meta_data_220909_1829_temp_data.rds')

#display
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE)

# macaque cell type prediction --------------------------------------------
# load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')
gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(mouse_multiome_Seurat@assays$converted@counts))
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 32,resolution = 1,group.by = 'cell_type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#project
temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$converted@scale.data,object = macaque_multiome_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'converted',l2.norm = TRUE,dims = 1:32,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$cell_type,l2.norm = TRUE,dims = 1:32,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$macaque_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$macaque_cell_type_score <- mouse_multiome_Seurat$prediction.score.max

#display
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)

#save meta data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/macaque_cell_type_predicted_mouse_multiome_Seurat_meta_data_220909_1851_temp_data.rds')

# further filter ----------------------------------------------------------
DefaultAssay(object = mouse_multiome_Seurat) <- 'RNA'

# InMGE and Meis2
FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Meis2','Sp8','Cxcl14','Htr3a','Neurod6'),pt.size = 0.1,reduction = 'umap',order = TRUE)
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 12 19 seems weird
meta_data <- readRDS(file = './res/step_55_fig_220907/mouse_multiome_meta_data_220909_1638_temp_data.rds')
mouse_multiome_Seurat$doublet_score <- meta_data[colnames(mouse_multiome_Seurat),"doublet_score"]
VlnPlot(object = mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','doublet_score'),group.by = 'RNA_snn_res.0.5',pt.size = 0)

#cluster 12 19 marker 
temp <- FindMarkers(object = mouse_multiome_Seurat,ident.1 = c('12','19'),group.by = 'RNA_snn_res.0.5',assay = 'RNA',slot = 'data',only.pos = TRUE,logfc.threshold = 0.5,verbose = TRUE)
FeaturePlot(object = mouse_multiome_Seurat,features = rownames(temp)[1:12])
FeaturePlot(object = mouse_multiome_Seurat,features = c('Lhx2','Lhx9'))

#filter cluster 19 12
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('12','19'))]]
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

# redo clustering and umap ---------------------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 12
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.7',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

# GAD2 expression cycling -------------------------------------------------
FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Meis2','Sp8','Cxcl14','Htr3a','Neurod6'),pt.size = 0.1,reduction = 'umap',order = TRUE)
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

table(mouse_multiome_Seurat[,mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('10')]$macaque_cell_type)
VlnPlot(object = mouse_multiome_Seurat,features = c('Gad2'),pt.size = 0,group.by = 'RNA_snn_res.0.5',assay = 'RNA')

#filter cluster 10
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('10'))]]

# redo clustering and umap ---------------------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 12
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.7',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 12 seems weird
FeaturePlot(object = mouse_multiome_Seurat,features = c('Eomes'),order = TRUE,reduction = 'umap')
#cluster 12 is cycling IP

# cluster 1 seems weird ---------------------------------------------------
#is cluster 1 InCGE?
FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Meis2','Sp8','Cxcl14','Htr3a','Nr2f2','Calb2','Neurod6'),pt.size = 0.1,reduction = 'umap',order = TRUE)
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 1 marker
temp <- FindMarkers(object = mouse_multiome_Seurat,ident.1 = c('1'),group.by = 'RNA_snn_res.0.5',assay = 'RNA',slot = 'data',only.pos = TRUE,logfc.threshold = 0.5,verbose = TRUE)

#filter cluster 1
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[!(mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('1'))]]
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

# redo clustering and umap ---------------------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 12
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1),group.by = 'RNA_snn_res.0.7',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

FeaturePlot(object = mouse_multiome_Seurat,features = c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Meis2','Sp8','Cxcl14','Htr3a','Nr2f2','Calb2','Neurod6'),pt.size = 0.1,reduction = 'umap',order = TRUE)

#seems ok
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/mouse_multiome_meta_data_220910_0038_temp_data.rds')

# mouse multiome ATAC filter ----------------------------------------------
#load data
E145_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix/')
E145_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E145_2/outs/filtered_feature_bc_matrix/')
E155_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix/')
E155_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix/')
E155_3 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E155_3/outs/filtered_feature_bc_matrix/')
E165_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E165_2/outs/filtered_feature_bc_matrix/')

E145_1 <- E145_1$`Gene Expression`
E145_2 <- E145_2$`Gene Expression`
E155_1 <- E155_1$`Gene Expression`
E155_2 <- E155_2$`Gene Expression`
E155_3 <- E155_3$`Gene Expression`
E165_2 <- E165_2$`Gene Expression`

colnames(E145_1) <- paste('E145_1#',colnames(E145_1),sep = '')
colnames(E145_2) <- paste('E145_2#',colnames(E145_2),sep = '')
colnames(E155_1) <- paste('E155_1#',colnames(E155_1),sep = '')
colnames(E155_2) <- paste('E155_2#',colnames(E155_2),sep = '')
colnames(E155_3) <- paste('E155_3#',colnames(E155_3),sep = '')
colnames(E165_2) <- paste('E165_2#',colnames(E165_2),sep = '')

E145_1 <- CreateSeuratObject(counts = E145_1,project = 'E145_1',assay = 'RNA',min.cells = 0,min.features = 0)
E145_2 <- CreateSeuratObject(counts = E145_2,project = 'E145_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_1 <- CreateSeuratObject(counts = E155_1,project = 'E155_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_2 <- CreateSeuratObject(counts = E155_2,project = 'E155_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_3 <- CreateSeuratObject(counts = E155_3,project = 'E155_3',assay = 'RNA',min.cells = 0,min.features = 0)
E165_2 <- CreateSeuratObject(counts = E165_2,project = 'E165_2',assay = 'RNA',min.cells = 0,min.features = 0)

MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
E145_1[['percent.mt']] <- PercentageFeatureSet(object = E145_1,features = MT_gene_list)
E145_2[['percent.mt']] <- PercentageFeatureSet(object = E145_2,features = MT_gene_list)
E155_1[['percent.mt']] <- PercentageFeatureSet(object = E155_1,features = MT_gene_list)
E155_2[['percent.mt']] <- PercentageFeatureSet(object = E155_2,features = MT_gene_list)
E155_3[['percent.mt']] <- PercentageFeatureSet(object = E155_3,features = MT_gene_list)
E165_2[['percent.mt']] <- PercentageFeatureSet(object = E165_2,features = MT_gene_list)

mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

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

mouse_multiome_Seurat$batch <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_1','E155_1','E155_2'),"batch"] <- '220707'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_2','E155_3','E165_2'),"batch"] <- '220907'

#load meta data
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)[meta_data$filted == 'preserved']]

#save memory
remove(list = c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2'))
gc()

#pre-process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 20
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1),group.by = 'RNA_snn_res.0.5',label = TRUE)

#display
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'Age',label = FALSE)
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'batch',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#seems ok

#doublet filter
E145_1_doublet <- read.csv(file = './res/step_42_fig_220707/E145_1_scrublet.csv')
E145_2_doublet <- read.csv(file = './res/step_55_fig_220907/E145_2_scrublet.csv')
E155_1_doublet <- read.csv(file = './res/step_42_fig_220707/E155_1_scrublet.csv')
E155_2_doublet <- read.csv(file = './res/step_42_fig_220707/E155_2_scrublet.csv')
E155_3_doublet <- read.csv(file = './res/step_55_fig_220907/E155_3_scrublet.csv')
E165_2_doublet <- read.csv(file = './res/step_55_fig_220907/E165_2_scrublet.csv')

E145_1_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E145_1_doublet$barcodes,fixed = TRUE)
E145_2_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E145_2_doublet$barcodes,fixed = TRUE)
E155_1_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_1_doublet$barcodes,fixed = TRUE)
E155_2_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_2_doublet$barcodes,fixed = TRUE)
E155_3_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E155_3_doublet$barcodes,fixed = TRUE)
E165_2_doublet$barcodes <- sub(pattern = '.',replacement = '-',x = E165_2_doublet$barcodes,fixed = TRUE)

E145_1_doublet$barcodes <- paste('E145_1#',E145_1_doublet$barcodes,sep = '')
E145_2_doublet$barcodes <- paste('E145_2#',E145_2_doublet$barcodes,sep = '')
E155_1_doublet$barcodes <- paste('E155_1#',E155_1_doublet$barcodes,sep = '')
E155_2_doublet$barcodes <- paste('E155_2#',E155_2_doublet$barcodes,sep = '')
E155_3_doublet$barcodes <- paste('E155_3#',E155_3_doublet$barcodes,sep = '')
E165_2_doublet$barcodes <- paste('E165_2#',E165_2_doublet$barcodes,sep = '')

rownames(E145_1_doublet) <- E145_1_doublet$barcodes
rownames(E145_2_doublet) <- E145_2_doublet$barcodes
rownames(E155_1_doublet) <- E155_1_doublet$barcodes
rownames(E155_2_doublet) <- E155_2_doublet$barcodes
rownames(E155_3_doublet) <- E155_3_doublet$barcodes
rownames(E165_2_doublet) <- E165_2_doublet$barcodes
temp <- rbind(E145_1_doublet,E145_2_doublet,E155_1_doublet,E155_2_doublet,E155_3_doublet,E165_2_doublet)

mouse_multiome_Seurat$doublet_score <- NA
mouse_multiome_Seurat$doublet <- NA

mouse_multiome_Seurat$doublet_score <- temp[colnames(mouse_multiome_Seurat),"score"]
mouse_multiome_Seurat$doublet <- temp[colnames(mouse_multiome_Seurat),"prediction"]

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'donor',label = FALSE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('doublet_score'),reduction = 'umap') + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#hard to tell doublet

#nFeature nCount filter
VlnPlot(object = mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA'),pt.size = 0,group.by = 'RNA_snn_res.1.1')

#cluster 21 may be doublet
temp <- FindMarkers(object = mouse_multiome_Seurat,ident.1 = '21',group.by = 'RNA_snn_res.1.1',assay = 'RNA',slot = 'data',only.pos = TRUE)

#do not know what it is, filter it
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[!(mouse_multiome_Seurat$RNA_snn_res.1.1 %in% c('21'))]]

p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'donor',label = FALSE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('doublet_score'),reduction = 'umap') + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)


# redo clustering and annotation ------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 25
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1),group.by = 'RNA_snn_res.0.5',label = TRUE)

FeaturePlot(object = mouse_multiome_Seurat,features = c('nFeature_RNA','nCount_RNA','doublet_score'))

#seems all right

# redo human region analysis and macaque cell type prediction -------------
#convert gene id
mouse_to_human <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
temp <- mouse_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)
mouse_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

#human region analysis
human_Whole_Brain <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_whole_brain_RNA_Seurat_nFeature_3000_dim_50_220709.rds')
gene_list <- dplyr::intersect(x = VariableFeatures(human_Whole_Brain),y = rownames(mouse_multiome_Seurat@assays$converted@counts))
human_Whole_Brain <- my_process_seurat(object = human_Whole_Brain,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA'),npcs = 50,preprocess = TRUE)
human_Whole_Brain <- my_process_seurat(object = human_Whole_Brain,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 50,resolution = 1,group.by = 'cell.type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$converted@scale.data,object = human_Whole_Brain,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

anchors <- my_FindTransferAnchors(reference = human_Whole_Brain,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'converted',l2.norm = TRUE,dims = 1:50,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = human_Whole_Brain$area,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$human_area <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$human_area_score <- mouse_multiome_Seurat$prediction.score.max
predictions <- TransferData(anchorset = anchors,refdata = human_Whole_Brain$cell.type,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$human_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$human_cell_type_score <- mouse_multiome_Seurat$prediction.score.max

#save meta data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/human_whole_brain_predicted_mouse_multiome_Seurat_meta_data_220911_1017_temp_data.rds')

#display
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE)

# macaque cell type prediction --------------------------------------------
# load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')
gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(mouse_multiome_Seurat@assays$converted@counts))
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 32,resolution = 1,group.by = 'cell_type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#project
temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$converted@scale.data,object = macaque_multiome_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'converted',l2.norm = TRUE,dims = 1:32,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$cell_type,l2.norm = TRUE,dims = 1:32,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$macaque_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$macaque_cell_type_score <- mouse_multiome_Seurat$prediction.score.max

#display
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)

#save meta data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/macaque_cell_type_predicted_mouse_multiome_Seurat_meta_data_220911_1141_temp_data.rds')

# finally check -----------------------------------------------------------
meta_data <- readRDS(file = './res/step_55_fig_220907/mouse_multiome_meta_data_220909_1638_temp_data.rds')
mouse_multiome_Seurat$doublet_score <- meta_data[colnames(mouse_multiome_Seurat),"doublet_score"]
VlnPlot(object = mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','doublet_score'),group.by = 'RNA_snn_res.0.5',pt.size = 0)

#cluster 6 9 14 15 16 may be doublet
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#OPC End Per Mic VLMC marker
DefaultAssay(mouse_multiome_Seurat) <- 'RNA'
FeaturePlot(object = mouse_multiome_Seurat,features = c('Col1a1','Lum','Cx3cr1','Cldn5','Pecam1','Pdgfrb','Sox10','Olig2','Egfr'),reduction = 'umap')
VlnPlot(object = mouse_multiome_Seurat,features = c('Col1a1','Lum','Cx3cr1','Cldn5','Pecam1','Pdgfrb','Sox10','Olig2','Egfr'),group.by = 'RNA_snn_res.0.5',pt.size = 0,assay = 'RNA')

#no need to filter

#brain region and cell type filter
p1 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p4 <- DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'Age',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+p4+plot_layout(ncol = 2)

FeaturePlot(object = mouse_multiome_Seurat,features = c('Sox2','Aqp4','Hes5','Hes1','Pax6','Eomes'),reduction = 'umap',order = TRUE)

#save meta_data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_55_fig_220907/mouse_multiome_meta_data_220911_1211_temp_data.rds')

#save data
saveRDS(object = mouse_multiome_Seurat,file = '/data/User/sunym/temp/mouse_multiome_Seurat_220911.rds')