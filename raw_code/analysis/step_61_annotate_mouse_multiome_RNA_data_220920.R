#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: annotate mouse multiome RNA data                                ##
## Data: 2022.09.20                                                                ##
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
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# load data ---------------------------------------------------------------

#express matrix
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

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

#filter cells
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')
meta_data <- meta_data[meta_data$filted == 'preserved',]
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]

#add meta data
mouse_multiome_Seurat$species <- 'mouse'
mouse_multiome_Seurat$donor <- meta_data$donor
mouse_multiome_Seurat$Age <- meta_data$Age
mouse_multiome_Seurat$batch <- meta_data$batch

# process data ------------------------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 20,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'RNA_snn_res.0.3',label = TRUE)

#filter liuyt found RG subtype
temp <- readRDS(file = '/data/User/sunym/trash/mouse_multiome_Seurat_220911_Converted.rds')
DimPlot(object = temp,group.by = 'RNA_snn_res.0.3',label = TRUE,repel = TRUE)

#filter cluster 6
temp <- temp[,temp$RNA_snn_res.0.3 != '6']
cell_list <- colnames(temp)

#saveRDS
saveRDS(object = cell_list,file = './res/step_61_fig_220921/liuyt_filted_mouse_multiome_Seurat_cell_list_220921.rds')

# redo process ------------------------------------------------------------
#express matrix
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

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

#filter cells
cell_list <- readRDS(file = './res/step_61_fig_220921/liuyt_filted_mouse_multiome_Seurat_cell_list_220921.rds')
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')
table(meta_data[cell_list,"filted"])
mouse_multiome_Seurat <- mouse_multiome_Seurat[,cell_list]

#add meta data
mouse_multiome_Seurat$species <- 'mouse'
mouse_multiome_Seurat$donor <- meta_data[colnames(mouse_multiome_Seurat),"donor"]
mouse_multiome_Seurat$Age <- meta_data[colnames(mouse_multiome_Seurat),"Age"]
mouse_multiome_Seurat$batch <- meta_data[colnames(mouse_multiome_Seurat),"batch"]

#process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 20,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'RNA_snn_res.0.3',label = TRUE)

# human region analysis and macaque cell type prediction ------------------
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
gc()
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

#display
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'human_area',label = TRUE,repel = TRUE)

#save meta data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_61_fig_220921/human_whole_brain_predicted_mouse_multiome_Seurat_meta_data_220921_1103_temp_data.rds')

#remove data
rm(list = c('human_Whole_Brain'))
gc()
rm(list = c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2'))
gc()

# macaque cell type prediction --------------------------------------------
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
DefaultAssay(object = mouse_multiome_Seurat) <- 'RNA'
FeaturePlot(object = mouse_multiome_Seurat,features = c('Gad2','Eomes'),slot = 'data')

#cycling IP really need further filter

#save meta data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_61_fig_220921/macaque_cell_type_predicted_mouse_multiome_Seurat_meta_data_220921_1122_temp_data.rds')

# filter cycling In lineage IP --------------------------------------------
temp <- readRDS(file = '/data/User/sunym/temp/mouse-multiome-interneuron-cyc-cellid.rds')
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE,cells.highlight = temp)

#filter
mouse_multiome_Seurat <- mouse_multiome_Seurat[,!(colnames(mouse_multiome_Seurat) %in% temp)]

#save cell list
cell_list <- colnames(mouse_multiome_Seurat)
saveRDS(object = cell_list,file = './res/step_61_fig_220921/liuyt_filted_mouse_multiome_Seurat_cell_list_220921_1126.rds')

# redo-process ------------------------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 20,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'RNA_snn_res.0.3',label = TRUE)
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)
FeaturePlot(object = mouse_multiome_Seurat,features = c('Sox10','Olig2'),slot = 'data')

#save data
saveRDS(object = mouse_multiome_Seurat,file = '/data/User/sunym/temp/mouse_multiome_Seurat_220921.rds')

# label transfer liuyt CCA integration ------------------------------------
#load data
Brain_RNA_Seurat <- readRDS(file = '/data/User/sunym/temp/Integration_threeSpecies-threeDatasets_cca_kanchor5ndim40_seurat_220921.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')

#separate
Brain_RNA_Seurat$species <- NA
Brain_RNA_Seurat@meta.data[Brain_RNA_Seurat$orig.ident %in% c('hft'),"species"] <- 'human'
Brain_RNA_Seurat@meta.data[Brain_RNA_Seurat$orig.ident %in% c('macaque'),"species"] <- 'macaque'
Brain_RNA_Seurat@meta.data[Brain_RNA_Seurat$orig.ident %in% c('E145','E155','E165'),"species"] <- 'mouse'
table(Brain_RNA_Seurat$species)

temp_macaque <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'macaque']
temp_human <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'human']
temp_mouse <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'mouse']

#label transfer human data
temp_macaque$cell_type <- macaque_multiome_Seurat@meta.data[colnames(temp_macaque),"cell_type"]
anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integrated',query_assay = 'integrated',l2.norm = FALSE,dims = 1:40,k.anchor = 5,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:40,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)
temp_human$macaque_cell_type <- temp_human$predicted.id
temp_human$macaque_cell_type_score <- temp_human$prediction.score.max
DimPlot(object = temp_human,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)

#label transfer mouse data
anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_mouse,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integrated',query_assay = 'integrated',l2.norm = FALSE,dims = 1:40,k.anchor = 5,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:40,verbose = TRUE)
temp_mouse <- AddMetaData(object = temp_mouse,metadata = predictions)
temp_mouse$macaque_cell_type <- temp_mouse$predicted.id
temp_mouse$macaque_cell_type_score <- temp_mouse$prediction.score.max
DimPlot(object = temp_mouse,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)

# harmony integration -----------------------------------------------------
#load mouse data
#express matrix
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

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)

#filter cells
cell_list <- readRDS(file = './res/step_61_fig_220921/liuyt_filted_mouse_multiome_Seurat_cell_list_220921_1126.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,cell_list]

#convert to human symbol
mouse_to_human <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
mouse_multiome_Seurat <- My_Convert_Homology_Gene_ID(express_matrix = mouse_multiome_Seurat,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)

#create mouse multiome human symbol
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)

#add meta data
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')
mouse_multiome_Seurat$species <- 'mouse'
mouse_multiome_Seurat$donor <- meta_data[colnames(mouse_multiome_Seurat),"donor"]
mouse_multiome_Seurat$Age <- meta_data[colnames(mouse_multiome_Seurat),"Age"]
mouse_multiome_Seurat$batch <- meta_data[colnames(mouse_multiome_Seurat),"batch"]

#process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 17,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'RNA_snn_res.0.3',label = TRUE)

#load data
human_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')

#harmony with pca
gene_list <- VariableFeatures(object = macaque_multiome_Seurat)
gene_list <- gene_list[gene_list %in% rownames(human_RNA_Seurat@assays$RNA@counts) & gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts)]
express_matrix <- cbind(human_RNA_Seurat@assays$RNA@data[gene_list,],macaque_multiome_Seurat@assays$RNA@data[gene_list,],mouse_multiome_Seurat@assays$RNA@counts[gene_list,])
temp <- c(rep('human',times = ncol(human_RNA_Seurat)),rep('macaque',times = ncol(macaque_multiome_Seurat)),rep('mouse',times = ncol(mouse_multiome_Seurat)))
names(temp) <- colnames(express_matrix)

embedding_matrix <- HarmonyMatrix(data_mat = express_matrix,meta_data = temp,do_pca = TRUE,npcs = 30,max.iter.harmony = 50,return_object = FALSE,verbose = TRUE,theta = 5)
rownames(embedding_matrix) <- names(temp)
colnames(embedding_matrix) <- paste('PC',1:30,sep = '_')

#create Brain_RNA_Seurat
express_matrix <- cbind(human_RNA_Seurat@assays$RNA@counts[gene_list,],macaque_multiome_Seurat@assays$RNA@counts[gene_list,],mouse_multiome_Seurat@assays$RNA@counts[gene_list,])
Brain_RNA_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'Brain',assay = 'RNA',min.cells = 0,min.features = 0)
Brain_RNA_Seurat$species <- temp
Brain_RNA_Seurat@reductions$pca <- CreateDimReducObject(embeddings = embedding_matrix,assay = 'RNA')
Brain_RNA_Seurat <- my_process_seurat(object = Brain_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'species',label = FALSE)

# Seurat CCA integration --------------------------------------------------
#load mouse data
#express matrix
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

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)

#filter cells
cell_list <- readRDS(file = './res/step_61_fig_220921/liuyt_filted_mouse_multiome_Seurat_cell_list_220921_1126.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,cell_list]

#convert to human symbol
mouse_to_human <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
mouse_multiome_Seurat <- My_Convert_Homology_Gene_ID(express_matrix = mouse_multiome_Seurat,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)

#create mouse multiome human symbol
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)

#add meta data
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')
mouse_multiome_Seurat$species <- 'mouse'
mouse_multiome_Seurat$donor <- meta_data[colnames(mouse_multiome_Seurat),"donor"]
mouse_multiome_Seurat$Age <- meta_data[colnames(mouse_multiome_Seurat),"Age"]
mouse_multiome_Seurat$batch <- meta_data[colnames(mouse_multiome_Seurat),"batch"]

#process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 17,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'RNA_snn_res.0.3',label = TRUE)

#load data
human_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')

#CCA integration
#select integration features
gene_list <- SelectIntegrationFeatures(object.list = list(human_RNA_Seurat,mouse_multiome_Seurat,macaque_multiome_Seurat),verbose = TRUE)

#find integration anchor
anchors <- FindIntegrationAnchors(object.list = list(human_RNA_Seurat,mouse_multiome_Seurat,macaque_multiome_Seurat),
                                  anchor.features = gene_list,reduction = 'cca',dims = 1:35,verbose = TRUE)

#integrate data
Brain_RNA_Seurat <- IntegrateData(anchorset = anchors,dims = 1:35,verbose = TRUE)

#process
Brain_RNA_Seurat$species <- NA
Brain_RNA_Seurat@meta.data[colnames(macaque_multiome_Seurat),"species"] <- 'macaque'
Brain_RNA_Seurat@meta.data[colnames(human_RNA_Seurat),"species"] <- 'human'
Brain_RNA_Seurat@meta.data[colnames(mouse_multiome_Seurat),"species"] <- 'mouse'

DefaultAssay(Brain_RNA_Seurat) <- 'integrated'
Brain_RNA_Seurat <- ScaleData(object = Brain_RNA_Seurat,features = rownames(Brain_RNA_Seurat@assays$integrated@data),assay = 'integrated',vars.to.regress = NULL,verbose = TRUE)
Brain_RNA_Seurat <- RunPCA(object = Brain_RNA_Seurat,assay = 'integrated',features = rownames(Brain_RNA_Seurat@assays$integrated@data),npcs = 50,verbose = TRUE)
Brain_RNA_Seurat <- FindNeighbors(object = Brain_RNA_Seurat,reduction = 'pca',dims = 1:35,assay = 'integrated',verbose = TRUE)
Brain_RNA_Seurat <- FindClusters(object = Brain_RNA_Seurat,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5),verbose = TRUE)
Brain_RNA_Seurat <- RunUMAP(object = Brain_RNA_Seurat,dims = 1:30,reduction = 'pca',assay = 'integrated',verbose = TRUE)

#plot
DimPlot(object = Brain_RNA_Seurat,group.by = 'species',label = FALSE)
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',split.by = 'species',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))

#run by for loop
for (i in c(28,29,30,31,32,36,38,40,45,50)) {
  
  ndim <- i
  
  #select integration features
  gene_list <- SelectIntegrationFeatures(object.list = list(human_RNA_Seurat,mouse_multiome_Seurat,macaque_multiome_Seurat),verbose = TRUE)
  
  #find integration anchor
  anchors <- FindIntegrationAnchors(object.list = list(human_RNA_Seurat,mouse_multiome_Seurat,macaque_multiome_Seurat),
                                    anchor.features = gene_list,reduction = 'cca',dims = 1:ndim,verbose = TRUE)
  
  #integrate data
  Brain_RNA_Seurat <- IntegrateData(anchorset = anchors,dims = 1:ndim,verbose = TRUE)
  
  #process
  Brain_RNA_Seurat$species <- NA
  Brain_RNA_Seurat@meta.data[colnames(macaque_multiome_Seurat),"species"] <- 'macaque'
  Brain_RNA_Seurat@meta.data[colnames(human_RNA_Seurat),"species"] <- 'human'
  Brain_RNA_Seurat@meta.data[colnames(mouse_multiome_Seurat),"species"] <- 'mouse'
  
  DefaultAssay(Brain_RNA_Seurat) <- 'integrated'
  Brain_RNA_Seurat <- ScaleData(object = Brain_RNA_Seurat,features = rownames(Brain_RNA_Seurat@assays$integrated@data),assay = 'integrated',vars.to.regress = NULL,verbose = TRUE)
  Brain_RNA_Seurat <- RunPCA(object = Brain_RNA_Seurat,assay = 'integrated',features = rownames(Brain_RNA_Seurat@assays$integrated@data),npcs = 50,verbose = TRUE)
  Brain_RNA_Seurat <- FindNeighbors(object = Brain_RNA_Seurat,reduction = 'pca',dims = 1:ndim,assay = 'integrated',verbose = TRUE)
  Brain_RNA_Seurat <- FindClusters(object = Brain_RNA_Seurat,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5),verbose = TRUE)
  Brain_RNA_Seurat <- RunUMAP(object = Brain_RNA_Seurat,dims = 1:ndim,reduction = 'pca',assay = 'integrated',verbose = TRUE)
  
  #save data
  char <- ''
  saveRDS(object = Brain_RNA_Seurat,file = char)
  gc()
}

# integrate using my harmony ----------------------------------------------
#load mouse data
#express matrix
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

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)

#filter cells
cell_list <- readRDS(file = './res/step_61_fig_220921/liuyt_filted_mouse_multiome_Seurat_cell_list_220921_1126.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,cell_list]

#convert to human symbol
mouse_to_human <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
mouse_multiome_Seurat <- My_Convert_Homology_Gene_ID(express_matrix = mouse_multiome_Seurat,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)

#create mouse multiome human symbol
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)

#add meta data
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')
mouse_multiome_Seurat$species <- 'mouse'
mouse_multiome_Seurat$donor <- meta_data[colnames(mouse_multiome_Seurat),"donor"]
mouse_multiome_Seurat$Age <- meta_data[colnames(mouse_multiome_Seurat),"Age"]
mouse_multiome_Seurat$batch <- meta_data[colnames(mouse_multiome_Seurat),"batch"]

#process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 17,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'RNA_snn_res.0.3',label = TRUE)

#load data
human_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')

#gene list
gene_list <- VariableFeatures(macaque_multiome_Seurat)
gene_list <- gene_list[gene_list %in% rownames(human_RNA_Seurat@assays$RNA@counts) & gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts)]
macaque_multiome_Seurat$species <- 'macaque'
human_RNA_Seurat$species <- 'human'
mouse_multiome_Seurat$species <- 'mouse'

#harmony integration
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque = macaque_multiome_Seurat,human = human_RNA_Seurat,mouse = mouse_multiome_Seurat),
                                           assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(macaque = c('nCount_RNA','donor'),human = NULL,mouse = NULL),
                                           npcs = 50,reference_loading = 'macaque',integration_var = 'species',harmony_input_dim = 35,max.iter.harmony = 50,UMAP_dim = 35,
                                           resolution = 1,lambda=1,theta=2,sigma=0.2)
my_send_sms('harmony done!')
#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#add meta_data
Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(macaque_multiome_Seurat),"cell_type"] <- macaque_multiome_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(human_RNA_Seurat),"cell_type"] <- as.character(human_RNA_Seurat$ReAnno_celltype)
Brain_RNA_Seurat@meta.data[colnames(mouse_multiome_Seurat),"cell_type"] <- 'mouse'

#plot
DimPlot(object = Brain_RNA_Seurat,group.by = 'dataset',label = FALSE)
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')

# label transfer liuyt CCA integration ------------------------------------
#load data
Brain_RNA_Seurat <- readRDS(file = '/data/User/sunym/temp/Integration_threeSpecies-threeDatasets_cca_kanchor5ndim30_HRMorder_seurat_220921.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')

#separate
Brain_RNA_Seurat$species <- NA
Brain_RNA_Seurat@meta.data[Brain_RNA_Seurat$orig.ident %in% c('hft'),"species"] <- 'human'
Brain_RNA_Seurat@meta.data[Brain_RNA_Seurat$orig.ident %in% c('macaque'),"species"] <- 'macaque'
Brain_RNA_Seurat@meta.data[Brain_RNA_Seurat$orig.ident %in% c('E145','E155','E165'),"species"] <- 'mouse'
table(Brain_RNA_Seurat$species)

temp_macaque <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'macaque']
temp_human <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'human']
temp_mouse <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'mouse']

#label transfer human data
temp_macaque$cell_type <- macaque_multiome_Seurat@meta.data[colnames(temp_macaque),"cell_type"]
anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integrated',query_assay = 'integrated',l2.norm = FALSE,dims = 1:30,k.anchor = 5,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:30,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)
temp_human$macaque_cell_type <- temp_human$predicted.id
temp_human$macaque_cell_type_score <- temp_human$prediction.score.max
DimPlot(object = temp_human,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)

#label transfer mouse data
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype
anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_mouse,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integrated',query_assay = 'integrated',l2.norm = FALSE,dims = 1:30,k.anchor = 5,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:30,verbose = TRUE)
temp_mouse <- AddMetaData(object = temp_mouse,metadata = predictions)
temp_mouse$macaque_cell_type <- temp_mouse$predicted.id
temp_mouse$macaque_cell_type_score <- temp_mouse$prediction.score.max
p1 <- DimPlot(object = temp_mouse,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE,cols = temp_col)
p2 <- DimPlot(object = temp_mouse,group.by = 'integrated_snn_res.1',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)
DefaultAssay(object = temp_mouse) <- 'RNA'
VlnPlot(object = temp_mouse,features = c('MKI67','TOP2A'),pt.size = 0,assay = 'RNA',group.by = 'integrated_snn_res.1')

#reanno
temp_mouse@meta.data[temp_mouse$integrated_snn_res.1 %in% c('18','23','11'),"macaque_cell_type"] <- 'RG-1'
temp_mouse@meta.data[temp_mouse$integrated_snn_res.1 %in% c('19','21','22'),"macaque_cell_type"] <- 'Cycling'
p1 <- DimPlot(object = temp_mouse,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE,cols = temp_col)
p2 <- DimPlot(object = temp_mouse,group.by = 'integrated_snn_res.1',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#create meta_data
meta_data <- data.frame(cell_id = c(colnames(temp_human),colnames(temp_macaque),colnames(temp_mouse)),
                        species = c(rep('human',times = ncol(temp_human)),rep('macaque',times = ncol(temp_macaque)),rep('mouse',times = ncol(temp_mouse))),
                        macaque_cell_type = c(temp_human$macaque_cell_type,temp_macaque$cell_type,temp_mouse$macaque_cell_type),
                        macaque_cell_type_score = c(temp_human$macaque_cell_type_score,rep(1,times = ncol(temp_macaque)),temp_mouse$macaque_cell_type_score))

#save data
saveRDS(object = meta_data,file = './res/step_61_fig_220921/liuyt_Brain_RNA_Seurat_integration_k_5_ndim_30_220922.rds')

# integrate with nature mouse data ----------------------------------------
#express matrix
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

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)

#filter cells
cell_list <- readRDS(file = './res/step_61_fig_220921/liuyt_filted_mouse_multiome_Seurat_cell_list_220921_1126.rds')
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,cell_list]

#create seurat object
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

#add meta data
mouse_multiome_Seurat$species <- 'mouse'
mouse_multiome_Seurat$donor <- meta_data[colnames(mouse_multiome_Seurat),"donor"]
mouse_multiome_Seurat$Age <- meta_data[colnames(mouse_multiome_Seurat),"Age"]
mouse_multiome_Seurat$batch <- meta_data[colnames(mouse_multiome_Seurat),"batch"]

meta_data <- readRDS(file = './res/step_61_fig_220921/liuyt_Brain_RNA_Seurat_integration_k_5_ndim_30_220922.rds')
mouse_multiome_Seurat$macaque_cell_type <- meta_data[colnames(mouse_multiome_Seurat),"macaque_cell_type"]

#process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 20,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'macaque_cell_type',label = TRUE)

#load public data
mouse_public_data <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/processed_data/mouse_E14_E15_E16_RNA_Seurat_220922.rds')


# harmony integration -----------------------------------------------------
gene_list <- VariableFeatures(mouse_public_data)
gene_list <- gene_list[gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts)]
mouse_multiome_Seurat$dataset <- 'multiome'
mouse_public_data$dataset <- 'public'

Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(multiome = mouse_multiome_Seurat,public = mouse_public_data),assay = 'RNA',
                                           variable_feature = gene_list,var_to_regress_list = list(multiome = NULL,public = c('nCount_RNA','biosample_id')),
                                           npcs = 50,reference_loading = 'public',integration_var = 'dataset',harmony_input_dim = 25,max.iter.harmony = 50,
                                           UMAP_dim = 25,resolution = 1)

Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(mouse_multiome_Seurat),"cell_type"] <- mouse_multiome_Seurat$macaque_cell_type
Brain_RNA_Seurat@meta.data[colnames(mouse_public_data),"cell_type"] <- mouse_public_data$New_cellType

#display
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
FeaturePlot(object = Brain_RNA_Seurat,features = 'Adarb2',split.by = 'dataset',order = TRUE)

# label transfer from public data -----------------------------------------
temp_multiome <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'multiome']
temp_public <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'public']
anchors <- my_FindTransferAnchors(reference = temp_public,query = temp_multiome,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integration',query_assay = 'integration',
                                  l2.norm = FALSE,dims = 1:25,features = gene_list,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = mouse_public_data$New_cellType,l2.norm = FALSE,dims = 1:25,verbose = TRUE)
temp_multiome <- AddMetaData(object = temp_multiome,metadata = predictions)
mouse_multiome_Seurat$mouse_cell_type <- temp_multiome@meta.data[colnames(mouse_multiome_Seurat),"predicted.id"]
mouse_multiome_Seurat$mouse_cell_type_score <- temp_multiome@meta.data[colnames(mouse_multiome_Seurat),"prediction.score.max"]
p1 <- DimPlot(temp_multiome,group.by = 'predicted.id',label = TRUE,repel = TRUE)
p2 <- DimPlot(temp_public,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#display
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'mouse_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

scibet::Confusion_heatmap(ori = mouse_multiome_Seurat$RNA_snn_res.0.7,prd = mouse_multiome_Seurat$mouse_cell_type)

#CPN marker
FeaturePlot(object = mouse_multiome_Seurat,features = 'Satb2')
FeaturePlot(object = mouse_multiome_Seurat,features = 'Eomes')

#annotate mouse multiome data
mouse_multiome_Seurat$cell_type <- NA

#End
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('14'),"cell_type"] <- 'End'

#Per
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('15'),"cell_type"] <- 'Per'

#VLMC
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('17'),"cell_type"] <- 'VLMC'

#Mic
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('16'),"cell_type"] <- 'Mic'

#In
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('11','5'),"cell_type"] <- 'In'

#Cyc
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('6'),"cell_type"] <- 'Cyc'

#Cyc-IP
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('13'),"cell_type"] <- 'Cyc-IP'

#RG
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('8'),"cell_type"] <- 'RG'

#IP
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('2'),"cell_type"] <- 'IP'

#Migrating Neuron
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('0','1','10','3'),"cell_type"] <- 'MGN'

#CPN
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('9'),"cell_type"] <- 'CPN'

#CThPN
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('7','12'),"cell_type"] <- 'CThPN'

#SCPN
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('4'),"cell_type"] <- 'SCPN'

#display
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'mouse_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#save meta_data
meta_data <- mouse_multiome_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_61_fig_220921/mouse_multiome_Seurat_meta_data_220922_2109.rds')

# finally process mouse multiome Seurat -----------------------------------
#express matrix
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

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)

#filter cells
meta_data <- readRDS(file = './res/step_61_fig_220921/mouse_multiome_Seurat_meta_data_220922_2109.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]

#create seurat object
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

#add meta data
mouse_multiome_Seurat$species <- 'mouse'
mouse_multiome_Seurat$donor <- meta_data[colnames(mouse_multiome_Seurat),"donor"]
mouse_multiome_Seurat$Age <- meta_data[colnames(mouse_multiome_Seurat),"Age"]
mouse_multiome_Seurat$batch <- meta_data[colnames(mouse_multiome_Seurat),"batch"]
mouse_multiome_Seurat$macaque_cell_type <- meta_data[colnames(mouse_multiome_Seurat),"macaque_cell_type"]
mouse_multiome_Seurat$mouse_cell_type <- meta_data[colnames(mouse_multiome_Seurat),"mouse_cell_type"]

#process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 20,resolution = seq(from = 0.3,to = 1.5,by = 0.1),group.by = 'macaque_cell_type',label = TRUE)

#display
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'mouse_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

scibet::Confusion_heatmap(ori = mouse_multiome_Seurat$RNA_snn_res.0.7,prd = mouse_multiome_Seurat$mouse_cell_type)

#annotate mouse multiome data
mouse_multiome_Seurat$cell_type <- NA

#End
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('14'),"cell_type"] <- 'End'

#Per
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('15'),"cell_type"] <- 'Per'

#VLMC
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('17'),"cell_type"] <- 'VLMC'

#Mic
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('16'),"cell_type"] <- 'Mic'

#In
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('11','5'),"cell_type"] <- 'In'

#Cyc
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('6'),"cell_type"] <- 'Cyc'

#Cyc-IP
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('13'),"cell_type"] <- 'Cyc-IP'

#RG
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('8'),"cell_type"] <- 'RG'

#IP
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('2'),"cell_type"] <- 'IP'

#Migrating Neuron
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('0','1','10','3'),"cell_type"] <- 'MGN'

#CPN
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('9'),"cell_type"] <- 'CPN'

#CThPN
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('7','12'),"cell_type"] <- 'CThPN'

#SCPN
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.7 %in% c('4'),"cell_type"] <- 'SCPN'

#display
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'mouse_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

#save data
saveRDS(object = mouse_multiome_Seurat,file = './processed_data/220802_summary/mouse_multiome_Seurat_220922.rds')
