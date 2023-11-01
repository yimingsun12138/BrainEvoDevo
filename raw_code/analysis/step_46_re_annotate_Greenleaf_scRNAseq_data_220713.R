#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: re_annotate Greenleaf scRNAseq data                             ##
## Data: 2022.07.13                                                                ##
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

# integrate macaque multiome and RNA data ---------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_multiome_Seurat_220712.rds')
macaque_RNA_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_RNA_Seurat_220712.rds')

#merge data
dim(macaque_multiome_Seurat@assays$RNA@counts)
dim(macaque_RNA_Seurat@assays$RNA@counts)
gene_list <- dplyr::intersect(x = rownames(macaque_multiome_Seurat@assays$RNA@counts),y = rownames(macaque_RNA_Seurat@assays$RNA@counts))

macaque_integration_Seurat <- cbind(macaque_multiome_Seurat@assays$RNA@counts[gene_list,],macaque_RNA_Seurat@assays$RNA@counts[gene_list,])
macaque_integration_Seurat <- CreateSeuratObject(counts = macaque_integration_Seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)

#add meta data
macaque_integration_Seurat$tech <- NA
macaque_integration_Seurat$batch <- NA
macaque_integration_Seurat$donor <- NA
macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat$sub_cell_type <- NA

macaque_integration_Seurat@meta.data[colnames(macaque_multiome_Seurat),"tech"] <- macaque_multiome_Seurat$tech
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_Seurat),"tech"] <- macaque_RNA_Seurat$tech
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_Seurat),"batch"] <- macaque_multiome_Seurat$batch
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_Seurat),"batch"] <- macaque_RNA_Seurat$batch
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_Seurat),"donor"] <- macaque_multiome_Seurat$donor
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_Seurat),"donor"] <- macaque_RNA_Seurat$donor
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_Seurat),"cell_type"] <- macaque_multiome_Seurat$cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_Seurat),"cell_type"] <- macaque_RNA_Seurat$cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_Seurat),"sub_cell_type"] <- macaque_multiome_Seurat$sub_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_Seurat),"sub_cell_type"] <- macaque_RNA_Seurat$sub_cell_type

#process macaque integration Seurat
macaque_integration_Seurat <- my_process_seurat(object = macaque_integration_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','tech','batch','donor'),npcs = 50,preprocess = TRUE)
macaque_integration_Seurat <- my_process_seurat(object = macaque_integration_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 28,resolution = 1,group.by = 'cell_type',label = TRUE)
DimPlot(macaque_integration_Seurat,group.by = 'tech',label = TRUE,repel = TRUE)
DimPlot(macaque_integration_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1 <- DimPlot(macaque_integration_Seurat,group.by = 'RNA_snn_res.1',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear))) + theme(aspect.ratio = 1)
p2 <- DimPlot(macaque_integration_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear))) + theme(aspect.ratio = 1)

pdf(file = './res/step_46_fig_220713/macaque_integration_Seurat_dimplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#save data
saveRDS(object = macaque_integration_Seurat,file = './processed_data/220711_summary/macaque_integration_Seurat_220713.rds')

# harmony macaque multiome data and Greenleaf data round 1 ------------------------
#load data
macaque_integration_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_integration_Seurat_220713.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$Age != 'pcw24']

#convert gene name
temp <- macaque_integration_Seurat@assays$RNA@counts
macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
macaque_to_human_anno <- macaque_to_human_anno[,c(2,1,4,3)]
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 4)

macaque_integration_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_integration_Seurat <- my_process_seurat(object = macaque_integration_Seurat,assay = 'converted',reduction.name = 'PCA',nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(x = VariableFeatures(macaque_integration_Seurat),y = rownames(Greenleaf_RNA_Seurat@assays$converted@counts))

#harmony integration
macaque_integration_Seurat$species <- 'macaque'
Greenleaf_RNA_Seurat$species <- 'human'
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque = macaque_integration_Seurat,human = Greenleaf_RNA_Seurat),assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque = c('nCount_converted','tech','batch','donor'),human = c('nCount_converted','Age','Batch')),npcs = 50,
                                           reference_loading = 'macaque',integration_var = 'species',harmony_input_dim = 30,max.iter.harmony = 50,UMAP_dim = 30,resolution = 1)

Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"cell_type"] <- Greenleaf_RNA_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(macaque_integration_Seurat),"cell_type"] <- macaque_integration_Seurat$cell_type

p1 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
p2 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 1)

#Greenleaf data
Greenleaf_RNA_Seurat$harmony_cluster <- Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"seurat_clusters"]
Greenleaf_RNA_Seurat@reductions$harmony_UMAP <- CreateDimReducObject(embeddings = Brain_RNA_Seurat@reductions$umap@cell.embeddings[colnames(Greenleaf_RNA_Seurat),],assay = 'converted',key = 'harmony_')

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('nCount_originalexp','nFeature_originalexp'),reduction = 'harmony_UMAP',max.cutoff = 'q99')
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'harmony_cluster',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear,ArchRPalettes$ironMan)))

pdf(file = './res/step_46_fig_220713/Greenleaf_RNA_Seurat_harmony_round_1_show_nFeature_nCount_dimplot.pdf',width = 16,height = 12)
p1+p2+p3+plot_layout(ncol = 2)
dev.off()

#filter all other cells in In
cell_list <- colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$harmony_cluster %in% c('5','3','31','16','7','13','18','33') & !(Greenleaf_RNA_Seurat$cell_type %in% c('InCGE','InMGE'))]
table(Greenleaf_RNA_Seurat@meta.data[cell_list,"cell_type"]) #mostly GluN-3
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,!(colnames(Greenleaf_RNA_Seurat) %in% cell_list)]

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('nCount_originalexp','nFeature_originalexp'),reduction = 'harmony_UMAP',max.cutoff = 'q99')
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'harmony_cluster',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear,ArchRPalettes$ironMan)))
p1+p2+p3+plot_layout(ncol = 2)


# harmony macaque multiome data and Greenleaf data round 2 ----------------
macaque_integration_Seurat$species <- 'macaque'
Greenleaf_RNA_Seurat$species <- 'human'
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque = macaque_integration_Seurat,human = Greenleaf_RNA_Seurat),assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque = c('nCount_converted','tech','batch','donor'),human = c('nCount_converted','Age','Batch')),npcs = 50,
                                           reference_loading = 'macaque',integration_var = 'species',harmony_input_dim = 30,max.iter.harmony = 50,UMAP_dim = 30,resolution = 1)

Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"cell_type"] <- Greenleaf_RNA_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(macaque_integration_Seurat),"cell_type"] <- macaque_integration_Seurat$cell_type

p1 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
p2 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 1)

#Greenleaf data
Greenleaf_RNA_Seurat$harmony_cluster <- Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"seurat_clusters"]
Greenleaf_RNA_Seurat@reductions$harmony_UMAP <- CreateDimReducObject(embeddings = Brain_RNA_Seurat@reductions$umap@cell.embeddings[colnames(Greenleaf_RNA_Seurat),],assay = 'converted',key = 'harmony_')

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('nCount_originalexp','nFeature_originalexp'),reduction = 'harmony_UMAP',max.cutoff = 'q99')
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'harmony_cluster',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear,ArchRPalettes$ironMan)))

pdf(file = './res/step_46_fig_220713/Greenleaf_RNA_Seurat_harmony_round_2_show_nFeature_nCount_dimplot.pdf',width = 16,height = 12)
p1+p2+p3+plot_layout(ncol = 2)
dev.off()

#filter unwanted cells
table(Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$harmony_cluster == '26']$cell_type)
cell_list <- colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$harmony_cluster %in% c('26') & !(Greenleaf_RNA_Seurat$cell_type %in% c('InCGE','InMGE'))]
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,!(colnames(Greenleaf_RNA_Seurat) %in% cell_list)]

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('nCount_originalexp','nFeature_originalexp'),reduction = 'harmony_UMAP',max.cutoff = 'q99')
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'harmony_cluster',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear,ArchRPalettes$ironMan)))
p1+p2+p3+plot_layout(ncol = 2)

# harmony macaque multiome data and Greenleaf data round 3 ----------------
macaque_integration_Seurat$species <- 'macaque'
Greenleaf_RNA_Seurat$species <- 'human'
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque = macaque_integration_Seurat,human = Greenleaf_RNA_Seurat),assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque = c('nCount_converted','tech','batch','donor'),human = c('nCount_converted','Age','Batch')),npcs = 50,
                                           reference_loading = 'macaque',integration_var = 'species',harmony_input_dim = 25,max.iter.harmony = 50,UMAP_dim = 25,resolution = 1)

Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"cell_type"] <- Greenleaf_RNA_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(macaque_integration_Seurat),"cell_type"] <- macaque_integration_Seurat$cell_type

p1 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
p2 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 1)

#label transfer
temp_macaque <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'macaque']
temp_human <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'human']
temp_macaque$sub_cell_type <- macaque_integration_Seurat@meta.data[colnames(temp_macaque),"sub_cell_type"]

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:25,features = gene_list,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:25,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)
Greenleaf_RNA_Seurat$macaque_cell_type <- temp_human@meta.data[colnames(Greenleaf_RNA_Seurat),"predicted.id"]
Greenleaf_RNA_Seurat$macaque_cell_type_score <- temp_human@meta.data[colnames(Greenleaf_RNA_Seurat),"prediction.score.max"]

predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$sub_cell_type,l2.norm = FALSE,dims = 1:25,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)
Greenleaf_RNA_Seurat$macaque_sub_cell_type <- temp_human@meta.data[colnames(Greenleaf_RNA_Seurat),"predicted.id"]
Greenleaf_RNA_Seurat$macaque_sub_cell_type_score <- temp_human@meta.data[colnames(Greenleaf_RNA_Seurat),"prediction.score.max"]

Greenleaf_RNA_Seurat@reductions$harmony_PCA <- CreateDimReducObject(embeddings = temp_human@reductions$pca@cell.embeddings[colnames(Greenleaf_RNA_Seurat),],assay = 'converted')
Greenleaf_RNA_Seurat@reductions$harmony_UMAP <- CreateDimReducObject(embeddings = temp_human@reductions$umap@cell.embeddings[colnames(Greenleaf_RNA_Seurat),],assay = 'converted')

#check marker
DimPlot(object = temp_human,group.by = 'predicted.id',label = TRUE,repel = TRUE) + 
  DimPlot(object = temp_human,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  FeaturePlot(object = temp_human,features = c('EOMES'),slot = 'data') + 
  FeaturePlot(object = temp_human,features = c('PPP1R17'),slot = 'data') + 
  plot_layout(ncol = 2)
#checked

scibet::Confusion_heatmap(ori = temp_human$predicted.id,prd = temp_human$cell_type)
DimPlot(Greenleaf_RNA_Seurat,reduction = 'UMAP',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)

#filter GluN miss assigned as End
table(Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$macaque_cell_type == 'End']$cell_type)
hist(Greenleaf_RNA_Seurat$macaque_cell_type_score)
hist(Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$macaque_cell_type == 'End' & Greenleaf_RNA_Seurat$cell_type != 'End']$macaque_cell_type_score)

cell_list <- colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$macaque_cell_type == 'End' & Greenleaf_RNA_Seurat$cell_type != 'End']
table(Greenleaf_RNA_Seurat@meta.data[cell_list,"cell_type"])
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,!(colnames(Greenleaf_RNA_Seurat) %in% cell_list)]

p1 <- DimPlot(object = Greenleaf_RNA_Seurat,reduction = 'harmony_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,reduction = 'harmony_UMAP',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = Greenleaf_RNA_Seurat,reduction = 'harmony_UMAP',group.by = 'macaque_sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)

DimPlot(object = Greenleaf_RNA_Seurat,reduction = 'UMAP',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)

#save meta data
meta_data <- Greenleaf_RNA_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_46_fig_220713/Greenleaf_integrated_Seurat_meta_data.rds')

pca_embedding <- Greenleaf_RNA_Seurat@reductions$harmony_PCA@cell.embeddings
umap_embedding <- Greenleaf_RNA_Seurat@reductions$harmony_UMAP@cell.embeddings

saveRDS(object = pca_embedding,file = './res/step_46_fig_220713/Greenleaf_integrated_Seurat_pca_embedding.rds')
saveRDS(object = umap_embedding,file = './res/step_46_fig_220713/Greenleaf_integrated_Seurat_umap_embedding.rds')

#save Brain_RNA_Seurat
cell_list <- c(colnames(macaque_integration_Seurat),colnames(Greenleaf_RNA_Seurat))
Brain_RNA_Seurat <- Brain_RNA_Seurat[,cell_list]

Brain_RNA_Seurat@meta.data[colnames(macaque_integration_Seurat),"cell_type"] <- macaque_integration_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"cell_type"] <- Greenleaf_RNA_Seurat$macaque_cell_type
Brain_RNA_Seurat$sub_cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(macaque_integration_Seurat),"sub_cell_type"] <- macaque_integration_Seurat$sub_cell_type
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"sub_cell_type"] <- Greenleaf_RNA_Seurat$macaque_sub_cell_type

DimPlot(Brain_RNA_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
DimPlot(Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
DimPlot(Brain_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)

DimPlot(Brain_RNA_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + 
  DimPlot(Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)

saveRDS(object = Brain_RNA_Seurat,file = './processed_data/220711_summary/macaque_integration_Seurat_Greenleaf_RNA_Seurat_harmony_integration_Seurat.rds')

# save Greenleaf RNA data -----------------------------------------------------
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
meta_data <- readRDS(file = './res/step_46_fig_220713/Greenleaf_integrated_Seurat_meta_data.rds')
pca_embedding <- readRDS(file = './res/step_46_fig_220713/Greenleaf_integrated_Seurat_pca_embedding.rds')
umap_embedding <- readRDS(file = './res/step_46_fig_220713/Greenleaf_integrated_Seurat_umap_embedding.rds')

Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,rownames(meta_data)]
Greenleaf_RNA_Seurat@meta.data[1:3,]
Greenleaf_RNA_Seurat$macaque_cell_type <- meta_data[colnames(Greenleaf_RNA_Seurat),"macaque_cell_type"]
Greenleaf_RNA_Seurat$macaque_cell_type_score <- meta_data[colnames(Greenleaf_RNA_Seurat),"macaque_cell_type_score"]
Greenleaf_RNA_Seurat$macaque_sub_cell_type <- meta_data[colnames(Greenleaf_RNA_Seurat),"macaque_sub_cell_type"]
Greenleaf_RNA_Seurat$macaque_sub_cell_type_score <- meta_data[colnames(Greenleaf_RNA_Seurat),"macaque_sub_cell_type_score"]

Greenleaf_RNA_Seurat@reductions$harmony_PCA <- CreateDimReducObject(embeddings = pca_embedding[colnames(Greenleaf_RNA_Seurat),],assay = 'converted')
Greenleaf_RNA_Seurat@reductions$harmony_UMAP <- CreateDimReducObject(embeddings = umap_embedding[colnames(Greenleaf_RNA_Seurat),],assay = 'converted',key = 'harmony_')

p1 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'macaque_sub_cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
pdf(file = './res/step_46_fig_220713/Greenleaf_RNA_Seurat_cell_type_dimplot_on_ori_UMAP.pdf',width = 18,height = 9)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'macaque_sub_cell_type',label = TRUE,repel = TRUE,reduction = 'harmony_UMAP') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
pdf(file = './res/step_46_fig_220713/Greenleaf_RNA_Seurat_cell_type_dimplot_on_harmony_UMAP.pdf',width = 18,height = 9)
p1+p2+plot_layout(ncol = 2)
dev.off()

#save data
saveRDS(object = Greenleaf_RNA_Seurat,file = './processed_data/220711_summary/Greenleaf_RNA_Seurat_220714.rds')

#compare original annotation
scibet::Confusion_heatmap(ori = Greenleaf_RNA_Seurat$cell_type,prd = Greenleaf_RNA_Seurat$macaque_cell_type)
temp <- readRDS(file = './processed_data/220305_summary/macaque_integration_greenleaf_RNA_harmony_label_transfer_seurat_220309.rds')
temp <- temp[,temp$dataset == 'human']

cell_list <- dplyr::intersect(x = colnames(temp),y = colnames(Greenleaf_RNA_Seurat))
scibet::Confusion_heatmap(ori = temp@meta.data[cell_list,"cell_type"],prd = Greenleaf_RNA_Seurat@meta.data[cell_list,"macaque_sub_cell_type"])
