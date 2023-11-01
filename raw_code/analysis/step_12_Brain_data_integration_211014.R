#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Brain data integration                                          ##
## Data: 2021.10.14                                                                ##
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
library(SummarizedExperiment)
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# data preparation --------------------------------------------------------
# macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
# DefaultAssay(macaque_old_seurat) <- 'RNA'
# 
# macaque_new_seurat <- readRDS(file = './res/step_11_fig_210926/macaque_RNA_210929_express_matrix_210929.rds')
# macaque_new_seurat <- CreateSeuratObject(counts = macaque_new_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
# temp <- base::lapply(colnames(macaque_new_seurat),function(c){
#   return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
# })
# temp <- unlist(temp)
# macaque_new_seurat$donor <- temp
# 
# table(duplicated(colnames(macaque_new_seurat),colnames(macaque_old_seurat)))
# table(rownames(macaque_new_seurat) == rownames(macaque_old_seurat))
# macaque_RNA_seurat <- cbind(macaque_new_seurat@assays$RNA@counts,macaque_old_seurat@assays$RNA@counts)
# macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
# temp <- base::lapply(colnames(macaque_RNA_seurat),function(c){
#   return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
# })
# temp <- unlist(temp)
# macaque_RNA_seurat$donor <- temp
# macaque_RNA_seurat@meta.data[,'batch'] <- NA
# macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"batch"] <- '200919'
# macaque_RNA_seurat@meta.data[colnames(macaque_new_seurat),"batch"] <- '210922'
# table(macaque_RNA_seurat$batch)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 40,resolution = c(0.5,1,1.5),group.by = 'batch',label = FALSE)
# my_send_sms(body = 'UMAP done!')
# 
# DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
# macaque_new_seurat <- readRDS(file = './processed_data/macaque_RNA_210929_seurat_annotated_211013.rds')
# macaque_RNA_seurat$cell_type <- NA
# macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"cell_type"] <- as.character(macaque_old_seurat$sub_cell_type)
# macaque_RNA_seurat@meta.data[colnames(macaque_new_seurat),"cell_type"] <- as.character(macaque_new_seurat$merge_anno)
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('20','24'),'cell_type'] <- 'Unknown_Ex-3'
# saveRDS(macaque_RNA_seurat,file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211014.rds')

# integrate with atac gene activity ---------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211014.rds')
macaque_ATAC <- readRDS(file = './res/step_12_fig_211014/scATAC-GeneScoreMatrix.rds')

temp <- assay(macaque_ATAC)
gene_list <- rowData(macaque_ATAC)
gene_list$id <- paste('a',1:17321,sep = '-')
rownames(temp) <- as.character(gene_list$id)
gene_list <- gene_list[,c('name','id')]
gene_list <- as.data.frame(gene_list)
gene_list <- cbind(gene_list,gene_list)
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = gene_list,filter_anno = FALSE)

macaque_ATAC_seurat <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- colData(macaque_ATAC)
temp <- as.data.frame(temp)
gc()
temp <- temp[colnames(macaque_ATAC_seurat),]

macaque_ATAC_seurat$sample <- as.character(temp[colnames(macaque_ATAC_seurat),"Sample"])
macaque_ATAC_seurat$cell_type <- as.character(temp[colnames(macaque_ATAC_seurat),"Celltype2"])
macaque_ATAC_seurat$Cluster <- as.character(temp[colnames(macaque_ATAC_seurat),"Cluster_HarmonyPeak_res.1"])

temp <- macaque_RNA_seurat
temp <- FindVariableFeatures(temp,selection.method = 'vst',nfeatures = 4000)
gene_list <- VariableFeatures(temp)
gene_list <- gene_list[gene_list %in% rownames(macaque_ATAC_seurat)]
rm(temp)
gc()
#harmony integration
macaque_RNA_seurat$dataset <- 'RNA'
macaque_ATAC_seurat$dataset <- 'ATAC'
test_seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,ATAC=macaque_ATAC_seurat),assay = 'RNA',variable_feature = gene_list,
                                      var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),ATAC=NULL),npcs = 50,reference_loading = 'RNA',
                                      integration_var = 'dataset',harmony_input_dim = 35,max.iter.harmony = 50,UMAP_dim = 35,resolution = 1,kmeans_init_iter_max = 200,
                                      yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/')

test_seurat$cell_type <- NA
test_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)
test_seurat@meta.data[colnames(macaque_ATAC_seurat),"cell_type"] <- as.character(macaque_ATAC_seurat$cell_type)
my_send_sms(body = 'harmony done!')


pdf(file = './res/step_12_fig_211014/macaque_RNA_ATAC_harmony_umap.pdf',width = 18,height = 8)
DimPlot(object = test_seurat,group.by = 'cell_type',split.by = 'dataset',label = TRUE,repel = TRUE)
dev.off()

test_seurat$Cluster <- NA
test_seurat@meta.data[colnames(macaque_ATAC_seurat),"Cluster"] <- as.character(macaque_ATAC_seurat$Cluster)
DimPlot(object = test_seurat[,test_seurat$dataset == 'ATAC'],group.by = 'Cluster',label = TRUE,repel = TRUE)

test_seurat$sample <- NA
test_seurat@meta.data[colnames(macaque_ATAC_seurat),"sample"] <- as.character(macaque_ATAC_seurat$sample)
test_seurat@meta.data[colnames(macaque_RNA_seurat),"sample"] <- as.character(macaque_RNA_seurat$donor)
DimPlot(object = test_seurat,group.by = 'sample',label = FALSE,split.by = 'dataset')

# integrate with greenleaf ------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211014.rds')
greenleaf_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
greenleaf_RNA_seurat <- greenleaf_RNA_seurat[,greenleaf_RNA_seurat$Age != 'pcw24']

macaque_RNA_seurat$species <- 'macaque'
greenleaf_RNA_seurat$species <- 'human'
#convert gene id
temp <- macaque_RNA_seurat@assays$RNA@counts
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')

temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp
rm(temp)
gc()

#find variable gene
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 5000,vars.to.regress = c('nCount_converted','donor','batch'),npcs = 50,preprocess = TRUE)
greenleaf_RNA_seurat <- my_process_seurat(object = greenleaf_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 5000,vars.to.regress = c('Age','Batch','nCount_converted'),npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),VariableFeatures(greenleaf_RNA_seurat))

#harmony
cl <- makeCluster(1)
clusterExport(cl,c('macaque_RNA_seurat','greenleaf_RNA_seurat','gene_list'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/','/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))
Brain_RNA_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_harmony_integration(named_seurat_list = list(macaque=macaque_RNA_seurat,human=greenleaf_RNA_seurat),assay = 'converted',variable_feature = gene_list,
                              var_to_regress_list = list(macaque=c('nCount_converted','donor','batch'),human=c('Age','Batch','nCount_converted')),npcs = 50,reference_loading = 'macaque',
                              integration_var = 'species',harmony_input_dim = 30,max.iter.harmony = 50,UMAP_dim = 30,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                              yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/')
  return(a)
})
stopCluster(cl)
Brain_RNA_seurat <- Brain_RNA_seurat[[1]]
gc()
my_send_sms('harmony done!')

Brain_RNA_seurat@meta.data[,'cell_type'] <- NA
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)
Brain_RNA_seurat@meta.data[colnames(greenleaf_RNA_seurat),"cell_type"] <- as.character(greenleaf_RNA_seurat$cell_type)
DimPlot(Brain_RNA_seurat,group.by = 'dataset',label = FALSE)

pdf(file = './res/step_12_fig_211014/macaque_RNA_greanleaf_RNA_harmony_umap.pdf',width = 16,height = 8)
DimPlot(Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')
dev.off()

#using marker gene to validate the unknown cells in greenleaf data
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']@reductions$umap
greenleaf_RNA_seurat@reductions$umap <- temp
rm(temp)
gc()

DimPlot(greenleaf_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'umap')
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque']@reductions$umap
macaque_RNA_seurat@reductions$UMAP <- temp
rm(temp)
gc()
#InPSB marker
DefaultAssay(greenleaf_RNA_seurat) <- 'converted'
greenleaf_RNA_seurat <- NormalizeData(object = greenleaf_RNA_seurat,normalization.method = "LogNormalize",scale.factor = 10000)
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = "LogNormalize",scale.factor = 10000)

p1 <- geneSetAveragePlot(genes = c('GAD1','GAD2','MEIS2'),object = greenleaf_RNA_seurat,object.class = 'seurat',assay = 'converted',
                         embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human InPSB marker')

p2 <- geneSetAveragePlot(genes = c('GAD1','GAD2','MEIS2'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque InPSB marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_greanleaf_RNA_harmony_InPSB_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#Unknown_InPSB marker
p1 <- geneSetAveragePlot(genes = c('RXRG','RARB','SV2C','RGS9'),object = greenleaf_RNA_seurat,object.class = 'seurat',assay = 'converted',
                         embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human Unknown InPSB marker')

p2 <- geneSetAveragePlot(genes = c('RXRG','RARB','SV2C','RGS9'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque Unknown InPSB marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_greanleaf_RNA_harmony_Unknown_InPSB_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#Unknow_Ex marker
p1 <- geneSetAveragePlot(genes = c('TLL1','NR4A2','B3GAT2','NTNG2'),object = greenleaf_RNA_seurat,object.class = 'seurat',assay = 'converted',
                         embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human Unknown Ex marker')

p2 <- geneSetAveragePlot(genes = c('TLL1','NR4A2','B3GAT2','NTNG2'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque Unknown Ex marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_greanleaf_RNA_harmony_Unknown_Ex_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#Unknown Ex-3 marker
p1 <- geneSetAveragePlot(genes = c('SLIT3','DGKG','ADAMTS3'),object = greenleaf_RNA_seurat,object.class = 'seurat',assay = 'converted',
                         embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human Unknown Ex-3 marker')

p2 <- geneSetAveragePlot(genes = c('SLIT3','DGKG','ADAMTS3'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque Unknown Ex-3 marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_greanleaf_RNA_harmony_Unknown_Ex-3_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#suggest: delete all the unknown cell type!

# integrate with neuron ---------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211014.rds')
human_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')

#convert data
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp

temp <- human_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno,filter_anno = FALSE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
human_RNA_seurat[['converted']] <- temp

rm(temp)
gc()

#pre-process
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 6000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
human_RNA_seurat <- my_process_seurat(object = human_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 6000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),VariableFeatures(human_RNA_seurat))

macaque_RNA_seurat$species <- 'macaque'
human_RNA_seurat$species <- 'human'
#harmony
cl <- makeCluster(1)
clusterExport(cl,c('macaque_RNA_seurat','human_RNA_seurat','gene_list'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/','/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))
Brain_RNA_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_harmony_integration(named_seurat_list = list(macaque=macaque_RNA_seurat,human=human_RNA_seurat),assay = 'converted',variable_feature = gene_list,
                              var_to_regress_list = list(macaque=c('nCount_converted','donor','batch'),human=c('Donor','Library','nCount_converted')),npcs = 50,reference_loading = 'macaque',
                              integration_var = 'species',harmony_input_dim = 35,max.iter.harmony = 50,UMAP_dim = 35,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                              yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/',sigma = 0.2)
  return(a)
})
stopCluster(cl)
Brain_RNA_seurat <- Brain_RNA_seurat[[1]]
gc()
my_send_sms('harmony done!')

DimPlot(object = Brain_RNA_seurat,group.by = 'dataset',label = FALSE)
Brain_RNA_seurat$cell_type <- NA
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)
Brain_RNA_seurat@meta.data[colnames(human_RNA_seurat),"cell_type"] <- as.character(human_RNA_seurat$Cluster)

pdf(file = './res/step_12_fig_211014/macaque_RNA_PDhuman_RNA_harmony_umap.pdf',width = 16,height = 8)
DimPlot(object = Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')
dev.off()

#modify original data
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque']@reductions$umap
macaque_RNA_seurat@reductions$UMAP <- temp
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']@reductions$umap
human_RNA_seurat@reductions$UMAP <- temp

DefaultAssay(macaque_RNA_seurat) <- 'RNA'
DefaultAssay(human_RNA_seurat) <- 'RNA'
macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)
human_RNA_seurat <- NormalizeData(object = human_RNA_seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)

#InPSB marker
p1 <- geneSetAveragePlot(genes = c('GAD1','GAD2','MEIS2'),object = human_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human InPSB marker')

p2 <- geneSetAveragePlot(genes = c('GAD1','GAD2','MEIS2'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque InPSB marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_PDhuman_RNA_InPSB_marker_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#Unknown_InPSB marker
p1 <- geneSetAveragePlot(genes = c('RXRG','RARB','SV2C','RGS9'),object = human_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human Unknown InPSB marker')

p2 <- geneSetAveragePlot(genes = c('RXRG','RARB','SV2C','RGS9'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque Unknown InPSB marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_PDhuman_RNA_Unknown_InPSB_marker_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#Unknown Ex marker
p1 <- geneSetAveragePlot(genes = c('TLL1','NR4A2','B3GAT2','NTNG2'),object = human_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human Unknown Ex marker')

p2 <- geneSetAveragePlot(genes = c('TLL1','NR4A2','B3GAT2','NTNG2'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque Unknown Ex marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_PDhuman_RNA_Unknown_Ex_marker_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#Unknown Ex-3 marker
p1 <- geneSetAveragePlot(genes = c('SLIT3','DGKG','ADAMTS3'),object = human_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'human Unknown Ex-3 marker')

p2 <- geneSetAveragePlot(genes = c('SLIT3','DGKG','ADAMTS3'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = FALSE,plot.title = 'macaque Unknown Ex-3 marker')

pdf(file = './res/step_12_fig_211014/macaque_RNA_PDhuman_RNA_Unknown_Ex_3_marker_featureplot.pdf',width = 14,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

#suggest: clean all unknown