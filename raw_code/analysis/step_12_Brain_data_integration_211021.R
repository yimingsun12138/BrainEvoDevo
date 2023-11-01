#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Brain data integration                                          ##
## Data: 2021.10.21                                                                ##
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

# macaque RNA seurat ------------------------------------------------------
# #load data
# macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211014.rds')
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# temp <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('Unknown_Ex','Unknown_Ex-3','Unknown_InPSB'))]
# 
# meta_data <- temp@meta.data
# meta_data <- meta_data[,c('donor','batch',"cell_type")]
# meta_data$cell_id <- rownames(meta_data)
# 
# temp <- temp@assays$RNA@counts
# meta_data <- meta_data[colnames(temp),]
# macaque_RNA_seurat <- CreateSeuratObject(counts = temp,project = 'Brain',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
# 
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('donor','batch','nCount_RNA'),npcs = 50,preprocess = TRUE)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 29,resolution = c(0.5,1,1.5),group.by = 'cell_type',label = TRUE)
# 
# saveRDS(macaque_RNA_seurat,file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
# 
# color_parameter <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
# temp_col <- color_parameter$col
# names(temp_col) <- color_parameter$id
# 
# p1 <- DimPlot(object = macaque_RNA_seurat,group.by = 'cell_type',cols = temp_col[unique(macaque_RNA_seurat$cell_type)],label = TRUE,repel = TRUE) + 
#   theme(aspect.ratio = 1)
# 
# p2 <- DimPlot(object = macaque_RNA_seurat,group.by = 'batch',cols = c('#d8b365','#5ab4ac'),label = FALSE) + 
#   theme(aspect.ratio = 1)
# 
# pdf(file = './res/step_12_fig_211021/macaque_RNA_seurat_merge_umap.pdf',width = 16,height = 8)
# p1+p2+plot_layout(ncol = 2)
# dev.off()

# integrate with PDhuman --------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')

macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('InPSB'))]
#convert gene id
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp

temp <- PDhuman_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno,filter_anno = FALSE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
PDhuman_RNA_seurat[['converted']] <- temp

rm(temp)
gc()

#find variable feature
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),rownames(PDhuman_RNA_seurat@assays$converted@counts))

# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('batch','donor','nCount_converted'),npcs = 50,preprocess = TRUE)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 40,resolution = 1,group.by = 'cell_type',label = TRUE)

macaque_RNA_seurat$species <- 'macaque'
PDhuman_RNA_seurat$species <- 'human'

#harmony
cl <- makeCluster(1)
clusterExport(cl,c('macaque_RNA_seurat','PDhuman_RNA_seurat','gene_list'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/','/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))
Brain_RNA_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_harmony_integration(named_seurat_list = list(macaque=macaque_RNA_seurat,human=PDhuman_RNA_seurat),assay = 'converted',variable_feature = gene_list,
                              var_to_regress_list = list(macaque=c('nCount_converted','donor','batch'),human=c('Donor','Library','nCount_converted')),npcs = 50,reference_loading = 'macaque',
                              integration_var = 'species',harmony_input_dim = 35,max.iter.harmony = 50,UMAP_dim = 35,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                              yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/')
  return(a)
})
stopCluster(cl)
Brain_RNA_seurat <- Brain_RNA_seurat[[1]]
gc()
my_send_sms('harmony done!')

DimPlot(object = Brain_RNA_seurat,group.by = 'dataset',label = FALSE)
Brain_RNA_seurat$cell_type <- NA
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)
Brain_RNA_seurat@meta.data[colnames(PDhuman_RNA_seurat),"cell_type"] <- as.character(PDhuman_RNA_seurat$Cluster)
DimPlot(object = Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')

#label transfer
cl <- makeCluster(1)
clusterExport(cl,c('Brain_RNA_seurat'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/',
                            '/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))

temp_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_MNN_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                             query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                             reference_var = 'cell_type',reduction = 'pca',
                             mnn = 100,knn = 300,Iteration = 5,return_query = TRUE)
  return(a)
})
stopCluster(cl)
gc()
temp_seurat <- temp_seurat[[1]]
my_send_sms(body = 'harmony mnn label transfer done!')

DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE)
Brain_RNA_seurat$predict_label <- as.character(Brain_RNA_seurat$cell_type)
Brain_RNA_seurat@meta.data[colnames(temp_seurat),"predict_label"] <- as.character(temp_seurat$predict_label)
DimPlot(Brain_RNA_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE,split.by = 'dataset')
DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE,cells.highlight = colnames(temp_seurat[,temp_seurat$predict_label == 'unknown']))
temp_seurat$cell_type <- temp_seurat$predict_label
temp_seurat@meta.data <- temp_seurat@meta.data[,1:10]

temp <- my_KNN_label_transfer(data = temp_seurat[,temp_seurat$cell_type != 'unknown'],query = temp_seurat[,temp_seurat$cell_type == 'unknown'],
                              reference_var = 'cell_type',reduction = 'pca',knn = 5,return_query = TRUE)
Brain_RNA_seurat@meta.data[colnames(temp),"predict_label"] <- as.character(temp$predict_label)
table(Brain_RNA_seurat$cell_type)
table(Brain_RNA_seurat$predict_label)

#umap
pdf(file = './res/step_12_fig_211021/macaque_RNA_PDhuman_RNA_harmony_umap.pdf',width = 16,height = 8)
DimPlot(Brain_RNA_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE,split.by = 'dataset') + 
  theme(aspect.ratio = 1)
dev.off()

#confusion heatmap
PDhuman_RNA_seurat$predict_label <- as.character(Brain_RNA_seurat@meta.data[colnames(PDhuman_RNA_seurat),"predict_label"])

pdf(file = './res/step_12_fig_211021/macaque_RNA_PDhuman_RNA_harmony_MNN_label_confusion_heatmap.pdf',width = 10,height = 10)
scibet::Confusion_heatmap(ori = PDhuman_RNA_seurat$Cluster,prd = PDhuman_RNA_seurat$predict_label) + 
  geom_text(aes(label = Prob),col='red') + 
  theme(aspect.ratio = 18/16,axis.title = element_text(size = 15)) + 
  xlab('ori') + ylab('prd label')
dev.off()

#marker dotplot
PDhuman_RNA_seurat <- my_process_seurat(object = PDhuman_RNA_seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

dotplot_matrix <- my_dotplot(PDhuman_RNA_seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Ast=c('AQP4','APOE'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2'),
                                             PSB=c('MEIS2','ETV1')),
                             group.by = 'predict_label', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_12_fig_211021/macaque_RNA_PDhuman_RNA_harmony_MNN_label_marker_dotplot.pdf',width = 16,height = 8)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

saveRDS(Brain_RNA_seurat,file = './processed_data/211014_summary/macaque_200919_210922_PDhuman_RNA_harmony_MNN_seurat_211029.rds')

# integrate with greenleaf------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
greenleaf_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')

greenleaf_RNA_seurat <- greenleaf_RNA_seurat[,greenleaf_RNA_seurat$Age != 'pcw24']
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('InPSB'))]

#convert gene id
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp

#find variable gene
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),rownames(greenleaf_RNA_seurat@assays$converted@counts))

# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('donor','batch','nCount_converted'),npcs = 50,preprocess = TRUE)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'cell_type',label = TRUE)

macaque_RNA_seurat$species <- 'macaque'
greenleaf_RNA_seurat$species <- 'human'

#harmony
cl <- makeCluster(1)
clusterExport(cl,c('macaque_RNA_seurat','greenleaf_RNA_seurat','gene_list'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/','/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))
Brain_RNA_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_harmony_integration(named_seurat_list = list(macaque=macaque_RNA_seurat,human=greenleaf_RNA_seurat),assay = 'converted',variable_feature = gene_list,
                              var_to_regress_list = list(macaque=c('nCount_converted','donor','batch'),human=c('Age','Batch','nCount_converted')),npcs = 50,reference_loading = 'macaque',
                              integration_var = 'species',harmony_input_dim = 25,max.iter.harmony = 50,UMAP_dim = 25,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                              yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/')
  return(a)
})
stopCluster(cl)
Brain_RNA_seurat <- Brain_RNA_seurat[[1]]
gc()
my_send_sms('harmony done!')

DimPlot(object = Brain_RNA_seurat,group.by = 'dataset',label = FALSE)
Brain_RNA_seurat$cell_type <- NA
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)
Brain_RNA_seurat@meta.data[colnames(greenleaf_RNA_seurat),"cell_type"] <- as.character(greenleaf_RNA_seurat$cell_type)
DimPlot(object = Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')

#label transfer
cl <- makeCluster(1)
clusterExport(cl,c('Brain_RNA_seurat'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/',
                            '/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))

temp_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_MNN_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                             query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                             reference_var = 'cell_type',reduction = 'pca',
                             mnn = 100,knn = 300,Iteration = 5,return_query = TRUE)
  return(a)
})
stopCluster(cl)
gc()
temp_seurat <- temp_seurat[[1]]
my_send_sms(body = 'harmony mnn label transfer done!')

DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE)
Brain_RNA_seurat$predict_label <- as.character(Brain_RNA_seurat$cell_type)
Brain_RNA_seurat@meta.data[colnames(temp_seurat),"predict_label"] <- as.character(temp_seurat$predict_label)
DimPlot(Brain_RNA_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE,split.by = 'dataset')
DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE,cells.highlight = colnames(temp_seurat[,temp_seurat$predict_label == 'unknown']))
temp_seurat$cell_type <- temp_seurat$predict_label

DimPlot(temp_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
temp_seurat@meta.data[temp_seurat$seurat_clusters %in% c('23','17'),"cell_type"] <- 'Unknown'
temp_seurat@meta.data <- temp_seurat@meta.data[,1:10]

temp <- my_KNN_label_transfer(data = temp_seurat[,temp_seurat$cell_type != 'unknown'],query = temp_seurat[,temp_seurat$cell_type == 'unknown'],
                              reference_var = 'cell_type',reduction = 'pca',knn = 5,return_query = TRUE)
Brain_RNA_seurat@meta.data[colnames(temp_seurat),"predict_label"] <- as.character(temp_seurat$cell_type)
Brain_RNA_seurat@meta.data[colnames(temp),"predict_label"] <- as.character(temp$predict_label)
table(Brain_RNA_seurat$cell_type)
table(Brain_RNA_seurat$predict_label)

#umap
pdf(file = './res/step_12_fig_211021/macaque_RNA_greenleaf_RNA_harmony_umap.pdf',width = 16,height = 8)
DimPlot(Brain_RNA_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE,split.by = 'dataset')
dev.off()

#confusion matrix
greenleaf_RNA_seurat$predict_label <- as.character(Brain_RNA_seurat@meta.data[colnames(greenleaf_RNA_seurat),"predict_label"])

pdf(file = './res/step_12_fig_211021/macaque_RNA_greenleaf_RNA_harmony_mnn_confusion_heatmap.pdf',width = 10,height = 8)
scibet::Confusion_heatmap(ori = greenleaf_RNA_seurat$cell_type,prd = greenleaf_RNA_seurat$predict_label) + 
  geom_text(aes(label = Prob),col='red') + 
  theme(aspect.ratio = 15/19,axis.title = element_text(size = 15)) + 
  xlab('ori') + ylab('prd label')
dev.off()

#marker plot
greenleaf_RNA_seurat <- my_process_seurat(object = greenleaf_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

dotplot_matrix <- my_dotplot(greenleaf_RNA_seurat,assay = 'converted', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Ast=c('AQP4','APOE'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'predict_label', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('Unknown','InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_12_fig_211021/macaque_RNA_greenleaf_RNA_harmony_MNN_label_marker_dotplot.pdf',width = 16,height = 8)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

temp_seurat@meta.data[colnames(greenleaf_RNA_seurat),'nCount'] <- as.numeric(greenleaf_RNA_seurat$nCount_originalexp)
temp_seurat@meta.data[colnames(greenleaf_RNA_seurat),'nFeature'] <- as.numeric(greenleaf_RNA_seurat$nFeature_originalexp)

pdf(file = './res/step_12_fig_211021/macaque_RNA_greenleaf_RNA_mnn_label_Unknown_nCount_nFeature_featureplot.pdf',width = 12,height = 6)
FeaturePlot(object = temp_seurat,features = c('nCount','nFeature'),reduction = 'umap')
dev.off()

p1 <- VlnPlot(greenleaf_RNA_seurat[,!(greenleaf_RNA_seurat$predict_label %in% c('End','Per','Mic','Astrocyte'))],features = c('nCount_originalexp'),group.by = 'predict_label',pt.size = 0) + ylim(c(0,30000)) + theme(aspect.ratio = 0.5)
p2 <- VlnPlot(greenleaf_RNA_seurat[,!(greenleaf_RNA_seurat$predict_label %in% c('End','Per','Mic','Astrocyte'))],features = c('nFeature_originalexp'),group.by = 'predict_label',pt.size = 0) + theme(aspect.ratio = 0.5)
pdf(file = './res/step_12_fig_211021/macaque_RNA_greenleaf_RNA_mnn_label_Unknown_nCount_nFeature_vlnplot.pdf',width = 16,height = 16)
p1+p2+plot_layout(ncol = 1)
dev.off()

Unknown_marker <- FindMarkers(object = greenleaf_RNA_seurat,ident.1 = 'Unknown',group.by = 'predict_label',assay = 'converted',slot = 'data',test.use = 'bimod',only.pos = TRUE)

pdf(file = './res/step_12_fig_211021/macaque_RNA_greenleaf_RNA_harmony_mnn_label_Unknown_marker_vlnplot.pdf',width = 18,height = 9)
VlnPlot(object = greenleaf_RNA_seurat,features = rownames(Unknown_marker)[1:12],group.by = 'predict_label',pt.size = 0)
dev.off()

saveRDS(Brain_RNA_seurat,file = './processed_data/211014_summary/macaque_200919_210922_greenleaf_RNA_harmony_MNN_seurat_211102.rds')
