#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: test harmony SNN pipeline                                       ##
## Data: 2021.09.23                                                                ##
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
library(rliger)
library(factoextra)
library(circlize)
library(ggsci)
library(SeuratWrappers)
library(magrittr)
library(harmony)
library(networkD3)
library(htmlwidgets)
library(SeuratData)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

#load data
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
DefaultAssay(PDhuman_RNA_seurat) <- 'RNA'
PDhuman_RNA_seurat <- RenameCells(object = PDhuman_RNA_seurat,new.names = paste('human',colnames(PDhuman_RNA_seurat),sep = '_'))
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$sub_cell_type %in% c('InPSB','Astrocyte'))]

color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

#pre-process
#convert gene name
dim(macaque_RNA_seurat)
dim(PDhuman_RNA_seurat)
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')
table(rownames(PDhuman_RNA_seurat) %in% c(human_to_human_anno[,1],human_to_human_anno[,2]))
temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
dim(temp)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp
temp <- PDhuman_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno,filter_anno = FALSE)
dim(temp)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
PDhuman_RNA_seurat[['converted']] <- temp

#find variable gene
macaque_RNA_seurat$dataset <- 'macaque'
PDhuman_RNA_seurat$dataset <- 'human'
#find variable gene
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',nfeatures = 6000,npcs = 50,preprocess = TRUE)
PDhuman_RNA_seurat <- my_process_seurat(object = PDhuman_RNA_seurat,assay = 'converted',nfeatures = 6000,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),VariableFeatures(PDhuman_RNA_seurat))

#integration
Brain_RNA_seurat <- my_harmony_integration(named_seurat_list = list(macaque=macaque_RNA_seurat,human=PDhuman_RNA_seurat),
                                           assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque=c('donor','nCount_converted'),human=c('Donor','Library','nCount_converted')),
                                           npcs = 50,reference_loading = 'macaque',integration_var = 'dataset',harmony_input_dim = 40,max.iter.harmony = 50,
                                           reference_dataset = NULL,UMAP_dim = 40,resolution = 1,kmeans_init_nstart = 10,kmeans_init_iter_max = 200,sigma = 0.3)

Brain_RNA_seurat@meta.data[,'cell_type'] <- NA
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$sub_cell_type
Brain_RNA_seurat@meta.data[colnames(PDhuman_RNA_seurat),"cell_type"] <- PDhuman_RNA_seurat$Cluster

# test KNN_n --------------------------------------------------------------
DimPlot(Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')

#set standard knn_k = 200, knn_n = 200
predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                         query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 200,knn_n = 200,n_core = 4,replace = TRUE,SNN = FALSE)
saveRDS(predicted_table,file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_200.rds')
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']
predicted_table <- predicted_table[colnames(temp),]
temp <- AddMetaData(object = temp,metadata = predicted_table[,1:2])

p1 <- DimPlot(temp,group.by = 'predicted_label',label = TRUE,repel = TRUE,cols = temp_col[unique(temp$predicted_label)]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))
pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_200_PDhuman_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

#validate by marker express
temp[['RNA']] <- CreateAssayObject(counts = PDhuman_RNA_seurat@assays$RNA@counts,min.cells = 0,min.features = 0)
dotplot_matrix <- my_dotplot(temp,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'predicted_label', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Mic','Per','End'))

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_200_PDhuman_marker_dotplot.pdf',width = 15,height = 6)
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

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_200_PDhuman_confusion_heatmap.pdf',width = 6,height = 6)
scibet::Confusion_heatmap(ori = temp$cell_type,prd = temp$predicted_label) + 
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 15)) + 
  xlab('original annotation') + ylab('predicted label')
dev.off()

#knn_n=10
predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                         query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 200,knn_n = 10,n_core = 4,replace = TRUE,SNN = FALSE)
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']
predicted_table <- predicted_table[colnames(temp),]
temp <- AddMetaData(object = temp,metadata = predicted_table[,1:2])

p1 <- DimPlot(temp,group.by = 'predicted_label',label = TRUE,repel = TRUE,cols = temp_col[unique(temp$predicted_label)]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_10_PDhuman_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

#knn_n=200,replace=FALSE
predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                         query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 200,knn_n = 200,n_core = 4,replace = FALSE,SNN = FALSE)
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']
predicted_table <- predicted_table[colnames(temp),]
temp <- AddMetaData(object = temp,metadata = predicted_table[,1:2])

p1 <- DimPlot(temp,group.by = 'predicted_label',label = TRUE,repel = TRUE,cols = temp_col[unique(temp$predicted_label)]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_200_no_replace_PDhuman_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

#test knn_n = 10,20,30,50,100,300,400,600,1000
test_n <- base::lapply(X = c(10,20,30,50,100,300,400,600,1000),FUN = function(knn_n){
  predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                           query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                           reference_var = 'cell_type',query_var = 'predicted_label',
                                           reduction = 'pca',knn_k = 200,knn_n = knn_n,n_core = 4,replace = TRUE,SNN = FALSE)
  ori <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_200.rds')
  predicted_table <- predicted_table[rownames(ori),]
  print(paste0('knn_n ',as.character(knn_n),' test done!'))
  return(sum(predicted_table$predicted_label == ori$predicted_label)/length(rownames(ori)))
})

test_n <- unlist(test_n)
test_n <- data.frame(knn_n=c(10,20,30,50,100,300,400,600,1000),accuracy=test_n)

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_test_n.pdf',width = 4,height = 3)
ggplot(test_n,aes(x=knn_n,y=accuracy))+
  geom_point(size=2,colour='red') + 
  geom_line() + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)
dev.off()



# test KNN_k --------------------------------------------------------------
test_k <- base::lapply(X = c(10,20,30,50,100,300,400,600,1000),FUN = function(knn_k){
  predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                           query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                           reference_var = 'cell_type',query_var = 'predicted_label',
                                           reduction = 'pca',knn_k = knn_k,knn_n = 200,n_core = 4,replace = TRUE,SNN = FALSE)
  ori <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_200_n_200.rds')
  predicted_table <- predicted_table[rownames(ori),]
  print(paste0('knn_k ',as.character(knn_k),' test done!'))
  return(sum(predicted_table$predicted_label == ori$predicted_label)/length(rownames(ori)))
})

test_k <- unlist(test_k)
test_k <- data.frame(knn_k=c(10,20,30,50,100,300,400,600,1000),accuracy=test_k)

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_n_200_test_k.pdf',width = 4,height = 3)
ggplot(test_k,aes(x=knn_k,y=accuracy))+
  geom_point(size=2,colour='red') + 
  geom_line() + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)
dev.off()

#knn_k = 10,knn_n = 200
predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                         query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 10,knn_n = 200,n_core = 4,replace = TRUE,SNN = FALSE)
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']
predicted_table <- predicted_table[colnames(temp),]
temp <- AddMetaData(object = temp,metadata = predicted_table[,1:2])

p1 <- DimPlot(temp,group.by = 'predicted_label',label = TRUE,repel = TRUE,cols = temp_col[unique(temp$predicted_label)]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_KNN_label_transfer_k_10_n_200_PDhuman_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()


# test SNN ---------------------------------------------------------------------
#knn_k = 500,knn_n = 200
predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                         query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 500,knn_n = 200,n_core = 4,replace = TRUE,SNN = TRUE)
saveRDS(predicted_table,file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_500_n_200.rds')
temp <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']
predicted_table <- predicted_table[colnames(temp),]
temp <- AddMetaData(object = temp,metadata = predicted_table[,1:2])

p1 <- DimPlot(temp,group.by = 'predicted_label',label = TRUE,repel = TRUE,cols = temp_col[unique(temp$predicted_label)]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_500_n_200_PDhuman_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

#test knn_k = c(20,100,300,800,1500)
for (knn_k in c(20,100,300,800,1500)) {
  predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                           query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                           reference_var = 'cell_type',query_var = 'predicted_label',
                                           reduction = 'pca',knn_k = knn_k,knn_n = 200,n_core = 4,replace = TRUE,SNN = TRUE)
  gc()
  char <- paste0('./res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_',as.character(knn_k),'_n_200.rds')
  saveRDS(predicted_table,file = char)
  print(paste0('knn_k ',as.character(knn_k),' test done!'))
}

knn_20 <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_20_n_200.rds')
knn_20 <- knn_20[colnames(PDhuman_RNA_seurat),]
knn_20$predicted_label <- paste(knn_20$predicted_label,'20',sep = '_')

knn_100 <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_100_n_200.rds')
knn_100 <- knn_100[colnames(PDhuman_RNA_seurat),]
knn_100$predicted_label <- paste(knn_100$predicted_label,'100',sep = '_')

knn_300 <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_300_n_200.rds')
knn_300 <- knn_300[colnames(PDhuman_RNA_seurat),]
knn_300$predicted_label <- paste(knn_300$predicted_label,'300',sep = '_')

knn_500 <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_500_n_200.rds')
knn_500 <- knn_500[colnames(PDhuman_RNA_seurat),]
knn_500$predicted_label <- paste(knn_500$predicted_label,'500',sep = '_')

knn_800 <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_800_n_200.rds')
knn_800 <- knn_800[colnames(PDhuman_RNA_seurat),]
knn_800$predicted_label <- paste(knn_800$predicted_label,'800',sep = '_')

knn_1500 <- readRDS(file = './res/step_10_fig_210923/macaque_PDhuman_SNN_label_transfer_k_1500_n_200.rds')
knn_1500 <- knn_1500[colnames(PDhuman_RNA_seurat),]
knn_1500$predicted_label <- paste(knn_1500$predicted_label,'1500',sep = '_')

confusion_matrix <- my_confusion_matrix(ori = knn_20$predicted_label,prd = knn_100$predicted_label)
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_100$predicted_label,prd = knn_300$predicted_label))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_300$predicted_label,prd = knn_500$predicted_label))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_500$predicted_label,prd = knn_800$predicted_label))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_800$predicted_label,prd = knn_1500$predicted_label))
#plot
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)
saveWidget(p1,file = './res/step_10_fig_210923/SNN_sankey_plot.html')


# public data validation --------------------------------------------------
#full data test
SeuratData::InstalledData()
data('panc8')
table(panc8$dataset)

panc8 <- my_process_seurat(object = panc8,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
panc8 <- my_process_seurat(object = panc8,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'celltype',label = TRUE)

p1 <- DimPlot(panc8,group.by = 'celltype',label = FALSE) + theme_classic() +theme(aspect.ratio = 1)
p2 <- DimPlot(panc8,group.by = 'tech',label = FALSE) + theme_classic() +theme(aspect.ratio = 1)
pdf(file = './res/step_10_fig_210923/panc8_raw_process_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

panc8_celseq <- panc8[,panc8$tech == 'celseq']
panc8_celseq2 <- panc8[,panc8$tech == 'celseq2']
panc8_fluidigmc1 <- panc8[,panc8$tech == 'fluidigmc1']
panc8_indrop <- panc8[,panc8$tech == 'indrop']
panc8_smartseq2 <- panc8[,panc8$tech == 'smartseq2']

panc8_smartseq2 <- my_process_seurat(object = panc8_smartseq2,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list_panc8 <- VariableFeatures(panc8_smartseq2)

test <- my_harmony_integration(named_seurat_list = list(celseq=panc8_celseq,celseq2=panc8_celseq2,fluidigmc1=panc8_fluidigmc1,
                                                        indrop=panc8_indrop,smartseq2=panc8_smartseq2),assay = 'RNA',
                               variable_feature = gene_list_panc8,
                               var_to_regress_list = list(celseq=NULL,celseq2=NULL,fluidigmc1=NULL,indrop=NULL,smartseq2=NULL),
                               npcs = 50,reference_loading = 'smartseq2',integration_var = 'tech',harmony_input_dim = 20,max.iter.harmony = 50,
                               UMAP_dim = 20,resolution = 1,kmeans_init_iter_max = 200)

DimPlot(test,group.by = 'dataset')

test@meta.data[,'cell_type'] <- NA
test@meta.data[,'cell_type'] <- as.character(panc8@meta.data[colnames(test),"celltype"])

p1 <- DimPlot(test,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_classic() +theme(aspect.ratio = 1)
p2 <- DimPlot(test,group.by = 'dataset',label = FALSE) + theme_classic() +theme(aspect.ratio = 1)

pdf(file = './res/step_10_fig_210923/panc8_harmony_integration_umap.pdf',width = 16,height = 8)
p1+p2 + plot_layout()
dev.off()

#missing data test
table(panc8[,panc8$celltype == 'activated_stellate']$tech)

#smartseq2 as reference
#smartseq2 missing stellate cell
#label transfer from smartseq2 to indrop
panc8_smartseq2 <- panc8_smartseq2[,!(panc8_smartseq2$celltype %in% c('activated_stellate','quiescent_stellate','schwann'))]

test <- my_harmony_integration(named_seurat_list = list(celseq=panc8_celseq,celseq2=panc8_celseq2,fluidigmc1=panc8_fluidigmc1,
                                                        indrop=panc8_indrop,smartseq2=panc8_smartseq2),assay = 'RNA',
                               variable_feature = gene_list_panc8,
                               var_to_regress_list = list(celseq=NULL,celseq2=NULL,fluidigmc1=NULL,indrop=NULL,smartseq2=NULL),
                               npcs = 50,reference_loading = 'smartseq2',integration_var = 'tech',harmony_input_dim = 20,max.iter.harmony = 50,
                               UMAP_dim = 20,resolution = 1,kmeans_init_iter_max = 200)

test@meta.data[,'cell_type'] <- NA
test@meta.data[,'cell_type'] <- as.character(panc8@meta.data[colnames(test),"celltype"])
test <- test[,test$dataset %in% c('smartseq2','indrop')]

p1 <- DimPlot(test,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_classic() +theme(aspect.ratio = 1)
p2 <- DimPlot(test,group.by = 'dataset',label = FALSE) + theme_classic() +theme(aspect.ratio = 1)

pdf(file = './res/step_10_fig_210923/panc8_without_stellate_harmony_integration_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

#SNN label transfer
predicted_table <- my_knn_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                         query = test[,test$dataset == 'indrop'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 50,knn_n = 50,n_core = 4,replace = TRUE,SNN = TRUE)

temp <- test[,test$dataset == 'indrop']
predicted_table <- predicted_table[colnames(temp),]
temp <- AddMetaData(object = temp,metadata = predicted_table[,c(1,2)])

p1 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_classic() +theme(aspect.ratio = 1)
p2 <- DimPlot(temp,group.by = 'predicted_label',label = TRUE,repel = TRUE) + theme_classic() +theme(aspect.ratio = 1)

pdf(file = './res/step_10_fig_210923/panc8_without_stellate_harmony_integration_SNN_label_transfer_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

#KNN label transfer
predicted_table <- my_knn_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                         query = test[,test$dataset == 'indrop'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 50,knn_n = 50,n_core = 4,replace = TRUE,SNN = FALSE)

temp <- test[,test$dataset == 'indrop']
predicted_table <- predicted_table[colnames(temp),]
temp <- AddMetaData(object = temp,metadata = predicted_table[,c(1,2)])

p1 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_classic() +theme(aspect.ratio = 1)
p2 <- DimPlot(temp,group.by = 'predicted_label',label = TRUE,repel = TRUE) + theme_classic() +theme(aspect.ratio = 1)

pdf(file = './res/step_10_fig_210923/panc8_without_stellate_harmony_integration_KNN_label_transfer_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()



# default harmony methods -------------------------------------------------
gene_list <- dplyr::intersect(rownames(macaque_RNA_seurat@assays$converted@counts),rownames(PDhuman_RNA_seurat@assays$converted@counts))
Brain_RNA_seurat <- cbind(macaque_RNA_seurat@assays$converted@counts[gene_list,],PDhuman_RNA_seurat@assays$converted@counts[gene_list,])
Brain_RNA_seurat <- CreateSeuratObject(counts = Brain_RNA_seurat,project = 'Brain',assay = 'RNA',min.cells = 0,min.features = 0)

Brain_RNA_seurat@meta.data[,'dataset'] <- NA
Brain_RNA_seurat@meta.data[,'cell_type'] <- NA

Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"dataset"] <- 'macaque'
Brain_RNA_seurat@meta.data[colnames(PDhuman_RNA_seurat),"dataset"] <- 'human'
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$sub_cell_type
Brain_RNA_seurat@meta.data[colnames(PDhuman_RNA_seurat),"cell_type"] <- PDhuman_RNA_seurat$Cluster

Brain_RNA_seurat <- my_process_seurat(object = Brain_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
Brain_RNA_seurat <- RunHarmony(object = Brain_RNA_seurat,group.by.vars = 'dataset')
Brain_RNA_seurat <- RunUMAP(object = Brain_RNA_seurat,reduction = 'harmony',dims = 1:40)

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_harmony_default_method_integration.pdf',width = 16,height = 8)
DimPlot(Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA,size = 1,colour = 'black'),
        axis.line = element_blank(),
        plot.title = element_blank())
dev.off()

pdf(file = './res/step_10_fig_210923/macaque_PDhuman_merge_pcaplot.pdf',width = 6,height = 6)
DimPlot(Brain_RNA_seurat,reduction = 'pca',dims = c(1,2),group.by = 'dataset',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1)
dev.off()
# integrate with greenleaf ------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$sub_cell_type %in% c('InPSB','Astrocyte'))]

greenleaf_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
DefaultAssay(greenleaf_RNA_seurat) <- 'converted'
# greenleaf_RNA_seurat <- greenleaf_RNA_seurat[,greenleaf_RNA_seurat$Age != 'pcw24']
#try integrate with whole greenleaf data
#convert data
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
macaque_RNA_seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

#find variable gene
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'pca',variable.feature = NULL,nfeatures = 5000,vars.to.regress = c('nCount_converted','donor'),npcs = 50,preprocess = TRUE)
greenleaf_RNA_seurat <- my_process_seurat(object = greenleaf_RNA_seurat,assay = 'converted',reduction.name = 'pca',variable.feature = NULL,nfeatures = 5000,vars.to.regress = c('nCount_converted','Batch'),npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),VariableFeatures(greenleaf_RNA_seurat))

#integration
macaque_RNA_seurat$dataset <- 'macaque'
greenleaf_RNA_seurat$dataset <- 'human'

Brain_RNA_seurat <- my_harmony_integration(named_seurat_list = list(macaque=macaque_RNA_seurat,human=greenleaf_RNA_seurat),
                                           assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque=c('nCount_converted','donor'),human=c('nCount_converted','Batch')),
                                           npcs = 50,reference_loading = 'macaque',integration_var = 'dataset',harmony_input_dim = 40,
                                           max.iter.harmony = 50,UMAP_dim = 40,resolution = 1,kmeans_init_iter_max = 200)

Brain_RNA_seurat@meta.data[,'cell_type'] <- NA
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$sub_cell_type)
Brain_RNA_seurat@meta.data[colnames(greenleaf_RNA_seurat),"cell_type"] <- greenleaf_RNA_seurat$cell_type

DimPlot(Brain_RNA_seurat,group.by = 'dataset',label = FALSE)
DimPlot(Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')
