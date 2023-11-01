#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: merge sample 50A and 82B                                        ##
## Data: 2021.09.27                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# load data ---------------------------------------------------------------
A50A_RNA_seurat <- Read10X(data.dir = './data/scRNA-seq/210922/A50A/outs/filtered_feature_bc_matrix/')
A82B_RNA_seurat <- Read10X(data.dir = './data/scRNA-seq/210922/A82B/outs/filtered_feature_bc_matrix/')

#create seurat object
A50A_RNA_seurat <- CreateSeuratObject(counts = A50A_RNA_seurat,project = 'A50A',assay = 'RNA',min.cells = 0,min.features = 0)
A82B_RNA_seurat <- CreateSeuratObject(counts = A82B_RNA_seurat,project = 'A82B',assay = 'RNA',min.cells = 0,min.features = 0)

# QC ----------------------------------------------------------------------
MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(A50A_RNA_seurat))

A50A_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = A50A_RNA_seurat,features = MT_gene_list)
A82B_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = A82B_RNA_seurat,features = MT_gene_list)

VlnPlot(A50A_RNA_seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(A50A_RNA_seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A50A_RNA_seurat <- subset(A50A_RNA_seurat,subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & nFeature_RNA > 200)

VlnPlot(A82B_RNA_seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(A82B_RNA_seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A82B_RNA_seurat <- subset(A82B_RNA_seurat,subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & nFeature_RNA > 200)

# merge -------------------------------------------------------------------
A50A_RNA_seurat <- RenameCells(object = A50A_RNA_seurat,new.names = paste('A50A',colnames(A50A_RNA_seurat),sep = '_'))
A82B_RNA_seurat <- RenameCells(object = A82B_RNA_seurat,new.names = paste('A82B',colnames(A82B_RNA_seurat),sep = '_'))

table(rownames(A50A_RNA_seurat) == rownames(A82B_RNA_seurat))

macaque_RNA_seurat <- cbind(A50A_RNA_seurat@assays$RNA@counts,A82B_RNA_seurat@assays$RNA@counts)
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_RNA_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_RNA_seurat$donor <- temp

# process seurat object ---------------------------------------------------
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

# filter doublet ----------------------------------------------------------
A50A_doublet <- read.csv(file = './res/step_11_fig_210926/A50A_scrublet_210926.csv')
A50A_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A50A_doublet$barcodes))
A50A_doublet$barcodes <- paste('A50A',as.character(A50A_doublet$barcodes),sep = '_')
table(duplicated(A50A_doublet$barcodes))
rownames(A50A_doublet) <- as.character(A50A_doublet$barcodes)
table(rownames(A50A_doublet) %in% colnames(macaque_RNA_seurat))

A82B_doublet <- read.csv(file = './res/step_11_fig_210926/A82B_scrublet_210926.csv')
A82B_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A82B_doublet$barcodes))
A82B_doublet$barcodes <- paste('A82B',as.character(A82B_doublet$barcodes),sep = '_')
table(duplicated(A82B_doublet$barcodes))
rownames(A82B_doublet) <- as.character(A82B_doublet$barcodes)
table(rownames(A82B_doublet) %in% colnames(macaque_RNA_seurat))

double_list <- rbind(A50A_doublet,A82B_doublet)
double_list <- double_list[,-1]
colnames(double_list) <- c('doublet_score','doublet')
double_list <- double_list[colnames(macaque_RNA_seurat),]
macaque_RNA_seurat@meta.data <- cbind(macaque_RNA_seurat@meta.data,double_list)

DimPlot(object = macaque_RNA_seurat,group.by = 'doublet',label = FALSE,cols = c('lightgrey','red'))

cell_proportion_table <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'doublet',split.by = 'seurat_clusters')
temp_col <- c('lightgrey','black')
names(temp_col) <- c('False','True')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = temp_col) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#cluster 20 and cluster 26
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
temp <- FindMarkers(object = macaque_RNA_seurat,ident.1 = '20',group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,min.pct = 0.5)
FeaturePlot(object = macaque_RNA_seurat,features = rownames(temp)[1:12])
VlnPlot(object = macaque_RNA_seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'seurat_clusters')

#filter
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('26','20'))]

# filter donor specific ---------------------------------------------------
p1 <- DimPlot(object = macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_RNA_seurat,group.by = 'donor',label = FALSE)
p1+p2

temp <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(temp,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)
#basically ok

# redo filter -------------------------------------------------------------
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_RNA_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_RNA_seurat$donor <- temp


macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = 1,group.by = 'seurat_clusters',label = TRUE)

#doublet
A50A_doublet <- read.csv(file = './res/step_11_fig_210926/A50A_scrublet_210926.csv')
A50A_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A50A_doublet$barcodes))
A50A_doublet$barcodes <- paste('A50A',as.character(A50A_doublet$barcodes),sep = '_')
table(duplicated(A50A_doublet$barcodes))
rownames(A50A_doublet) <- as.character(A50A_doublet$barcodes)
table(rownames(A50A_doublet) %in% colnames(macaque_RNA_seurat))

A82B_doublet <- read.csv(file = './res/step_11_fig_210926/A82B_scrublet_210926.csv')
A82B_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A82B_doublet$barcodes))
A82B_doublet$barcodes <- paste('A82B',as.character(A82B_doublet$barcodes),sep = '_')
table(duplicated(A82B_doublet$barcodes))
rownames(A82B_doublet) <- as.character(A82B_doublet$barcodes)
table(rownames(A82B_doublet) %in% colnames(macaque_RNA_seurat))

double_list <- rbind(A50A_doublet,A82B_doublet)
double_list <- double_list[,-1]
colnames(double_list) <- c('doublet_score','doublet')
double_list <- double_list[colnames(macaque_RNA_seurat),]
macaque_RNA_seurat@meta.data <- cbind(macaque_RNA_seurat@meta.data,double_list)

DimPlot(macaque_RNA_seurat,group.by = 'doublet',label = FALSE,cols = c('lightgrey','black'),pt.size = 0.5)
cell_proportion_table <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'doublet',split.by = 'seurat_clusters')
temp_col <- c('lightgrey','black')
names(temp_col) <- c('False','True')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = temp_col) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#donor specific
p1 <- DimPlot(object = macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_RNA_seurat,group.by = 'donor',label = FALSE)
p1+p2

cell_proportion_table <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#low quality
FeaturePlot(object = macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'))
VlnPlot(object = macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'),pt.size = 0,group.by = 'seurat_clusters')
#cluster 0 marker 
temp <- FindMarkers(object = macaque_RNA_seurat,ident.1 = '0',group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE)
VlnPlot(object = macaque_RNA_seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'seurat_clusters')
#filter clusters 0
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$seurat_clusters != '0']

# redo filter round 2 -----------------------------------------------------
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_RNA_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_RNA_seurat$donor <- temp

macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = 1,group.by = 'seurat_clusters',label = TRUE)

#doublet
A50A_doublet <- read.csv(file = './res/step_11_fig_210926/A50A_scrublet_210926.csv')
A50A_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A50A_doublet$barcodes))
A50A_doublet$barcodes <- paste('A50A',as.character(A50A_doublet$barcodes),sep = '_')
table(duplicated(A50A_doublet$barcodes))
rownames(A50A_doublet) <- as.character(A50A_doublet$barcodes)
table(rownames(A50A_doublet) %in% colnames(macaque_RNA_seurat))

A82B_doublet <- read.csv(file = './res/step_11_fig_210926/A82B_scrublet_210926.csv')
A82B_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A82B_doublet$barcodes))
A82B_doublet$barcodes <- paste('A82B',as.character(A82B_doublet$barcodes),sep = '_')
table(duplicated(A82B_doublet$barcodes))
rownames(A82B_doublet) <- as.character(A82B_doublet$barcodes)
table(rownames(A82B_doublet) %in% colnames(macaque_RNA_seurat))

double_list <- rbind(A50A_doublet,A82B_doublet)
double_list <- double_list[,-1]
colnames(double_list) <- c('doublet_score','doublet')
double_list <- double_list[colnames(macaque_RNA_seurat),]
macaque_RNA_seurat@meta.data <- cbind(macaque_RNA_seurat@meta.data,double_list)

DimPlot(macaque_RNA_seurat,group.by = 'doublet',label = FALSE,cols = c('lightgrey','black'),pt.size = 0.5)
cell_proportion_table <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'doublet',split.by = 'seurat_clusters')
temp_col <- c('lightgrey','black')
names(temp_col) <- c('False','True')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = temp_col) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#donor specific
p1 <- DimPlot(object = macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_RNA_seurat,group.by = 'donor',label = FALSE)
p1+p2

cell_proportion_table <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#cluster 25 26
FeaturePlot(object = macaque_RNA_seurat,features = c('PDGFRB','AQP4','APOE'))

#low quality
FeaturePlot(object = macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'))
VlnPlot(object = macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'),pt.size = 0,group.by = 'seurat_clusters')
#cluster 0 1 2
dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
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
                             group.by = 'seurat_clusters', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

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

#cluster0 1 2 are In and Ex
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
saveRDS(macaque_RNA_seurat,file = './res/step_11_fig_210926/macaque_RNA_210929_express_matrix_210929.rds')

# try merge with old macaque data -----------------------------------------
macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
DefaultAssay(macaque_old_seurat) <- 'RNA'

macaque_new_seurat <- readRDS(file = './res/step_11_fig_210926/macaque_RNA_210929_express_matrix_210929.rds')
macaque_new_seurat <- CreateSeuratObject(counts = macaque_new_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_new_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_new_seurat$donor <- temp

#preprocess
macaque_new_seurat <- my_process_seurat(object = macaque_new_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_new_seurat <- my_process_seurat(object = macaque_new_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 34,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

#merge
table(duplicated(colnames(macaque_new_seurat),colnames(macaque_old_seurat)))
table(rownames(macaque_new_seurat) == rownames(macaque_old_seurat))

macaque_RNA_seurat <- cbind(macaque_new_seurat@assays$RNA@counts,macaque_old_seurat@assays$RNA@counts)
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_RNA_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_RNA_seurat$donor <- temp

macaque_RNA_seurat@meta.data[,'batch'] <- NA
macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"batch"] <- '200919'
macaque_RNA_seurat@meta.data[colnames(macaque_new_seurat),"batch"] <- '210922'
table(macaque_RNA_seurat$batch)

#process
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 40,resolution = c(0.5,1,1.5),group.by = 'batch',label = FALSE)

macaque_RNA_seurat@meta.data[,'cell_type'] <- NA
macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"cell_type"] <- as.character(macaque_old_seurat$sub_cell_type)
table(macaque_RNA_seurat$cell_type)
#label transfer using MNN
cl <- makeCluster(1)
clusterExport(cl,c('macaque_RNA_seurat'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/',
                            '/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))

temp_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_MNN_label_transfer(data = macaque_RNA_seurat[,macaque_RNA_seurat$batch == '200919'],
                             query = macaque_RNA_seurat[,macaque_RNA_seurat$batch == '210922'],
                             reference_var = 'cell_type',reduction = 'pca',
                             mnn = 100,knn = 300,Iteration = 5,return_query = TRUE)
  return(a)
})

stopCluster(cl)
gc()

my_send_sms(body = 'merge mnn label transfer done!')
temp_seurat <- temp_seurat[[1]]
#label transfer done
DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE)

temp_seurat@meta.data[,'Clusters'] <- as.character(macaque_new_seurat@meta.data[colnames(temp_seurat),"seurat_clusters"])

p1 <- DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new predict label')
p2 <- DimPlot(temp_seurat,group.by = 'Clusters',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new Clusters')
p3 <- DimPlot(macaque_RNA_seurat,group.by = 'batch',label = FALSE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque batch')
p4 <- DimPlot(macaque_RNA_seurat[,colnames(macaque_old_seurat)],group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque old cell type')

pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_label_transfer.pdf',width = 16,height = 14)
p1+p2+p4+p3+plot_layout(ncol = 2)
dev.off()

#using cluster information to reannotate
macaque_new_seurat <- temp_seurat
rm(temp_seurat)
gc()

dotplot_matrix <- my_dotplot(macaque_new_seurat,assay = 'RNA', 
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
                             group.by = 'Clusters', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

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

#cluster 23 7 28 9 Unknown_InPSB
#cluster 24 21 17 Unknown_Ex_3
#cluster 13 26 0 InMGE
#cluster 3 11 InCGE
#cluster 20 InPSB
#cluster 15 10 Ex-4
#cluster 8 16 1 5 Ex-3
#cluster 18 Unknown_Ex

macaque_new_seurat$cell_type <- as.character(macaque_new_seurat$predict_label)
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('23','7','28','9'),"cell_type"] <- 'Unknown_InPSB'
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('24','21','17'),"cell_type"] <- 'Unknown_Ex-3'
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('13','26','0'),"cell_type"] <- 'InMGE'
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('3','11'),"cell_type"] <- 'InCGE'
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('20'),"cell_type"] <- 'InPSB'
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('15','10'),"cell_type"] <- 'Ex-4'
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('8','16','1','5'),"cell_type"] <- 'Ex-3'
macaque_new_seurat@meta.data[macaque_new_seurat$Clusters %in% c('18'),"cell_type"] <- 'Unknown_Ex'

pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_label_transfer_reannotation.pdf',width = 12,height = 8)
DimPlot(macaque_new_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new annotation')
dev.off()

#deal with unknown
pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_label_transfer_reannotation_highlight_unknown.pdf',width = 10,height = 10)
DimPlot(macaque_new_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,
        cells.highlight = colnames(macaque_new_seurat[,macaque_new_seurat$cell_type == 'unknown'])) + 
  theme(aspect.ratio = 1,legend.position = 'none') + labs(title = 'macaque new hightlight unknown')
dev.off()

temp_seurat <- my_KNN_label_transfer(data = macaque_new_seurat[,macaque_new_seurat$cell_type != 'unknown'],
                                     query = macaque_new_seurat[,macaque_new_seurat$cell_type == 'unknown'],
                                     reference_var = 'cell_type',reduction = 'pca',knn = 5,return_query = TRUE)

#marker gene express
dotplot_matrix <- my_dotplot(macaque_new_seurat,assay = 'RNA', 
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
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('unknown','Unknown_InPSB','InPSB','InCGE','InMGE','Ex-4','Unknown_Ex','Unknown_Ex-3','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_label_transfer_reannotation_marker_dotplot.pdf',width = 16,height = 8)
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


dotplot_matrix <- my_dotplot(macaque_old_seurat,assay = 'RNA', 
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
                             group.by = 'sub_cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InPSB','InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_11_fig_210926/macaque_old_RNA_seurat_label_transfer_reannotation_marker_dotplot.pdf',width = 16,height = 8)
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

#eliminate unknown
macaque_new_seurat@meta.data[colnames(temp_seurat),"cell_type"] <- as.character(temp_seurat$predict_label)
table(macaque_new_seurat$cell_type)
DimPlot(macaque_new_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

macaque_new_seurat$merge_cell_type <- macaque_new_seurat$cell_type

temp <- macaque_new_seurat@meta.data
temp$cell_id <- rownames(temp)
temp <- temp[,c('donor','batch','Clusters','predict_label','merge_cell_type','cell_id')]
saveRDS(temp,file = './res/step_11_fig_210926/macaque_RNA_210929_annotation_using_merge.rds')

#final fig
pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_label_transfer_reannotation_without_unknown.pdf',width = 12,height = 8)
DimPlot(macaque_new_seurat,group.by = 'merge_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new annotation')
dev.off()

dotplot_matrix <- my_dotplot(macaque_new_seurat,assay = 'RNA', 
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
                             group.by = 'merge_cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('Unknown_InPSB','InPSB','InCGE','InMGE','Ex-4','Unknown_Ex','Unknown_Ex-3','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_label_transfer_reannotation_without_unknown_marker_dotplot.pdf',width = 16,height = 8)
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
# try harmony integration -------------------------------------------------
macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
DefaultAssay(macaque_old_seurat) <- 'RNA'

macaque_new_seurat <- readRDS(file = './res/step_11_fig_210926/macaque_RNA_210929_express_matrix_210929.rds')
macaque_new_seurat <- CreateSeuratObject(counts = macaque_new_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_new_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_new_seurat$donor <- temp

gene_list <- VariableFeatures(macaque_old_seurat)
macaque_old_seurat$batch <- '200919'
macaque_new_seurat$batch <- '210922'

macaque_RNA_seurat <- my_harmony_integration(named_seurat_list = list(old=macaque_old_seurat,new=macaque_new_seurat),
                                             assay = 'RNA',variable_feature = gene_list,
                                             var_to_regress_list = list(old=c('nCount_RNA','donor'),new=c('nCount_RNA','donor')),
                                             npcs = 50,reference_loading = 'old',integration_var = 'batch',harmony_input_dim = 33,
                                             max.iter.harmony = 50,UMAP_dim = 33,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                             yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/',sigma = 0.8)

macaque_RNA_seurat@meta.data[,'cell_type'] <- NA
macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"cell_type"] <- as.character(macaque_old_seurat$sub_cell_type)
temp <- readRDS(file = './res/step_11_fig_210926/macaque_RNA_210929_annotation_using_merge.rds')
macaque_RNA_seurat@meta.data[rownames(temp),"cell_type"] <- as.character(temp$merge_cell_type)
DimPlot(macaque_RNA_seurat,group.by = 'dataset',label = FALSE)
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#label transfer
cl <- makeCluster(1)
clusterExport(cl,c('macaque_RNA_seurat'),envir = environment())
clusterEvalQ(cl,.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.1/',
                            '/data/User/sunym/software/R_lib/R_4.1.1/')))
clusterEvalQ(cl,source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R'))

temp_seurat <- parallel::parLapply(cl = cl,X = c(1),fun = function(x){
  a <- my_MNN_label_transfer(data = macaque_RNA_seurat[,macaque_RNA_seurat$dataset == '200919'],
                             query = macaque_RNA_seurat[,macaque_RNA_seurat$dataset == '210922'],
                             reference_var = 'cell_type',reduction = 'pca',
                             mnn = 100,knn = 300,Iteration = 5,return_query = TRUE)
  return(a)
})
stopCluster(cl)
gc()
temp_seurat <- temp_seurat[[1]]
my_send_sms(body = 'harmony mnn label transfer done!')
#label transfer done!
DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE)

p1 <- DimPlot(macaque_RNA_seurat[,colnames(macaque_old_seurat)],group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque old seurat cell type')

p2 <- DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new seurat predict label')

p3 <- DimPlot(temp_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'Integrated macaque data')

p1+p2+p3+plot_layout()

#modify annotation
macaque_new_seurat@meta.data[colnames(temp_seurat),'Clusters'] <- as.character(temp_seurat$seurat_clusters)
macaque_new_seurat <- NormalizeData(macaque_new_seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)

dotplot_matrix <- my_dotplot(macaque_new_seurat,assay = 'RNA', 
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
                             group.by = 'Clusters', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

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

#cluster 21 InPSB
#cluster 4 31 Unknown_InPSB
#cluster 9 Unknown_Ex-3
#cluster 8 22 5 6 18 Ex-3
#cluster 11 20 Ex-4
temp_seurat$cell_type <- as.character(temp_seurat$predict_label)
temp_seurat@meta.data[temp_seurat$seurat_clusters %in% c('21'),"cell_type"] <- 'InPSB'
temp_seurat@meta.data[temp_seurat$seurat_clusters %in% c('4','31'),"cell_type"] <- 'Unknown_InPSB'
temp_seurat@meta.data[temp_seurat$seurat_clusters %in% c('9'),"cell_type"] <- 'Unknown_Ex-3'
temp_seurat@meta.data[temp_seurat$seurat_clusters %in% c('8','22','5','6','18'),"cell_type"] <- 'Ex-3'
temp_seurat@meta.data[temp_seurat$seurat_clusters %in% c('11','20'),"cell_type"] <- 'Ex-4'
table(temp_seurat$cell_type)

#eliminate unknown
DimPlot(temp_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cells.highlight = colnames(temp_seurat[,temp_seurat$cell_type == 'unknown']))
temp <- my_KNN_label_transfer(data = temp_seurat[,temp_seurat$cell_type != 'unknown'],
                              query = temp_seurat[,temp_seurat$cell_type == 'unknown'],
                              reference_var = 'cell_type',reduction = 'pca',knn = 5,return_query = TRUE)
temp_seurat@meta.data[colnames(temp),"cell_type"] <- as.character(temp$predict_label)
table(temp_seurat$cell_type)
rm(temp)
gc()

#plot
p1 <- DimPlot(macaque_RNA_seurat[,colnames(macaque_old_seurat)],group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque old seurat cell type')

p2 <- DimPlot(temp_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new seurat predict label')

p3 <- DimPlot(temp_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new seurat clusters')

p4 <- DimPlot(temp_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque new seurat cell type')

pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_harmony_integration_label_transfer.pdf',width = 16,height = 14)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()


#marker gene express
macaque_new_seurat@meta.data[colnames(temp_seurat),'cell_type'] <- as.character(temp_seurat$cell_type)

dotplot_matrix <- my_dotplot(macaque_new_seurat,assay = 'RNA', 
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
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('Unknown_InPSB','InPSB','InCGE','InMGE','Ex-4','Unknown_Ex-3','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_11_fig_210926/macaque_RNA_210929_harmony_label_transfer_reannotation_marker_dotplot.pdf',width = 16,height = 8)
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

macaque_new_seurat@meta.data[colnames(temp_seurat),'predict_label'] <- as.character(temp_seurat$predict_label)
temp <- macaque_new_seurat@meta.data
temp$cell_id <- rownames(temp)
saveRDS(temp,file = './res/step_11_fig_210926/macaque_RNA_210929_annotation_using_harmony.rds')

# final annotation --------------------------------------------------------
#load data
macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
DefaultAssay(macaque_old_seurat) <- 'RNA'

macaque_new_seurat <- readRDS(file = './res/step_11_fig_210926/macaque_RNA_210929_express_matrix_210929.rds')
macaque_new_seurat <- CreateSeuratObject(counts = macaque_new_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_new_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_new_seurat$donor <- temp

#load annotation
merge_anno <- readRDS(file = './res/step_11_fig_210926/macaque_RNA_210929_annotation_using_merge.rds')
merge_anno <- merge_anno[colnames(macaque_new_seurat),]
harmony_anno <- readRDS(file = './res/step_11_fig_210926/macaque_RNA_210929_annotation_using_harmony.rds')
harmony_anno <- harmony_anno[colnames(macaque_new_seurat),]

macaque_new_seurat$merge_anno <- merge_anno$merge_cell_type
macaque_new_seurat$harmony_anno <- harmony_anno$cell_type

#compare merge anno and harmony anno
pdf(file = './res/step_11_fig_210926/macaque_new_seurat_merge_to_harmony_annotation_confusion_heatmap.pdf',width = 12,height = 10)
scibet::Confusion_heatmap(ori = macaque_new_seurat$merge_anno,prd = macaque_new_seurat$harmony_anno) + 
  geom_text(aes(label = Prob),col='red') + 
  theme(aspect.ratio = 20/21,axis.title = element_text(size = 15)) + 
  xlab('annotation using merge') + ylab('annotation using harmony')
dev.off()

confusion_matrix <- my_confusion_matrix(ori = paste(macaque_new_seurat$merge_anno,'merge',sep = '_'),prd = paste(macaque_new_seurat$harmony_anno,'harmony',sep = '_'))
nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

saveWidget(p1,file = './res/step_11_fig_210926/macaque_RNA_210929_annotation_merge_to_harmony.html')

#using merge annotation as final annotation
macaque_new_seurat <- my_process_seurat(object = macaque_new_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_new_seurat <- my_process_seurat(object = macaque_new_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 34,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

#get merged UMAP
table(duplicated(colnames(macaque_new_seurat),colnames(macaque_old_seurat)))
table(rownames(macaque_new_seurat) == rownames(macaque_old_seurat))
macaque_RNA_seurat <- cbind(macaque_new_seurat@assays$RNA@counts,macaque_old_seurat@assays$RNA@counts)
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(colnames(macaque_RNA_seurat),function(c){
  return(strsplit(x = c,split = '_',fixed = TRUE)[[1]][1])
})
temp <- unlist(temp)
macaque_RNA_seurat$donor <- temp
macaque_RNA_seurat@meta.data[,'batch'] <- NA
macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"batch"] <- '200919'
macaque_RNA_seurat@meta.data[colnames(macaque_new_seurat),"batch"] <- '210922'
table(macaque_RNA_seurat$batch)
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 40,resolution = c(0.5,1,1.5),group.by = 'batch',label = FALSE)
my_send_sms(body = 'UMAP done!')

temp <- macaque_RNA_seurat[,colnames(macaque_new_seurat)]
temp <- temp@reductions$umap

macaque_new_seurat@reductions$UMAP <- temp

#save data
saveRDS(macaque_new_seurat,file = './processed_data/macaque_RNA_210929_seurat_annotated_211013.rds')

#Unknown_Ex-3 corresponding cell in macaque_old_seurat
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
temp <- macaque_RNA_seurat[,colnames(macaque_old_seurat)]
temp <- temp[,temp$seurat_clusters %in% c('20','24')]

macaque_old_seurat$new_cell_type <- as.character(macaque_old_seurat$sub_cell_type)
macaque_old_seurat@meta.data[colnames(temp),"new_cell_type"] <- 'Unknown_Ex-3'

pdf(file = './res/step_11_fig_210926/macaque_old__seurat_show_unknown_Ex_3_umap.pdf',width = 8,height = 6)
DimPlot(macaque_old_seurat,group.by = 'new_cell_type',label = TRUE,repel = TRUE)
dev.off()
temp <- macaque_old_seurat@meta.data
saveRDS(temp,file = './res/step_11_fig_210926/macaque_old_seurat_meta_data_with_unkonwn_Ex_3.rds')

# new cell type in macaque new seurat -------------------------------------
#load data
macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
DefaultAssay(macaque_old_seurat) <- 'RNA'
macaque_new_seurat <- readRDS(file = './processed_data/macaque_RNA_210929_seurat_annotated_211013.rds')
DefaultAssay(macaque_new_seurat) <- 'RNA'

#find all markers
Idents(macaque_new_seurat) <- 'merge_anno'
macaque_new_marker <- FindAllMarkers(object = macaque_new_seurat,assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,
                                     logfc.threshold = 0.25,min.pct = 0.1)
Idents(macaque_old_seurat) <- 'sub_cell_type'
macaque_old_marker <- FindAllMarkers(object = macaque_old_seurat,assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,
                                     logfc.threshold = 0.25,min.pct = 0.1)
my_send_sms(body = 'macaque new and old seurat data find markers done!')

table(macaque_new_marker$cluster)
table(macaque_old_marker$cluster)

#intersect proportion
new_list <- unique(as.character(macaque_new_seurat$merge_anno))
old_list <- unique(as.character(macaque_old_seurat$sub_cell_type))

temp <- do.call(rbind,base::lapply(X = new_list,function(x){
  temp_list <- base::lapply(old_list,function(y){
    gene_list_a <- macaque_new_marker[macaque_new_marker$cluster == x,"gene"]
    gene_list_b <- macaque_old_marker[macaque_old_marker$cluster == y,"gene"]
    gene_list_overlap <- dplyr::intersect(gene_list_a,gene_list_b)
    return(length(gene_list_overlap)/(length(gene_list_a)+length(gene_list_b)-length(gene_list_overlap)))
  })
  temp_list <- unlist(temp_list)
  names(temp_list) <- old_list
  return(temp_list)
}))

colnames(temp) <- old_list
rownames(temp) <- new_list
temp <- as.matrix(temp)

#plot

#cluster
col_cluste_modified <- readRDS(file = './res/step_11_fig_210926/macaque_old_new_marker_overlap_col_cluste_modified.rds')
row_cluste_modified <- readRDS(file = './res/step_11_fig_210926/macaque_old_new_marker_overlap_row_cluste_modified.rds')

col_fun = colorRamp2(c(0,0.3,0.6), c('#3B4992FF','white','#EE0000FF'))

pdf(file = './res/step_11_fig_210926/macaque_old_new_marker_overlap_jaccard_index_hratmap.pdf',width = 8,height = 7)
Heatmap(temp,cluster_columns = col_cluste_modified,cluster_rows = row_cluste_modified,show_row_names = TRUE,show_column_names = TRUE,
        height = unit(5.25,'inches'),width = unit(4.5,'inches'),name = 'Jaccard index',col = col_fun)
dev.off()


## unknown Ex-3 marker -----------------------------------------------------
#anno
temp <- readRDS(file = './res/step_11_fig_210926/macaque_old_seurat_meta_data_with_unkonwn_Ex_3.rds')
macaque_old_seurat$new_cell_type <- temp[colnames(macaque_old_seurat),"new_cell_type"]

Unknown_Ex_3_marker <- FindMarkers(object = macaque_new_seurat,ident.1 = 'Unknown_Ex-3',group.by = 'merge_anno',assay = 'RNA',slot = 'data',
                                   test.use = 'bimod',logfc.threshold = 0.25,only.pos = TRUE,min.diff.pct = 0.2)

Unknown_Ex_3_marker_old <- FindMarkers(object = macaque_old_seurat,ident.1 = 'Unknown_Ex-3',group.by = 'new_cell_type',assay = 'RNA',slot = 'data',
                                       test.use = 'bimod',logfc.threshold = 0.25,only.pos = TRUE,min.diff.pct = 0.2)


gene_list <- dplyr::intersect(rownames(Unknown_Ex_3_marker),rownames(Unknown_Ex_3_marker_old))

p1 <- my_dotplot(object = macaque_new_seurat,assay = 'RNA',features = gene_list,cols = c('#3B4992FF','white','#EE0000FF'),
                 return_data_plot = FALSE,col.min = -2.5,col.max = 2.5,scale = TRUE,group.by = 'merge_anno') + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        legend.position = 'right',
        panel.grid = element_line(color="grey",size = 0.1),
        aspect.ratio = 0.5) + 
  xlab('') + ylab('macaque new data')

p2 <- my_dotplot(object = macaque_old_seurat,assay = 'RNA',features = gene_list,cols = c('#3B4992FF','white','#EE0000FF'),
                 return_data_plot = FALSE,col.min = -2.5,col.max = 2.5,scale = TRUE,group.by = 'new_cell_type') + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        legend.position = 'right',
        panel.grid = element_line(color="grey",size = 0.1),
        aspect.ratio = 0.5) + 
  xlab('') + ylab('macaque old data')

pdf(file = './res/step_11_fig_210926/macaque_old_new_unknow_Ex_3_marker_dotplot.pdf',width = 12,height = 12)
p1+p2+plot_layout(ncol = 1)
dev.off()

p1 <- geneSetAveragePlot(genes = c('SLIT3','DGKG','ADAMTS3'),object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panel',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = TRUE,plot.title = 'macaque new data')

p2 <- geneSetAveragePlot(genes = c('SLIT3','DGKG','ADAMTS3'),object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panel',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = TRUE,plot.title = 'macaque old data')

pdf(file = './res/step_11_fig_210926/macaque_old_new_unknow_Ex_3_marker_featureplot.pdf')
p1+p2+plot_layout(ncol = 1)
dev.off()


## Unknown Ex ---------------------------------------------------------------
Unknown_Ex_marker <- FindMarkers(object = macaque_new_seurat,ident.1 = 'Unknown_Ex',group.by = 'merge_anno',assay = 'RNA',slot = 'data',
                                 test.use = 'bimod',logfc.threshold = 1,min.diff.pct = 0.2,only.pos = TRUE,min.pct = 0.5)

pdf(file = './res/step_11_fig_210926/macaque_new_seurat_unknown_Ex_marker_vlnplot.pdf',width = 24,height = 10)
VlnPlot(macaque_new_seurat,features = rownames(Unknown_Ex_marker)[1:12],pt.size = 0,group.by = 'merge_anno')
dev.off()

pdf(file = './res/step_11_fig_210926/macaque_new_seurat_unknown_Ex_marker_dotplot.pdf',width = 14,height = 7)
my_dotplot(object = macaque_new_seurat,assay = 'RNA',features = rownames(Unknown_Ex_marker),cols = c('#3B4992FF','white','#EE0000FF'),
           return_data_plot = FALSE,col.min = -2.5,col.max = 2.5,scale = TRUE,group.by = 'merge_anno') + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        legend.position = 'right',
        panel.grid = element_line(color="grey",size = 0.1),
        aspect.ratio = 0.5) + 
  xlab('') + ylab('macaque new data')
dev.off()

p1 <- geneSetAveragePlot(genes = c('TLL1','NR4A2','B3GAT2','NTNG2'),object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panel',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = TRUE,plot.title = 'macaque new data')

p2 <- geneSetAveragePlot(genes = c('TLL1','NR4A2','B3GAT2','NTNG2'),object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panel',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = TRUE,plot.title = 'macaque old data')

pdf(file = './res/step_11_fig_210926/macaque_old_new_unknow_Ex_featureplot.pdf',width = 18,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

## Unknown_InPSB -----------------------------------------------------------
Unknown_InPSB_marker <- FindMarkers(object = macaque_new_seurat,ident.1 = 'Unknown_InPSB',group.by = 'merge_anno',assay = 'RNA',slot = 'data',
                                    test.use = 'bimod',logfc.threshold = 0.25,min.pct = 0.5,min.diff.pct = 0.4,only.pos = TRUE)

pdf(file = './res/step_11_fig_210926/macaque_new_seurat_unknown_InPSB_vlnplot.pdf',width = 24,height = 12)
VlnPlot(macaque_new_seurat,features = rownames(Unknown_InPSB_marker),pt.size = 0)
dev.off()

p1 <- geneSetAveragePlot(genes = c('RXRG','RARB','SV2C','RGS9'),object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'UMAP',reduction_key = 'UMAP',plot.type = 'panel',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = TRUE,plot.title = 'macaque new data')

p2 <- geneSetAveragePlot(genes = c('RXRG','RARB','SV2C','RGS9'),object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panel',scaled = FALSE,color.palette = c('lightgrey','blue'),
                         aspectratio = 1,trim = NULL,print = TRUE,plot.title = 'macaque old data')

pdf(file = './res/step_11_fig_210926/macaque_old_new_unknow_InPSB_featureplot.pdf',width = 18,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()