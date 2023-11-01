#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: annotate cell types in 68 and 84                                ##
## Data: 2021.06.18                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_filted_express_matrix_210618.rds')

#preprocess
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',min.cells = 0,min.features = 0)

MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(macaque_RNA_seurat))
macaque_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_seurat,features = MT_gene_list)
VlnPlot(macaque_RNA_seurat,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)

macaque_RNA_seurat$species <- 'macaque'
donor_list <- rownames(macaque_RNA_seurat@meta.data)
donor_list <- base::lapply(donor_list,function(x){
  temp <- strsplit(x,split = '_',fixed = TRUE)
  return(temp[[1]][1])
})
donor_list <- unlist(donor_list)
macaque_RNA_seurat$donor <- donor_list
table(macaque_RNA_seurat$donor)

macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

#reduce dimension
use_dims <- 37
macaque_RNA_seurat <- FindNeighbors(macaque_RNA_seurat,dims = 1:use_dims)
macaque_RNA_seurat <- FindClusters(macaque_RNA_seurat,resolution = 1.6)
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:use_dims,reduction = 'pca')
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'seurat_clusters',label = TRUE,repel = TRUE)


#annotate
macaque_RNA_seurat$cell_type <- macaque_RNA_seurat$seurat_clusters
macaque_RNA_seurat@meta.data$cell_type <- as.character(macaque_RNA_seurat@meta.data$cell_type)
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
#Mic cluster 31
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('CX3CR1','AIF1'))
Mic_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = 'Mic',assay = 'RNA',slot = 'data',only.pos = TRUE,logfc.threshold = 0,min.pct = 0.2)
Mic_DEG <- Mic_DEG[order(Mic_DEG$pct.1,decreasing = TRUE),]
Mic_DEG <- Mic_DEG[order(Mic_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(Mic_DEG)[1:12])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters == '31',"cell_type"] <- 'Mic'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#Per cluster33
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('RGS5'))
Per_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = 'Per',assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
Per_DEG <- Per_DEG[order(Per_DEG$pct.1,decreasing = TRUE),]
Per_DEG <- Per_DEG[order(Per_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(Per_DEG)[1:12])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters == '33',"cell_type"] <- 'Per'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#End cluster28
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('CLDN5','ITM2A'))
End_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = 'End',assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
End_DEG <- End_DEG[order(End_DEG$pct.1,decreasing = TRUE),]
End_DEG <- End_DEG[order(End_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(End_DEG)[1:12])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters == '28',"cell_type"] <- 'End'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#OPC cluster25
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('ALDH1L1','APOE','SLC1A3','CST3','MBP','AIF1','P2RY12','C1QB','PDGFRA','OLIG2'))
OPC_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = 'OPC',assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
OPC_DEG <- OPC_DEG[order(OPC_DEG$pct.1,decreasing = TRUE),]
OPC_DEG <- OPC_DEG[order(OPC_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(OPC_DEG)[1:12])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters == '25',"cell_type"] <- 'OPC'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#In cluster0,10,11,27,29,3,32,5,7
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('CALB2','SST','DLX1','DLX2','GAD1','GAD2'))
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#InCGE cluster10 32 5 7
InCGE_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = 'InCGE',assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
InCGE_DEG <- InCGE_DEG[order(InCGE_DEG$pct.1,decreasing = TRUE),]
InCGE_DEG <- InCGE_DEG[order(InCGE_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('0','10','11','27','29','3','32','5','7')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(InCGE_DEG)[1:12])
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('0','10','11','27','29','3','32','5','7')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('HTR3A','PROX1','CXCL14'))
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('10','32','5','7'),"cell_type"] <- 'InCGE'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#InMGE cluster0 11 27 3
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('0','11','27','29','3','InCGE')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('NPY','SST','LHX6','NXPH1'))
InMGE_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = 'InMGE',assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
InMGE_DEG <- InMGE_DEG[order(InMGE_DEG$pct.1,decreasing = TRUE),]
InMGE_DEG <- InMGE_DEG[order(InMGE_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('0','11','27','29','3','InCGE')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(InMGE_DEG)[1:12])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('0','11','27','3'),"cell_type"] <- 'InMGE'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#RG cluster13 26
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('CRYAB','HOPX','SOX2','HES1','PAX6','VIM'))
RG_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('oRG','vRG'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
RG_DEG <- RG_DEG[order(RG_DEG$pct.1,decreasing = TRUE),]
RG_DEG <- RG_DEG[order(RG_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','26','34','23','21','13','17','30')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(RG_DEG)[13:24])

#oRG cluster13
oRG_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('oRG'),ident.2 = c('vRG','OPC','PgG2M','PgS','IP'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
oRG_DEG <- oRG_DEG[order(oRG_DEG$pct.2,decreasing = FALSE),]
oRG_DEG <- oRG_DEG[order(oRG_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','26','34','23','21','13','17','30')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(oRG_DEG)[13:24])
#donor validate
donor_cluster_proportion <- My_Cell_Proportion(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('13','26')],split.by = 'cell_type',group.by = 'donor')
ggplot(data = donor_cluster_proportion, mapping = aes(x = cell_type,y = Proportion,fill = donor))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of donor contribution to each cluster',fill = 'donor')+
  scale_fill_manual(values = c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF'))+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()
#liuyt validate
liuyt_meta <- readRDS(file = './processed_data/macaque_scRNA_seq_annotation_from_liuyt_210612.rds')
rownames(liuyt_meta) <- paste('A',rownames(liuyt_meta),'-1',sep = '')
table(colnames(macaque_RNA_seurat) %in% rownames(liuyt_meta))
cell_list <- dplyr::intersect(rownames(liuyt_meta),colnames(macaque_RNA_seurat))
cluster_representation <- data.frame(cell_id=cell_list,
                                     sunym_cluster=as.character(macaque_RNA_seurat@meta.data[cell_list,"seurat_clusters"]),
                                     liuyt_cluster=as.character(liuyt_meta[cell_list,"SCT_snn_res.1.8"]))
temp <- My_confusion_heatmap(ori = cluster_representation$liuyt_cluster,prd = cluster_representation$sunym_cluster,return_data = TRUE)
temp_matrix <- matrix(nrow = length(unique(temp$prd)),ncol = length(unique(temp$ori)))
colnames(temp_matrix) <- unique(temp$ori)
rownames(temp_matrix) <- unique(temp$prd)
for (i in 1:dim(temp_matrix)[1]) {
  for (j in 1:dim(temp_matrix)[2]) {
    temp_matrix[i,j] <- as.numeric(temp[temp$ori == colnames(temp_matrix)[j] & temp$prd == rownames(temp_matrix)[i],"Prob"])
  }
}
temp_matrix <- temp_matrix[names(sort(base::apply(temp_matrix,1,max),decreasing = TRUE)),
                           names(sort(base::apply(temp_matrix,2,max),decreasing = TRUE))]
Heatmap(temp_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_row_names = TRUE,show_column_names = TRUE,
        column_title = 'liuyt cluster',row_title = 'sunym cluster',
        height = 4.9,width = 5.2,
        heatmap_width = unit(8,'inches'),heatmap_height = unit(8,'inches'),
        name = 'Prop')
# cluster13 is oRG in liuyt
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('13'),"cell_type"] <- 'oRG'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#vRG
vRG_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('vRG'),ident.2 = c('oRG','OPC','PgG2M','PgS','IP'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
vRG_DEG <- vRG_DEG[order(vRG_DEG$pct.1,decreasing = TRUE),]
vRG_DEG <- vRG_DEG[order(vRG_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','26','34','23','21','oRG','17','30')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(vRG_DEG)[13:24])
#liuyt validate
liuyt_meta <- readRDS(file = './processed_data/macaque_scRNA_seq_annotation_from_liuyt_210612.rds')
rownames(liuyt_meta) <- paste('A',rownames(liuyt_meta),'-1',sep = '')
table(colnames(macaque_RNA_seurat) %in% rownames(liuyt_meta))
cell_list <- dplyr::intersect(rownames(liuyt_meta),colnames(macaque_RNA_seurat))
cluster_representation <- data.frame(cell_id=cell_list,
                                     sunym_cluster=as.character(macaque_RNA_seurat@meta.data[cell_list,"seurat_clusters"]),
                                     liuyt_cluster=as.character(liuyt_meta[cell_list,"SCT_snn_res.1.8"]))
temp <- My_confusion_heatmap(ori = cluster_representation$liuyt_cluster,prd = cluster_representation$sunym_cluster,return_data = TRUE)
temp_matrix <- matrix(nrow = length(unique(temp$prd)),ncol = length(unique(temp$ori)))
colnames(temp_matrix) <- unique(temp$ori)
rownames(temp_matrix) <- unique(temp$prd)
for (i in 1:dim(temp_matrix)[1]) {
  for (j in 1:dim(temp_matrix)[2]) {
    temp_matrix[i,j] <- as.numeric(temp[temp$ori == colnames(temp_matrix)[j] & temp$prd == rownames(temp_matrix)[i],"Prob"])
  }
}
temp_matrix <- temp_matrix[names(sort(base::apply(temp_matrix,1,max),decreasing = TRUE)),
                           names(sort(base::apply(temp_matrix,2,max),decreasing = TRUE))]
Heatmap(temp_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_row_names = TRUE,show_column_names = TRUE,
        column_title = 'liuyt cluster',row_title = 'sunym cluster',
        height = 4.9,width = 5.2,
        heatmap_width = unit(8,'inches'),heatmap_height = unit(8,'inches'),
        name = 'Prop')
# cluster26 is a part of cluster14 in liuyt meta
# cluster26 and cluster23 is cluster14 in liuyt meta

#cycling progenitor cluster21
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes[!(s.genes %in% rownames(macaque_RNA_seurat))]
#ENSMMUG00000004522
s.genes <- c(s.genes[s.genes %in% rownames(macaque_RNA_seurat)],'ENSMMUG00000004522')
table((s.genes %in% rownames(macaque_RNA_seurat)))
table((g2m.genes %in% rownames(macaque_RNA_seurat)))
g2m.genes <- g2m.genes[g2m.genes %in% rownames(macaque_RNA_seurat)]
macaque_RNA_seurat <- CellCycleScoring(macaque_RNA_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
p1 <- DimPlot(macaque_RNA_seurat,group.by = 'Phase',label = FALSE)
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2

Cyc_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('PgG2M','PgS'),ident.2 = c('oRG','OPC','vRG','IP'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
Cyc_DEG <- Cyc_DEG[order(Cyc_DEG$pct.1,decreasing = TRUE),]
Cyc_DEG <- Cyc_DEG[order(Cyc_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','26','34','23','21','oRG','17','30')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(Cyc_DEG)[13:24])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('21'),"cell_type"] <- 'Cycling'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('26','23'),"cell_type"] <- 'vRG'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#IP cluster17
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('oRG','Cycling','17','30','9','1')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('SOX2','EOMES','NEUROD2','TUBB3','NRP1'))
IP_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('IP'),ident.2 = c('oRG','vRG','PgG2M','PgS','ExN','ExM'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
IP_DEG <- IP_DEG[order(IP_DEG$pct.1,decreasing = TRUE),]
IP_DEG <- IP_DEG[order(IP_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('oRG','Cycling','17','30','9','1')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(IP_DEG)[13:24])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('17'),"cell_type"] <- 'IP'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

###############################################################
##take a rest and validate above results using greenleaf data##
###############################################################
#Cyc
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','34','Cycling','IP','30','9')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('TOP2A','MKI67','CLSPN','AURKA'))
#oRG
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','34','Cycling','IP','30','9')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('SOX9','HES1','ATP1A2','MOXD1','HOPX','FAM107A','MT3'))
#vRG
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$seurat_clusters %in% c('26','23')],group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',features = c('SOX9','HES1','ATP1A2','FBXO32','CNN1','CNN2'))
#WTF!
#tRG
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','34','Cycling','IP','30','9')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('CRYAB','NR4A1','FOXJ1'))
#OPC
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','34','Cycling','IP','30','9')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('SOX10','ENSMMUG00000056728','MBP'))
#ENSMMUG00000056728 refers to NKX2-2
#oIP
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','34','Cycling','IP','30','9')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('ASCL1','OLIG2','PDGFRA','EGFR'))
#Astrocyte
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','34','Cycling','IP','30','9','Astrocyte')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('AQP4','APOE','AGT'))
donor_cluster_proportion <- My_Cell_Proportion(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('34')],split.by = 'cell_type',group.by = 'donor')
ggplot(data = donor_cluster_proportion, mapping = aes(x = cell_type,y = Proportion,fill = donor))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of donor contribution to each cluster',fill = 'donor')+
  scale_fill_manual(values = c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF'))+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()
#cluster34 very likely to be Astrocyte
cluster_34_DEG <- FindMarkers(macaque_RNA_seurat,group.by = 'cell_type',ident.1 = c('34'),ident.2 = c('30','9','Cycling','IP','OPC','oRG','vRG'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('34'),"cell_type"] <- 'Astrocyte'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
#nIP
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','Astrocyte','Cycling','IP','30','9','1')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('EOMES','PPP1R17','PENK','NEUROG1','NEUROG2'))
#GluN cluster 30 9 1 4 12 22 2 14 8 16 15 24 6 18 20 19
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('OPC','vRG','oRG','Astrocyte','Cycling','IP','30','9','1','12','22','4','2','14','8','16','6','15','18','19','20','24','InCGE','InMGE')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('NEUROD2','TBR1','BCL11B','SATB2','SLC17A7'))
#In InCGE,InMGE,29
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$seurat_clusters %in% c('10','32','7','5','29','27','0','3','11')],group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',features = c('DLX2','DLX5','GAD2'))
#InCGE perfect
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$seurat_clusters %in% c('10','32','7','5','29','27','0','3','11')],group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',features = c('GAD2','SP8','NR2F2'))
#InMGE perfect
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$seurat_clusters %in% c('10','32','7','5','29','27','0','3','11')],group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',features = c('LHX6','SST','NKX2-1'))
#InPSB
#Transcriptomic and anatomic parcellation of 5-HT3AR expressing cortical interneuron subtypes revealed by single-cell RNA sequencing
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$seurat_clusters %in% c('10','32','7','5','29','27','0','3','11')],group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',features = c('MEIS2','PAX6','ETV1'))
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('29'),"cell_type"] <- 'InPSB'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
#Mic
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('AIF1','CCL3'))
#End perfect
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('CLDN5','PECAM1'))
#Per perfect
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('FOXC2','PDGFRB'))
#VLMC
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('FOXC2','COL1A1','LUM'))
#RBC
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('HEMGN'))


#The hardest part: ExN
#SP ExN cluster19 20
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('30','9','1','12','22','4','2','14','8','16','6','15','18','24','19','20')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('NR4A2','CRYM','ST18','CDH18','SOX5'))
ExDp_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('ExDp1','ExDp2'),ident.2 = c('ExN','ExM','ExM-U'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
ExDp_DEG <- ExDp_DEG[order(ExDp_DEG$pct.2,decreasing = FALSE),]
ExDp_DEG <- ExDp_DEG[order(ExDp_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('30','9','1','12','22','4','2','14','8','16','6','15','18','19','20','24')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(ExDp_DEG)[1:12])
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('19','20'),"cell_type"] <- 'Ex-SP'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
#FEZF2?

#Upper layer ExN
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('30','9','1','12','22','4','2','14','8','16','6','15','18','24','Ex-SP')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('FEZF2','TLE4','FOXP2','SOX5','TCERG1L','LDB2','CRYM','PCP4','NR4A2','TSHZ2'))
#TLE?
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('30','9','1','12','22','4','2','14','8','16','6','15','18','24','Ex-SP')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = c('SATB2','PTN','CUX1','CUX2','RORB','LMO4','LPL','LDB2','DKK3','PANTR1'))
ExM_U_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('ExM-U'),ident.2 = c('ExN','ExM','IP','ExDp1','ExDp2'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
ExM_U_DEG <- ExM_U_DEG[order(ExM_U_DEG$pct.2,decreasing = FALSE),]
ExM_U_DEG <- ExM_U_DEG[order(ExM_U_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('30','9','1','12','22','4','2','14','8','16','6','15','18','24','Ex-SP')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(ExM_U_DEG)[13:24])
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('30','9','1','12','22','4','2','14','8','16','6','15','18','24','Ex-SP')],group.by = 'cell_type',assay = 'RNA',slot = 'data',
        features = c('CUX2','POU3F2','LMO4','SATB2','TLE3','LHX2','MARCKSL1','ID2','MEF2C','TLE1','CUX1','RORB','RAC3','IGFBP4'))
#'NR4A2','LMO3','LDB2','ETV1','CRYM','FEZF2','SOX5','DKK3','TLE4','SEMA3E','TBR1'

#layer marker
macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[,"gene_name"][which(is.na(macaque_annotation$gene_name))] <- macaque_annotation[,"gene_id"][which(is.na(macaque_annotation$gene_name))]

layer_marker <- readRDS(file = './res/step_2_fig_210622/macaque_layer_DE_list_210622.rds')
CPo_marker <- layer_marker[['CPo']]
table(CPo_marker %in% macaque_annotation$gene_id)
CPo_marker <- CPo_marker[CPo_marker %in% macaque_annotation$gene_id]
CPo_marker <- macaque_annotation$gene_name[base::match(CPo_marker,table = macaque_annotation$gene_id)]
CPo_marker <- CPo_marker[CPo_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(CPo_marker))
geneSetAveragePlot(genes = CPo_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))


CPi_marker <- layer_marker[['CPi']]
table(CPi_marker %in% macaque_annotation$gene_id)
CPi_marker <- CPi_marker[CPi_marker %in% macaque_annotation$gene_id]
CPi_marker <- macaque_annotation$gene_name[base::match(CPi_marker,table = macaque_annotation$gene_id)]
CPi_marker <- CPi_marker[CPi_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(CPi_marker))
geneSetAveragePlot(genes = CPi_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))

SP_marker <- layer_marker[['SP']]
table(SP_marker %in% macaque_annotation$gene_id)
SP_marker <- SP_marker[SP_marker %in% macaque_annotation$gene_id]
SP_marker <- macaque_annotation$gene_name[base::match(SP_marker,table = macaque_annotation$gene_id)]
SP_marker <- SP_marker[SP_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(SP_marker))
geneSetAveragePlot(genes = SP_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))

IZ_marker <- layer_marker[['IZ']]
table(IZ_marker %in% macaque_annotation$gene_id)
IZ_marker <- IZ_marker[IZ_marker %in% macaque_annotation$gene_id]
IZ_marker <- macaque_annotation$gene_name[base::match(IZ_marker,table = macaque_annotation$gene_id)]
IZ_marker <- IZ_marker[IZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(IZ_marker))
geneSetAveragePlot(genes = IZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))

OSVZ_marker <- layer_marker[['OSVZ']]
table(OSVZ_marker %in% macaque_annotation$gene_id)
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% macaque_annotation$gene_id]
OSVZ_marker <- macaque_annotation$gene_name[base::match(OSVZ_marker,table = macaque_annotation$gene_id)]
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(OSVZ_marker))
geneSetAveragePlot(genes = OSVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))

ISVZ_marker <- layer_marker[['ISVZ']]
table(ISVZ_marker %in% macaque_annotation$gene_id)
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% macaque_annotation$gene_id]
ISVZ_marker <- macaque_annotation$gene_name[base::match(ISVZ_marker,table = macaque_annotation$gene_id)]
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(ISVZ_marker))
geneSetAveragePlot(genes = ISVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))

VZ_marker <- layer_marker[['VZ']]
table(VZ_marker %in% macaque_annotation$gene_id)
VZ_marker <- VZ_marker[VZ_marker %in% macaque_annotation$gene_id]
VZ_marker <- macaque_annotation$gene_name[base::match(VZ_marker,table = macaque_annotation$gene_id)]
VZ_marker <- VZ_marker[VZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(VZ_marker))
geneSetAveragePlot(genes = VZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))

MZ_marker <- layer_marker[['MZ']]
table(MZ_marker %in% macaque_annotation$gene_id)
MZ_marker <- MZ_marker[MZ_marker %in% macaque_annotation$gene_id]
MZ_marker <- macaque_annotation$gene_name[base::match(MZ_marker,table = macaque_annotation$gene_id)]
MZ_marker <- MZ_marker[MZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(MZ_marker))
geneSetAveragePlot(genes = MZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,object.class = 'matrix',plot.type = 'average',projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,reduction_key = 'UMAP',color.palette = c('grey','red'),trim = c(0,0.99))

#final annotation
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
#upper layer
ExM_U_DEG <- FindMarkers(PDhuman_RNA_seurat,group.by = 'Cluster',ident.1 = c('ExM-U'),ident.2 = c('ExN','ExM','ExDp1','ExDp2'),assay = 'RNA',slot = 'data',only.pos = TRUE,test.use = 'bimod',logfc.threshold = 0,min.pct = 0.2)
ExM_U_DEG <- ExM_U_DEG[order(ExM_U_DEG$pct.2,decreasing = FALSE),]
ExM_U_DEG <- ExM_U_DEG[order(ExM_U_DEG$p_val_adj,decreasing = FALSE),]
VlnPlot(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('30','9','1','12','22','4','2','14','8','16','6','15','18','24','Ex-SP')],group.by = 'cell_type',assay = 'RNA',slot = 'data',features = rownames(ExM_U_DEG)[13:24])
temp <- macaque_RNA_seurat
temp <- FindClusters(temp,resolution = 0.2)
DimPlot(temp,label = TRUE)
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('18','16','15','6','24'),"cell_type"] <- 'Ex-U'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('30','9','1','4','2','14','8','12','22'),"cell_type"] <- 'Ex'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
saveRDS(macaque_RNA_seurat,file = './processed_data/macaque_RNA_seurat_annotated_210623.rds')