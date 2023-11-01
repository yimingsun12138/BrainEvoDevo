#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: re_cluster macaque RNA_seq data                                 ##
## Data: 2022.06.26                                                                ##
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

# load macaque multiome data ----------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']

#recreate macaque multiome data
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech",'donor',"cell_type")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

# process of macaque multiome data ----------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,
                                             nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)

ndims <- 31
macaque_multiome_Seurat <- RunUMAP(object = macaque_multiome_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

macaque_multiome_Seurat <- FindNeighbors(object = macaque_multiome_Seurat,reduction = 'pca',dims = 1:ndims)
macaque_multiome_Seurat <- FindClusters(object = macaque_multiome_Seurat,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1))

#show clustree
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_clustree.pdf',width = 20,height = 10)
clustree(x = macaque_multiome_Seurat@meta.data,prefix = 'RNA_snn_res.')
dev.off()

for (i in c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1)) {
  col_name <- paste('RNA_snn_res.',i,sep = '')
  p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
  p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = col_name,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
  char <- paste0('./res/step_41_fig_220626/macaque_multiome_Seurat_',col_name,'_dimplot.pdf')
  pdf(file = char,width = 16,height = 7)
  print(p2+p1+plot_layout(ncol = 2))
  dev.off()
}

#show donor contribution
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'donor',label = FALSE) + theme(aspect.ratio = 1)
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_cell_type_donor_dimplot.pdf',width = 16,height = 7)
p1+p2+plot_layout(ncol = 2)
dev.off()

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_cell_type_split_by_donor_dimplot.pdf',width = 10,height = 10)
DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',split.by = 'donor',ncol = 2,label = FALSE) + theme(aspect.ratio = 1)
dev.off()

#save meta data
meta_data <- macaque_multiome_Seurat@meta.data
saveRDS(meta_data,file = './res/step_41_fig_220626/macaque_multiome_Seurat_PCA_clustering_meta_data.rds')

# refine macaque multiome data annotation ---------------------------------

# using scATAC cluster to facilicate annotation
macaque_multiome_ArchR <- loadArchRProject(path = './ArchR/processed_data/macaque_multiome_ArchR_220411')
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_Seurat_PCA_clustering_meta_data.rds')
meta_data <- meta_data[rownames(macaque_multiome_ArchR@cellColData),]

for (i in c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1)) {
  col_name <- paste('RNA_snn_res.',i,sep = '')
  p1 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = meta_data,group.by = 'cell_type',split.by = NULL,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
  p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = meta_data,group.by = col_name,split.by = NULL,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
  p3 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),group.by = 'Clusters',split.by = NULL,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
  char <- paste0('./res/step_41_fig_220626/macaque_multiome_ArchR_',col_name,'_dimplot.pdf')
  pdf(file = char,width = 24,height = 7)
  print(p2+p3+p1+plot_layout(ncol = 3))
  dev.off()
}

my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = meta_data,group.by = 'donor',split.by = NULL,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)

# #add cluster(run on 202.205.131.32)
# for (i in c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1)) {
#   #cluster
#   macaque_multiome_ArchR <- addClusters(
#     input = macaque_multiome_ArchR,
#     reducedDims = "IterativeLSI",
#     method = "Seurat",
#     name = paste0('ATAC_snn_res.',i),
#     resolution = i,
#     maxClusters = 50,
#     dimsToUse = 1:16,
#     force = TRUE
#   )
# }
# meta_data <- macaque_multiome_ArchR@cellColData
# saveRDS(meta_data,file = './res/step_41_fig_220626/macaque_multiome_ArchR_LSI_clustering_meta_data.rds')

meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_ArchR_LSI_clustering_meta_data.rds')
meta_data <- as.data.frame(meta_data)
for (i in c(0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,2.1)) {
  col_name <- paste('ATAC_snn_res.',i,sep = '')
  p1 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = meta_data,group.by = 'cell_type',split.by = NULL,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
  p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = meta_data,group.by = col_name,split.by = NULL,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
  char <- paste0('./res/step_41_fig_220626/macaque_multiome_ArchR_',col_name,'_dimplot.pdf')
  pdf(file = char,width = 16,height = 7)
  print(p2+p1+plot_layout(ncol = 2))
  dev.off()
}


## check Ex-2 sub-cluster marker -------------------------------------------
cluster_0_marker_resolution_0.3 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '0',group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)
temp <- cluster_0_marker_resolution_0.3[cluster_0_marker_resolution_0.3$pct.1 > 0.4 & cluster_0_marker_resolution_0.3$pct.2 < 0.2 & cluster_0_marker_resolution_0.3$avg_log2FC > 0.5,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data')

cluster_7_marker_resolution_0.3 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '7',group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)
temp <- cluster_7_marker_resolution_0.3[cluster_7_marker_resolution_0.3$pct.1 > 0.4 & cluster_7_marker_resolution_0.3$pct.2 < 0.2 & cluster_7_marker_resolution_0.3$avg_log2FC > 0.5,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data')

ggvenn(data = list(cluster7 = rownames(cluster_7_marker_resolution_0.3),cluster0 = rownames(cluster_0_marker_resolution_0.3)),
       c('cluster7','cluster0'))

temp <- rownames(cluster_0_marker_resolution_0.3)[!(rownames(cluster_0_marker_resolution_0.3) %in% rownames(cluster_7_marker_resolution_0.3))]
temp <- cluster_0_marker_resolution_0.3[temp,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data')

#seems no very specific marker for the sub-cluster of Ex-2, consider combine them

#check cluster 6 marker
cluster_6_marker_resolution_0.3 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '6',group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)
temp <- cluster_6_marker_resolution_0.3[cluster_6_marker_resolution_0.3$pct.1 > 0.4 & cluster_6_marker_resolution_0.3$pct.2 < 0.2 & cluster_6_marker_resolution_0.3$avg_log2FC > 0.5,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data')
#specific marker detected!

#cluster 0 and 7 should consider as whole
cluster_0_7_marker_resolution_0.3 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('0','7'),group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)
temp <- cluster_0_7_marker_resolution_0.3[cluster_0_7_marker_resolution_0.3$pct.1 > 0.4 & cluster_0_7_marker_resolution_0.3$pct.2 < 0.2 & cluster_0_7_marker_resolution_0.3$avg_log2FC > 0.5,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data')

cluster_0_marker_resolution_0.3 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '0',ident.2 = '7',group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)
cluster_7_marker_resolution_0.3 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '7',ident.2 = '0',group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)

#Ex-2 marker: CHRM3, ACTN2, SNX7, LMO4, SSX2IP, NWD2, UNC5C, LRFN2, KLHL4, SNTG1

temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '0',group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)
temp <- temp[dplyr::intersect(rownames(temp),rownames(cluster_0_marker_resolution_0.3)),]
temp <- temp[temp$pct.1 > 0.3 & temp$pct.2 < 0.3 & temp$avg_log2FC > 0.25,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data')

#cluster 0 marker RASGRF2, KCTD16, ADAMTSL3, KCNH5, ZDHHC2, LPL, NEBL, VWC2L, RORB, MYO5B, LPIN2

temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '7',group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data',test.use = 'wilcox',only.pos = TRUE,verbose = TRUE)
temp <- temp[dplyr::intersect(rownames(temp),rownames(cluster_7_marker_resolution_0.3)),]
temp <- temp[temp$pct.1 > 0.4 & temp$pct.2 < 0.2 & temp$avg_log2FC > 0.5,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,group.by = 'RNA_snn_res.0.3',assay = 'RNA',slot = 'data')

#cluster 7 marker FOXP1, CNTN3, GBE1, PROM1, ARAP2, ADAMTS3, COL25A1, POU6F2, ZNF804B, IL1RAPL2, POSTN, CA10

#feature plot check
FeaturePlot(object = macaque_multiome_Seurat,
            features = c('RASGRF2', 'KCTD16', 'ADAMTSL3', 'KCNH5', 'ZDHHC2', 'LPL', 'NEBL', 'VWC2L', 'RORB', 'MYO5B', 'LPIN2'),
            pt.size = 0.1,slot = 'data')
#KCNH5 VWC2L RORB checked

FeaturePlot(object = macaque_multiome_Seurat,
            features = c('FOXP1', 'CNTN3', 'GBE1', 'PROM1', 'ARAP2', 'ADAMTS3', 'COL25A1', 'POU6F2', 'ZNF804B', 'IL1RAPL2', 'POSTN', 'CA10'),
            pt.size = 0.1,slot = 'data')
#ADAMTS3 IL1RAPL2 POSTN checked

FeaturePlot(object = macaque_multiome_Seurat,
            features = c('CHRM3', 'ACTN2', 'SNX7', 'LMO4', 'SSX2IP', 'NWD2', 'UNC5C', 'LRFN2', 'KLHL4', 'SNTG1'),
            pt.size = 0.1,slot = 'data')

#CHRM3 NWD2 UNC5C LRFN2 KLHL4 SNTG1 checked

#plot
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_2_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = c('CHRM3','NWD2','UNC5C','LRFN2','KLHL4','SNTG1'),
                   object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',
                   reduction_key = 'UMAP',plot.title = 'Ex-2 marker',guide.name = 'Expression',plot.type = 'panels',
                   scaled = 'FALSE',color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,trim = c(0,0.99),
                   print = TRUE,lims = NULL,num.panel.rows = 2)
dev.off()

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_2_sub_cluster_1_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('KCNH5','VWC2L','RORB'),
                   object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',
                   reduction_key = 'UMAP',plot.title = 'Ex-2 sub cluster 1 marker',guide.name = 'Expression',plot.type = 'panels',
                   scaled = 'FALSE',color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,trim = c(0,0.99),
                   print = TRUE,lims = NULL,num.panel.rows = 1)
dev.off()

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_2_sub_cluster_2_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('ADAMTS3','IL1RAPL2','POSTN'),
                   object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',
                   reduction_key = 'UMAP',plot.title = 'Ex-2 sub cluster 2 marker',guide.name = 'Expression',plot.type = 'panels',
                   scaled = 'FALSE',color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,trim = c(0,0.99),
                   print = TRUE,lims = NULL,num.panel.rows = 1)
dev.off()

## check IP and Ex-1 -------------------------------------------------------
p1 <- DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.1.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('EOMES'),pt.size = 0.1,slot = 'data') + theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('PPP1R17'),pt.size = 0.1,slot = 'data') + theme(aspect.ratio = 1)

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_IP_marker_featureplot.pdf',width = 24,height = 7)
p2+p3+p1+plot_layout(ncol = 3)
dev.off()

p2 <- VlnPlot(object = macaque_multiome_Seurat,features = c('EOMES'),pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.1.7',slot = 'data') + theme(aspect.ratio = 1)
p3 <- VlnPlot(object = macaque_multiome_Seurat,features = c('PPP1R17'),pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.1.7',slot = 'data') + theme(aspect.ratio = 1)
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_IP_marker_vlnplot.pdf',width = 24,height = 7)
p2+p3+p1+plot_layout(ncol = 3)
dev.off()

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_IP_marker_dotplot.pdf',width = 5,height = 8)
DotPlot(object = macaque_multiome_Seurat,assay = 'RNA',features = c('EOMES','PPP1R17'),cols = c('lightgrey','blue'),
        col.min = -2.5,col.max = 2.5,group.by = 'RNA_snn_res.1.7',scale = FALSE) + 
  theme(aspect.ratio = 4) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank())
dev.off()

#seems cluster 5 is not IP, but this can due to the sparsity of our snRNA_seq

#do label transfer
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']

#recreate macaque multiome data
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech",'donor',"cell_type")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

#add meta data
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_Seurat_PCA_clustering_meta_data.rds')
meta_data <- meta_data[rownames(macaque_multiome_Seurat@meta.data),c(8:17)]
macaque_multiome_Seurat@meta.data <- cbind(macaque_multiome_Seurat@meta.data,meta_data)

#NEURON data
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')

temp <- macaque_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 200*(1024^3),workers = 6)
macaque_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
PDhuman_RNA_seurat[['converted']] <- CreateAssayObject(counts = PDhuman_RNA_seurat@assays$RNA@counts,min.cells = 0,min.features = 0)
macaque_multiome_Seurat$species <- 'macaque'
PDhuman_RNA_seurat$species <- 'human'
DefaultAssay(PDhuman_RNA_seurat) <- 'converted'

macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(PDhuman_RNA_seurat@assays$converted@counts))
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque = macaque_multiome_Seurat,human = PDhuman_RNA_seurat),
                                           assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque = c('nCount_converted','donor'),human = c('nCount_converted')),
                                           npcs = 50,reference_loading = 'macaque',integration_var = 'species',harmony_input_dim = 30,
                                           max.iter.harmony = 50,UMAP_dim = 30,resolution = 1)

Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[rownames(PDhuman_RNA_seurat@meta.data),"cell_type"] <- PDhuman_RNA_seurat$Cluster
Brain_RNA_Seurat@meta.data[rownames(macaque_multiome_Seurat@meta.data),"cell_type"] <- macaque_multiome_Seurat$cell_type

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_PDhuman_RNA_Seurat_harmony_integration.pdf',width = 14,height = 7)
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
dev.off()

macaque_multiome_Seurat@reductions$PCA <- CreateDimReducObject(embeddings = Brain_RNA_Seurat@reductions$pca@cell.embeddings[rownames(macaque_multiome_Seurat@meta.data),],assay = 'converted')
PDhuman_RNA_seurat@reductions$PCA <- CreateDimReducObject(embeddings = Brain_RNA_Seurat@reductions$pca@cell.embeddings[rownames(PDhuman_RNA_seurat@meta.data),],assay = 'converted')
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = PDhuman_RNA_seurat,ref_reduction = 'PCA',query_reduction = 'PCA',ref_assay = 'converted',query_assay = 'converted',l2.norm = FALSE,dims = 1:30,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$cell_type,l2.norm = FALSE,verbose = TRUE,slot = 'data',dims = 1:30)
PDhuman_RNA_seurat <- AddMetaData(object = PDhuman_RNA_seurat,metadata = predictions)

scibet::Confusion_heatmap(ori = PDhuman_RNA_seurat$Cluster,prd = PDhuman_RNA_seurat$predicted.id)
PDhuman_RNA_seurat@reductions$umap <- CreateDimReducObject(embeddings = Brain_RNA_Seurat@reductions$umap@cell.embeddings[rownames(PDhuman_RNA_seurat@meta.data),],assay = 'converted')
p1 <- DimPlot(object = PDhuman_RNA_seurat,group.by = 'Cluster',label = TRUE,repel = TRUE,reduction = 'umap') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = PDhuman_RNA_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'umap') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + theme(aspect.ratio = 1)
pdf(file = './res/step_41_fig_220626/PDhuman_RNA_Seurat_harmony_label_transfer_dimplot.pdf',width = 14,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#Greenleaf data
Greenleaf_RNA_Seurat <- readRDS(file = './data/(public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$Age != 'pcw24']
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')

temp <- macaque_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 200*(1024^3),workers = 6)
macaque_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_multiome_Seurat$species <- 'macaque'
Greenleaf_RNA_Seurat$species <- 'human'

macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(Greenleaf_RNA_Seurat@assays$converted@counts))
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque = macaque_multiome_Seurat,human = Greenleaf_RNA_Seurat),
                                           assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque = c('nCount_converted','donor'),human = c('nCount_converted','Age','Batch')),
                                           npcs = 50,reference_loading = 'macaque',integration_var = 'species',harmony_input_dim = 40,
                                           max.iter.harmony = 50,UMAP_dim = 40,resolution = 1)

Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[rownames(Greenleaf_RNA_Seurat@meta.data),"cell_type"] <- Greenleaf_RNA_Seurat$cell_type
Brain_RNA_Seurat@meta.data[rownames(macaque_multiome_Seurat@meta.data),"cell_type"] <- macaque_multiome_Seurat$cell_type

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Greenleaf_RNA_Seurat_harmony_integration.pdf',width = 14,height = 7)
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
dev.off()

## check Ex-1 sub-cluster --------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']

#recreate macaque multiome data
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech",'donor',"cell_type")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

#preprocess of macaque multiome Seurat
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,
                                             nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)

ndims <- 31
macaque_multiome_Seurat <- RunUMAP(object = macaque_multiome_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

#load meta data
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_Seurat_PCA_clustering_meta_data.rds')
meta_data <- meta_data[,8:17]
macaque_multiome_Seurat@meta.data <- cbind(macaque_multiome_Seurat@meta.data,meta_data[rownames(macaque_multiome_Seurat@meta.data),])

DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.1.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

#find marker
cluster_5_6_marker_resolution_1.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('5','6'),group.by = 'RNA_snn_res.1.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_5_marker_resolution_1.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '5',ident.2 = '6',group.by = 'RNA_snn_res.1.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_6_marker_resolution_1.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '6',ident.2 = '5',group.by = 'RNA_snn_res.1.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)

#cluster 5 marker
temp <- dplyr::intersect(x = rownames(cluster_5_6_marker_resolution_1.7),rownames(cluster_5_marker_resolution_1.7))
temp <- cluster_5_marker_resolution_1.7[temp,]
VlnPlot(object = macaque_multiome_Seurat[,!(macaque_multiome_Seurat$RNA_snn_res.1.7 %in% c('25','28','31','32'))],features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.1.7',slot = 'data')

#HS3ST1 DEPTOR SLC17A6 DPY19L1 NRP1 ZMAT4 RAB12
FeaturePlot(object = macaque_multiome_Seurat,features = c('HS3ST1','DEPTOR','SLC17A6','DPY19L1','NRP1','ZMAT4','RAB12'))

#SLC17A6 DPY19L1 NRP1 checked
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_1_sub_cluster_1_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('SLC17A6','DPY19L1','NRP1'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'Ex-1 sub cluster 1 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

#cluster 6 marker
temp <- dplyr::intersect(x = rownames(cluster_5_6_marker_resolution_1.7),y = rownames(cluster_6_marker_resolution_1.7))
temp <- cluster_6_marker_resolution_1.7[temp,]
VlnPlot(object = macaque_multiome_Seurat[,!(macaque_multiome_Seurat$RNA_snn_res.1.7 %in% c('25','28','31','32'))],features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.1.7',slot = 'data')
#KCNK2 EIF1B ARHGAP18 HS6ST2 HMCN1 MCUR1 RNF182 LRRN1 GALNT10 SNTG2 GREM2 KCNK2 PRSS12 RNF182 CCBE1

FeaturePlot(object = macaque_multiome_Seurat,features = c('KCNK2','EIF1B','ARHGAP18','HS6ST2','HMCN1','MCUR1','RNF182','LRRN1','GALNT10','SNTG2','GREM2','KCNK2','PRSS12','RNF182','CCBE1'),pt.size = 0.1)

#KCNK2 HMCN1 PRSS12 CCBE1 checked
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_1_sub_cluster_2_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('KCNK2','HMCN1','PRSS12','CCBE1'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'Ex-1 sub cluster 2 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

## check Ex-3 sub cluster --------------------------------------------------
cluster_10_15_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('10','15'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_10_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '10',ident.2 = '15',group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_15_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '15',ident.2 = '10',group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)

#Ex-3 marker
temp <- cluster_10_15_marker_resolution_0.7[cluster_10_15_marker_resolution_0.7$pct.1 > 0.4 & cluster_10_15_marker_resolution_0.7$pct.2 < 0.2 & cluster_10_15_marker_resolution_0.7$avg_log2FC > 0.5,]
VlnPlot(object = macaque_multiome_Seurat[,!(macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('16','17','18','19'))],features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')
#GRIK3 IQCF2 FRAS1 RXFP1 SEMA3E PRKAG2 PDZD2 EPB41L4A GABRA1 RYR3 DLC1 TRPM3 TLE4 LMO7 HS3ST4

FeaturePlot(object = macaque_multiome_Seurat,features = c('GRIK3','IQCF2','FRAS1','RXFP1','SEMA3E','PRKAG2','PDZD2','EPB41L4A','GABRA1','RYR3','DLC1','TRPM3','TLE4','LMO7','HS3ST4'),pt.size = 0.1)
#RXFP1 PDZD2 TRPM3 TLE4 HS3ST4 checked

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_3_marker_featureplot.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = c('RXFP1','PDZD2','TRPM3','TLE4','HS3ST4'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'Ex-3 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 2)
dev.off()

#Ex-3 sub cluster 1 marker
temp <- dplyr::intersect(x = rownames(cluster_10_15_marker_resolution_0.7),y = rownames(cluster_10_marker_resolution_0.7))
temp <- cluster_10_marker_resolution_0.7[temp,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',slot = 'data',group.by = 'RNA_snn_res.0.7')
#SYT6 SEMA3E FOXP2 CDH6 KLHL32 ADGRG6 RGS6 GLRA2 GDPD5 SSTR2 NPTX1 FAM135B KLHL3 

FeaturePlot(object = macaque_multiome_Seurat,features = c('SYT6','SEMA3E','FOXP2','CDH6','KLHL32','ADGRG6','RGS6','GLRA2','GDPD5','SSTR2','NPTX1','FAM135B','KLHL3'),pt.size = 0.1,slot = 'data')
#FOXP2 KLHL32 RGS6 GLRA2 checked

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_3_sub_cluster_1_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('FOXP2','KLHL32','RGS6','GLRA2'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'Ex-3 sub cluster 1 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

#Ex-3 sub cluster 2 marker
temp <- dplyr::intersect(x = rownames(cluster_10_15_marker_resolution_0.7),y = rownames(cluster_15_marker_resolution_0.7))
temp <- cluster_15_marker_resolution_0.7[temp,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',slot = 'data',group.by = 'RNA_snn_res.0.7')
#CLSTN2 SEPRINI1 IQCF2 FAM160A1 POU6F2 ANKRD33B CDH18 SLIT3 HS3ST5 GFRA1 NGEF TMEM178A ARHGAP44 GSG1L PTPRR ADAMTSL1 NR4A2 CDH9 RERG KCNN2 LIN28B HCRTR2 SLC4A8 ZMAT4 KCNAB1 TMEM65

FeaturePlot(object = macaque_multiome_Seurat,features = c('CLSTN2','SERPINI1','IQCF2','FAM160A1','POU6F2','ANKRD33B','CDH18','SLIT3','HS3ST5','GFRA1','NGEF','TMEM178A','ARHGAP44','GSG1L','PTPRR','ADAMTSL1','NR4A2','CDH9','RERG','KCNN2','LIN28B','HCRTR2','SLC4A8','ZMAT4','KCNAB1','TMEM65'))
#CLSTN2 SLIT3 TMEM178A ADAMTSL1

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_Ex_3_sub_cluster_2_marker_featureplot.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = c('CLSTN2','SLIT3','TMEM178A','ADAMTSL1'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'Ex-3 sub cluster 2 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

## check InCGE -------------------------------------------------------------
#according to the adult brain research, InCGE can be further divided into Lamp5/Sncg and Vip sub cluster
p1 <- geneSetAveragePlot(genes = c('LAMP5','SNCG','VIP'),object = macaque_multiome_Seurat,
                         object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         plot.title = 'macaque InCGE sub cell type classic marker',plot.type = 'panels',scaled = FALSE,
                         color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                         print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)

Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')

p2 <- VlnPlot(object = Greenleaf_RNA_Seurat,features = c('LAMP5','SNCG','VIP'),pt.size = 0,assay = 'converted',slot = 'data',group.by = 'cell_type')

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_adult_InCGE_sub_cell_type_marker_featureplot.pdf',width = 12,height = 4)
p1
dev.off()

pdf(file = './res/step_41_fig_220626/Greenleaf_RNA_Seurat_adult_InCGE_sub_cell_type_marker_vlnplot.pdf',width = 16,height = 4)
p2
dev.off()

#find markers
DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + 
  DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  plot_layout(ncol = 2)

#InCGE marker
cluster_2_4_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('2','4'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_2_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('2'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_4_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('4'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)

temp <- dplyr::intersect(x = rownames(cluster_2_4_marker_resolution_0.7),y = rownames(cluster_2_marker_resolution_0.7))
temp <- dplyr::intersect(x = temp,y = rownames(cluster_4_marker_resolution_0.7))
VlnPlot(object = macaque_multiome_Seurat,features = temp[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')

#THRB LRRC1 RIPOR2 ADARB2 SORCS3 ST8SIA5
FeaturePlot(object = macaque_multiome_Seurat,features = c('THRB','LRRC1','RIPOR2','ADARB2','SORCS3','ST8SIA5','GAD1','GAD2'),pt.size = 0.1,slot = 'data')
#ADARB2 GAD1 GAD2 checked

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_InCGE_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('ADARB2','GAD1','GAD2'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'macaque InCGE marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

#InCGE subcluster 1 marker
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('2'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_2_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('2'),ident.2 = c('4'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)

temp <- dplyr::intersect(x = rownames(cluster_2_marker_resolution_0.7),y = rownames(temp))
temp <- cluster_2_marker_resolution_0.7[temp,]
temp <- temp[temp$pct.1 > 0.2 & temp$pct.2 < 0.3,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')

#KITLG PIP5K1B RND3 BMPER TMEM123 FOXO3 LRRC1
FeaturePlot(object = macaque_multiome_Seurat,features = c('KITLG','PIP5K1B','RND3','BMPER','TMEM123','FOXO3','LRRC1'),pt.size = 0,slot = 'data')

#RND3 BMPER TMEM123 checked
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_InCGE_sub_cluster_1_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('RND3','BMPER','TMEM123'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'macaque InCGE sub cluster 1 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

p1 <- DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + 
  DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  plot_layout(ncol = 2)
p2 <- VlnPlot(object = macaque_multiome_Seurat,features = 'nFeature_RNA',pt.size = 0,group.by = 'RNA_snn_res.0.7')

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_RNA_snn_res.0.7_dimplot_nFeature_RNA_vlnplot.pdf',width = 6,height = 12)
p1+p2+plot_layout(ncol = 1)
dev.off()

pdf(file = './res/step_41_fig_220626/Greenleaf_RNA_Seurat_nFeature_originalexp_vlnplot.pdf',width = 8,height = 4)
VlnPlot(object = Greenleaf_RNA_Seurat,features = 'nFeature_originalexp',pt.size = 0,group.by = 'cell_type')
dev.off()

#InCGE sub cluster 2 marker
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('4'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
cluster_4_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('4'),ident.2 = c('2'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)

temp <- dplyr::intersect(x = rownames(cluster_4_marker_resolution_0.7),y = rownames(temp))
temp <- cluster_4_marker_resolution_0.7[temp,]
temp <- temp[temp$pct.1 > 0.4 & temp$pct.2 < 0.2,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')

#NR3C2 GALNTL6 CNR1 SYNPR KCNC2 UBASH3B PROX1
FeaturePlot(object = macaque_multiome_Seurat,features = c('NR3C2','GALNTL6','CNR1','SYNPR','KCNC2','UBASH3B','PROX1'),pt.size = 0.1)

#NR3C2 CNR1 PROX1
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_InCGE_sub_cluster_2_marker_featureplot.pdf',width = 12,height = 4)
geneSetAveragePlot(genes = c('NR3C2','CNR1','PROX1'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'macaque InCGE sub cluster 2 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

## check InMGE -------------------------------------------------------------
#according to the adult brain research, InMGE can be further divided into SST and Pvalb sub cell type
p1 <- geneSetAveragePlot(genes = c('SST','PVALB'),object = macaque_multiome_Seurat,
                         object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         plot.title = 'macaque InMGE sub cell type classic marker',plot.type = 'panels',scaled = FALSE,
                         color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                         print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)

p2 <- VlnPlot(object = Greenleaf_RNA_Seurat,features = c('SST','PVALB'),pt.size = 0,assay = 'converted',slot = 'data',group.by = 'cell_type')

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_adult_InMGE_sub_cell_type_marker_featureplot.pdf',width = 8,height = 4)
p1
dev.off()

pdf(file = './res/step_41_fig_220626/Greenleaf_RNA_Seurat_adult_InMGE_sub_cell_type_marker_vlnplot.pdf',width = 10,height = 4)
p2
dev.off()

#find markers
DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + 
  DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  plot_layout(ncol = 2)

#InMGE marker
cluster_8_9_1_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('8','9','1'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
temp <- cluster_8_9_1_marker_resolution_0.7[cluster_8_9_1_marker_resolution_0.7$pct.1 > 0.3 & cluster_8_9_1_marker_resolution_0.7$pct.2 < 0.3,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')

#KCTD8 NXPH1 PAM TRERF1 DACH2 EML6 LHX6 GAD2 GAD1
FeaturePlot(object = macaque_multiome_Seurat,features = c('KCTD8','NXPH1','PAM','TRERF1','DACH2','EML6','LHX6','GAD2','GAD1'),pt.size = 0.1)

#NXPH1 EML6 LHX6 GAD1 GAD2
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_InMGE_marker_featureplot.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = c('NXPH1','EML6','LHX6','GAD1','GAD2'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'macaque InMGE marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 2)
dev.off()

#InMGE sub cluster 1 marker
cluster_9_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('9'),ident.2 = c('8','1'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('9'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
temp <- dplyr::intersect(x = rownames(temp),y = rownames(cluster_9_marker_resolution_0.7))
temp <- cluster_9_marker_resolution_0.7[temp,]
temp <- temp[temp$pct.1 > 0.4 & temp$pct.2 < 0.2 & temp$avg_log2FC > 0.5,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')

#CACNA2D3 SYNPR KLHL5 GRIK1 KCNIP1 SAMD5 PTCHD4 TRHDE GRIN3A NETO1
FeaturePlot(object = macaque_multiome_Seurat,features = c('CACNA2D3','SYNPR','KLHL5','GRIK1','KCNIP1','SAMD5','PTCHD4','TRHDE','GRIN3A','NETO1'),pt.size = 0.1)

#PTCHD4 KLHL5 GRIN3A checked
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_InMGE_sub_cluster_1_marker_featureplot.pdf',width = 8,height = 4)
geneSetAveragePlot(genes = c('PTCHD4','KLHL5'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'macaque InMGE sub cluster 1 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

#InMGE sub cluster 2 marker
cluster_8_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('8'),ident.2 = c('9','1'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('8'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
temp <- dplyr::intersect(x = rownames(temp),y = rownames(cluster_8_marker_resolution_0.7))
temp <- cluster_8_marker_resolution_0.7[temp,]
temp <- temp[temp$pct.1 > 0.3 & temp$pct.2 < 0.3,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')

#BRINP2 STK32B NXPH2 LPAR1 SVEP1
FeaturePlot(object = macaque_multiome_Seurat,features = c('BRINP2','STK32B','NXPH2','LPAR1','SVEP1'),pt.size = 0.1)
#NXPH2 LPAR1 checked

pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_InMGE_sub_cluster_2_marker_featureplot.pdf',width = 8,height = 4)
geneSetAveragePlot(genes = c('NXPH2','LPAR1'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'macaque InMGE sub cluster 2 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

#InMGE sub cluster 3 marker
cluster_1_marker_resolution_0.7 <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('1'),ident.2 = c('9','8'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('1'),group.by = 'RNA_snn_res.0.7',assay = 'RNA',slot = 'data',test.use = 'wilcox',verbose = TRUE,only.pos = TRUE)
temp <- dplyr::intersect(x = rownames(temp),y = rownames(cluster_1_marker_resolution_0.7))
temp <- cluster_1_marker_resolution_0.7[temp,]
temp <- temp[temp$pct.1 > 0.2 & temp$pct.2 < 0.3,]
VlnPlot(object = macaque_multiome_Seurat,features = rownames(temp)[1:12],pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.0.7',slot = 'data')

#DACH2 LGI1
FeaturePlot(object = macaque_multiome_Seurat,features = c('DACH2','LGI1'),pt.size = 0.1)

#DACH2 LGI1 checked
pdf(file = './res/step_41_fig_220626/macaque_multiome_Seurat_InMGE_sub_cluster_3_marker_featureplot.pdf',width = 8,height = 4)
geneSetAveragePlot(genes = c('DACH2','LGI1'),object = macaque_multiome_Seurat,
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.title = 'macaque InMGE sub cluster 3 marker',plot.type = 'panels',scaled = FALSE,
                   color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,point.size = 0.3,
                   print = TRUE,trim = c(0,0.99),lims = NULL,num.panel.rows = 1)
dev.off()

## re-annotation --------------------------------------------------------

#choose cluster resolution 0.7
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']

meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech",'donor',"cell_type")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,
                                             nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)

ndims <- 31
macaque_multiome_Seurat <- FindNeighbors(object = macaque_multiome_Seurat,reduction = 'pca',dims = 1:ndims)
macaque_multiome_Seurat <- FindClusters(object = macaque_multiome_Seurat,resolution = c(0.7))
macaque_multiome_Seurat <- RunUMAP(object = macaque_multiome_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

p1+p2+plot_layout(ncol = 2)

#check End and Per
#End: CLDN5,PECAM1
FeaturePlot(object = macaque_multiome_Seurat,features = c('CLDN5','PECAM1'),pt.size = 0.1)
#cluster 16 End

#Per: FOXC2,PDGFRB
FeaturePlot(object = macaque_multiome_Seurat,features = c('FOXC2','PDGFRB'),pt.size = 0.1)
#cliuster 18 Per

#VLMC: COL1A1,LUM
FeaturePlot(object = macaque_multiome_Seurat,features = c('COL1A1','LUM'),pt.size = 0.1)
#cluster 19 VLMC

#sub-cluster for cluster 7
Idents(macaque_multiome_Seurat) <- 'RNA_snn_res.0.7'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = '7',graph.name = 'RNA_snn',resolution = 0.15)
DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'sub.cluster',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

#annotate cell_type
macaque_multiome_Seurat$cell_type <- NA
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('13'),"cell_type"] <- 'OPC'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('12'),"cell_type"] <- 'RG-2'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('14'),"cell_type"] <- 'RG-1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('11'),"cell_type"] <- 'Cycling'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('7_1'),"cell_type"] <- 'IP'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('7_0'),"cell_type"] <- 'Ex-1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('5'),"cell_type"] <- 'Ex-2'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('0','3'),"cell_type"] <- 'Ex-3'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('6'),"cell_type"] <- 'Ex-4'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('10','15'),"cell_type"] <- 'Ex-5'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('2','4'),"cell_type"] <- 'InCGE'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('1','8','9'),"cell_type"] <- 'InMGE'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('17'),"cell_type"] <- 'Mic'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('16'),"cell_type"] <- 'End'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('18'),"cell_type"] <- 'Per'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('19'),"cell_type"] <- 'VLMC'

DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

#annotate sub_cell_type
macaque_multiome_Seurat$sub_cell_type <- macaque_multiome_Seurat$cell_type
#Ex-5
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'Ex-5',graph.name = 'RNA_snn',resolution = 0.1)
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('Ex-5_0'),"sub_cell_type"] <- 'Ex-5_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('Ex-5_1'),"sub_cell_type"] <- 'Ex-5_2'

#InCGE
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'InCGE',graph.name = 'RNA_snn',resolution = 0.1)
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('InCGE_0'),"sub_cell_type"] <- 'InCGE_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('InCGE_1'),"sub_cell_type"] <- 'InCGE_2'

#InMGE
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'InMGE',graph.name = 'RNA_snn',resolution = 0.2)
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('InMGE_0'),"sub_cell_type"] <- 'InMGE_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('InMGE_1'),"sub_cell_type"] <- 'InMGE_2'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('InMGE_2'),"sub_cell_type"] <- 'InMGE_3'

#Cycling
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'Cycling',graph.name = 'RNA_snn',resolution = 0.1)
FeaturePlot(object = macaque_multiome_Seurat,features = c('TOP2A','MKI67','CLSPN','AURKA'),pt.size = 0.1)
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('Cycling_2'),"cell_type"] <- 'IP'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('Cycling_2'),"sub_cell_type"] <- 'IP'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('Cycling_0'),"sub_cell_type"] <- 'Cyc-G2M'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('Cycling_1'),"sub_cell_type"] <- 'Cyc-S'

#RG-2
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'RG-2',graph.name = 'RNA_snn',resolution = 0.08)
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('RG-2_0'),"sub_cell_type"] <- 'RG-2_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('RG-2_1'),"sub_cell_type"] <- 'RG-2_2'

#plot
p1 <- DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + 
  theme(aspect.ratio = 1)

p2 <- DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + 
  theme(aspect.ratio = 1)

p1+p2+plot_layout(ncol = 2)

macaque_multiome_Seurat$cell_id <- rownames(macaque_multiome_Seurat@meta.data)
meta_data <- macaque_multiome_Seurat@meta.data[,c("cell_id","cell_type","sub_cell_type")]
saveRDS(meta_data,file = './res/step_41_fig_220626/macaque_multiome_Seurat_reannotate_round_1_meta_data.rds')

# validate annotation on macaque_multiome_ArchR -------------------------------------------
macaque_multiome_ArchR <- loadArchRProject(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/')
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_Seurat_reannotate_round_1_meta_data.rds')

macaque_multiome_ArchR$new_cell_type <- NA
macaque_multiome_ArchR$sub_cell_type <- NA
macaque_multiome_ArchR@cellColData[rownames(meta_data),"new_cell_type"] <- meta_data$cell_type
macaque_multiome_ArchR@cellColData[rownames(meta_data),"sub_cell_type"] <- meta_data$sub_cell_type

p1 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),
                 group.by = 'Clusters',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),
                 group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

p3 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),
                 group.by = 'new_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

p4 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),
                 group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

p1+p2+p3+p4+plot_layout(ncol = 4)

#all cell type checked on scATAC

# validate annotation on macaque scRNA ------------------------------------
#process of macaque multiome data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']

meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech",'donor',"cell_type")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,
                                             nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)

ndims <- 31
macaque_multiome_Seurat <- FindNeighbors(object = macaque_multiome_Seurat,reduction = 'pca',dims = 1:ndims)
macaque_multiome_Seurat <- FindClusters(object = macaque_multiome_Seurat,resolution = c(0.7))
macaque_multiome_Seurat <- RunUMAP(object = macaque_multiome_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

p1+p2+plot_layout(ncol = 2)

#load meta data
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_Seurat_reannotate_round_1_meta_data.rds')
macaque_multiome_Seurat$new_cell_type <- meta_data[rownames(macaque_multiome_Seurat@meta.data),"cell_type"]
macaque_multiome_Seurat$sub_cell_type <- meta_data[rownames(macaque_multiome_Seurat@meta.data),"sub_cell_type"]

p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'new_cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

p1+p2+plot_layout(ncol = 2)

#load macaque_RNA_Seurat
macaque_RNA_Seurat <- readRDS(file = './processed_data/220504_summary/macaque_RNA_Seurat_220507.rds')
meta_data <- macaque_RNA_Seurat@meta.data
meta_data <- meta_data[,c("batch","donor","tech","cell_type")]
macaque_RNA_Seurat <- macaque_RNA_Seurat@assays$RNA@counts
macaque_RNA_Seurat <- CreateSeuratObject(counts = macaque_RNA_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = VariableFeatures(macaque_multiome_Seurat),vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)

#PCA projection
temp <- projectMatrix_SeuratUMAP(X_scaled = macaque_RNA_Seurat@assays$RNA@scale.data,object = macaque_multiome_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
macaque_RNA_Seurat@reductions$projected_PCA <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = macaque_RNA_Seurat,ref_reduction = 'pca',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'RNA',dims = 1:ndims,l2.norm = FALSE,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$new_cell_type,l2.norm = FALSE,dims = 1:ndims,verbose = TRUE,slot = 'data')
macaque_RNA_Seurat <- AddMetaData(object = macaque_RNA_Seurat,metadata = predictions)
macaque_RNA_Seurat$new_cell_type <- macaque_RNA_Seurat$predicted.id

predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$sub_cell_type,l2.norm = FALSE,dims = 1:ndims,verbose = TRUE,slot = 'data')
macaque_RNA_Seurat <- AddMetaData(object = macaque_RNA_Seurat,metadata = predictions)
macaque_RNA_Seurat$sub_cell_type <- macaque_RNA_Seurat$predicted.id

#macaque RNA Seurat UMAP
macaque_RNA_Seurat <- FindNeighbors(object = macaque_RNA_Seurat,reduction = 'pca',dims = 1:ndims,assay = 'RNA')
macaque_RNA_Seurat <- FindClusters(object = macaque_RNA_Seurat,resolution = c(0.3,0.5,0.6,0.7,0.9,1.1,1.3,1.5,1.7,2.1))
macaque_RNA_Seurat <- RunUMAP(object = macaque_RNA_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)

#sub-cluster
Idents(macaque_RNA_Seurat) <- 'RNA_snn_res.0.3'
macaque_RNA_Seurat <- FindSubCluster(object = macaque_RNA_Seurat,cluster = '0',graph.name = 'RNA_snn',resolution = 0.1)
macaque_RNA_Seurat$Cluster <- macaque_RNA_Seurat$sub.cluster
Idents(macaque_RNA_Seurat) <- 'Cluster'
macaque_RNA_Seurat <- FindSubCluster(object = macaque_RNA_Seurat,cluster = '2',graph.name = 'RNA_snn',resolution = 0.1)
macaque_RNA_Seurat$Cluster <- macaque_RNA_Seurat$sub.cluster
Idents(macaque_RNA_Seurat) <- 'Cluster'
macaque_RNA_Seurat <- FindSubCluster(object = macaque_RNA_Seurat,cluster = '5',graph.name = 'RNA_snn',resolution = 0.06)
macaque_RNA_Seurat$Cluster <- macaque_RNA_Seurat$sub.cluster

p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'Cluster',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + theme(aspect.ratio = 1,legend.position = 'none')
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p3 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'new_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p4 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Cluster_coordinance_of_multiome_cell_type.pdf',width = 24,height = 6)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear))) + theme(aspect.ratio = 1,legend.position = 'none')
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p3 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'new_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p4 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Cluster_coordinance_of_multiome_cell_type_RNA_snn_res.0.7.pdf',width = 24,height = 6)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

#Ex-1 Ex-2 Ex-3 can be further divided

#marker situation
#IP
p1 <- geneSetAveragePlot(genes = c('EOMES','PPP1R17'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat IP marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_IP_marker_featureplot.pdf',width = 18,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(2,1))
dev.off()

#Ex-1
p1 <- geneSetAveragePlot(genes = c('SLC17A6','DPY19L1','NRP1'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat Ex-1 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Ex_1_marker_featureplot.pdf',width = 24,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(3,1))
dev.off()

#Ex-2
p1 <- geneSetAveragePlot(genes = c('KCNK2','HMCN1','PRSS12','CCBE1'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat Ex-2 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 2)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Ex_2_marker_featureplot.pdf',width = 18,height = 12)
p1+p2+plot_layout(ncol = 2,widths = c(2,1),heights = c(1,1))
dev.off()

#Ex-3
p1 <- geneSetAveragePlot(genes = c('KCNH5','VWC2L','RORB'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat Ex-3 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Ex_3_marker_featureplot.pdf',width = 24,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(3,1))
dev.off()

#Ex-4
p1 <- geneSetAveragePlot(genes = c('ADAMTS3','IL1RAPL2','POSTN'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat Ex-4 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Ex_4_marker_featureplot.pdf',width = 24,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(3,1))
dev.off()

#Ex-5_1
p1 <- geneSetAveragePlot(genes = c('FOXP2','KLHL32','RGS6','GLRA2'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat Ex-5_1 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 2)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Ex_5_1_marker_featureplot.pdf',width = 18,height = 12)
p1+p2+plot_layout(ncol = 2,widths = c(2,1),heights = c(1,1))
dev.off()

#Ex-5_2
p1 <- geneSetAveragePlot(genes = c('CLSTN2','SLIT3','TMEM178A','ADAMTSL1'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat Ex-5_2 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 2)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_Ex_5_2_marker_featureplot.pdf',width = 18,height = 12)
p1+p2+plot_layout(ncol = 2,widths = c(2,1),heights = c(1,1))
dev.off()

#InCGE_1
p1 <- geneSetAveragePlot(genes = c('RND3','BMPER','TMEM123'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat InCGE_1 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_InCGE_1_marker_featureplot.pdf',width = 24,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(3,1))
dev.off()

#InCGE_2
p1 <- geneSetAveragePlot(genes = c('NR3C2','CNR1','PROX1'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat InCGE_2 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_InCGE_2_marker_featureplot.pdf',width = 24,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(3,1))
dev.off()

#InMGE_3
p1 <- geneSetAveragePlot(genes = c('PTCHD4','KLHL5'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat InMGE_3 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_InMGE_3_marker_featureplot.pdf',width = 18,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(2,1))
dev.off()

#InMGE_2
p1 <- geneSetAveragePlot(genes = c('NXPH2','LPAR1'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat InMGE_2 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_InMGE_2_marker_featureplot.pdf',width = 18,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(2,1))
dev.off()

#InMGE_1
p1 <- geneSetAveragePlot(genes = c('DACH2','LGI1'),object = macaque_RNA_Seurat,object.class = 'seurat',assay = 'RNA',
                         embedding = 'umap',reduction_key = 'UMAP',plot.title = 'macaque RNA Seurat InMGE_1 marker',guide.name = 'expression',
                         plot.type = 'panels',scaled = FALSE,color.palette = ArchRPalettes$whitePurple[2:9],aspectratio = 1,
                         point.size = 0.1,trim = c(0,0.99),lims = NULL,print = TRUE,num.panel.rows = 1)

p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

pdf(file = './res/step_41_fig_220626/macaque_RNA_Seurat_InMGE_1_marker_featureplot.pdf',width = 18,height = 6)
p1+p2+plot_layout(ncol = 2,widths = c(2,1))
dev.off()

# Greenleaf data projection validation ------------------------------------
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
macaque_to_human_anno <- macaque_to_human_anno[,c(2,1,4,3)]

temp <- macaque_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 200*(1024)^3,workers = 6)
macaque_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(Greenleaf_RNA_Seurat@assays$converted@counts))
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_converted','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'cell_type',label = TRUE)
DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE)

Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,!(Greenleaf_RNA_Seurat$Age %in% c('pcw24'))]
Greenleaf_RNA_Seurat <- my_process_seurat(object = Greenleaf_RNA_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_converted','Age','Batch'),npcs = 50,preprocess = TRUE)

#projection
temp <- projectMatrix_SeuratUMAP(X_scaled = Greenleaf_RNA_Seurat@assays$converted@scale.data,object = macaque_multiome_Seurat,assayUsed = 'converted',missing_gene = FALSE)
Greenleaf_RNA_Seurat@reductions$projected_PCA <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
Greenleaf_RNA_Seurat@reductions$projected_UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

p1 <- DimPlot(object = macaque_multiome_Seurat,pt.size = 0.1,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,pt.size = 0.1,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#label tansfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = Greenleaf_RNA_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'converted',query_assay = 'converted',l2.norm = FALSE,dims = 1:30,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$new_cell_type,l2.norm = FALSE,dims = 1:30,verbose = TRUE,slot = 'data')
Greenleaf_RNA_Seurat <- AddMetaData(object = Greenleaf_RNA_Seurat,metadata = predictions)
Greenleaf_RNA_Seurat$new_cell_type <- Greenleaf_RNA_Seurat$predicted.id

p1 <- DimPlot(object = Greenleaf_RNA_Seurat,pt.size = 0.1,reduction = 'projected_UMAP',group.by = 'new_cell_type',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,pt.size = 0.1,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#the dead end