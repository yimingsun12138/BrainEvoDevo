#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Greenleaf RNA Seurat subcluster analysis                        ##
## Data: 2022.04.04                                                                ##
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
library(parallel)
library(ArchR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ComplexHeatmap)
library(dendextend)
library(Matrix)
library(matrixStats)
library(patchwork)
library(cowplot)
library(harmony)
library(scibet)
library(networkD3)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# macaque integration label transfer --------------------------------------
#load data
macaque_integration_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']
meta_data <- readRDS(file = './res/step_22_fig_220330/macaque_multiome_Seurat_meta_data.rds')
macaque_multiome_Seurat$sub_cell_type <- meta_data[colnames(macaque_multiome_Seurat),"sub_cell_type"]
DimPlot(macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_RNA_Seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'scRNA']

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = macaque_RNA_Seurat,
                                  ref_reduction = 'harmonyPCA',query_reduction = 'harmonyPCA',
                                  ref_assay = 'RNA',query_assay = 'RNA',l2.norm = FALSE,dims = 1:28,
                                  features = NULL,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$sub_cell_type,
                            weight.reduction = 'pcaproject',
                            l2.norm = FALSE,dims = 1:28,verbose = TRUE,slot = 'data')
macaque_RNA_Seurat <- AddMetaData(object = macaque_RNA_Seurat,metadata = predictions)
DimPlot(macaque_RNA_Seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

scibet::Confusion_heatmap(ori = macaque_RNA_Seurat$cell_type,prd = macaque_RNA_Seurat$predicted.id)

confusion_matrix <- my_confusion_matrix(ori = paste('ori',macaque_RNA_Seurat$cell_type,sep = '_'),prd = paste('prd',macaque_RNA_Seurat$predicted.id,sep = '_'))
nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

#add meta data
macaque_integration_Seurat$sub_cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_Seurat),"sub_cell_type"] <- macaque_multiome_Seurat$sub_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_Seurat),"sub_cell_type"] <- macaque_RNA_Seurat$predicted.id

DimPlot(macaque_integration_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

#save data
saveRDS(object = macaque_integration_Seurat@meta.data,file = './res/step_23_fig_220404/macaque_integration_Seurat_meta_data_220404.rds')

# macaque multiome Seurat integration with greenleaf RNA Seurat -----------------------------------
#load data
#human
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
Brain_integration_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_greenleaf_RNA_harmony_label_transfer_seurat_220309.rds')
DimPlot(Brain_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
#macaque
macaque_integration_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
meta_data <- readRDS(file = './res/step_23_fig_220404/macaque_integration_Seurat_meta_data_220404.rds')
macaque_integration_Seurat$sub_cell_type <- meta_data[colnames(macaque_integration_Seurat),"sub_cell_type"]

#test sub cell type
temp <- Brain_integration_Seurat[,Brain_integration_Seurat$dataset == 'macaque']
table(colnames(temp) %in% colnames(macaque_integration_Seurat))
temp$sub_cell_type <- macaque_integration_Seurat@meta.data[colnames(temp),"sub_cell_type"]
DimPlot(temp,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = c(as.character(ArchRPalettes$stallion),'#7DD06F','#844081'))

#label transfer
temp_macaque <- Brain_integration_Seurat[,Brain_integration_Seurat$dataset == 'macaque']
temp_human <- Brain_integration_Seurat[,Brain_integration_Seurat$dataset == 'human']
temp_macaque$sub_cell_type <- macaque_integration_Seurat@meta.data[colnames(temp_macaque),"sub_cell_type"]

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:30,features = NULL,
                                  verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$sub_cell_type,weight.reduction = 'pcaproject',
                            l2.norm = FALSE,dims = 1:30,verbose = TRUE,slot = 'data')
temp_human <- AddMetaData(object = temp_human,metadata = predictions)
DimPlot(temp_human,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'umap') + 
  scale_color_manual(values = c(as.character(ArchRPalettes$stallion),'#7DD06F','#844081'))

#save data
saveRDS(temp_human@meta.data,file = './res/step_23_fig_220404/Greenleaf_RNA_Seurat_label_transfer_meta_data_220405.rds')

# cell type level marker plot -------------------------------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
meta_data <- readRDS(file = './res/step_23_fig_220404/Greenleaf_RNA_Seurat_label_transfer_meta_data_220405.rds')
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,rownames(meta_data)]
Greenleaf_RNA_Seurat$sub_cell_type <- meta_data$predicted.id
Greenleaf_RNA_Seurat$cell_type <- meta_data$cell_type
DimPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP')

#process data
Greenleaf_RNA_Seurat <- my_process_seurat(object = Greenleaf_RNA_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#RG marker
RG_marker <- read.csv(file = './res/step_22_fig_220401/RG_marker.csv')
RG_marker <- RG_marker[,1]
table(RG_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
RG_marker <- RG_marker[RG_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/RG_marker_vlnplot.pdf',width = 24,height = (round(length(RG_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = RG_marker,pt.size = 0,ncol = 3)
dev.off()

#RG_1 marker
RG_1_marker <- read.csv(file = './res/step_22_fig_220401/RG_1_marker.csv')
RG_1_marker <- RG_1_marker[,1]
table(RG_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
RG_1_marker <- RG_1_marker[RG_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/RG_1_marker_vlnplot.pdf',width = 24,height = (round(length(RG_1_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = RG_1_marker,pt.size = 0,ncol = 3)
dev.off()

#RG_2 marker
RG_2_marker <- read.csv(file = './res/step_22_fig_220401/RG_2_marker.csv')
RG_2_marker <- RG_2_marker[,1]
table(RG_2_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
RG_2_marker <- RG_2_marker[RG_2_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/RG_2_marker_vlnplot.pdf',width = 24,height = (round(length(RG_2_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = RG_2_marker,pt.size = 0,ncol = 3)
dev.off()

#RG_3 marker
RG_3_marker <- read.csv(file = './res/step_22_fig_220401/RG_3_marker.csv')
RG_3_marker <- RG_3_marker[,1]
table(RG_3_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
RG_3_marker <- RG_3_marker[RG_3_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/RG_3_marker_vlnplot.pdf',width = 24,height = (round(length(RG_3_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = RG_3_marker,pt.size = 0,ncol = 3)
dev.off()

#IP marker
IP_marker <- read.csv(file = './res/step_22_fig_220401/IP_marker.csv')
IP_marker <- IP_marker[,1]
table(IP_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
IP_marker <- IP_marker[IP_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/IP_marker_vlnplot.pdf',width = 24,height = (round(length(IP_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = IP_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_1 marker
Ex_1_marker <- read.csv(file = './res/step_22_fig_220401/Ex_1_marker.csv')
Ex_1_marker <- Ex_1_marker[,1]
table(Ex_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_1_marker <- Ex_1_marker[Ex_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_1_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_1_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = Ex_1_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_2 marker
Ex_2_marker <- read.csv(file = './res/step_22_fig_220401/Ex_2_marker.csv')
Ex_2_marker <- Ex_2_marker[,1]
table(Ex_2_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_2_marker <- Ex_2_marker[Ex_2_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_2_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_2_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = Ex_2_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_3 marker
Ex_3_marker <- read.csv(file = './res/step_22_fig_220401/Ex_3_marker.csv')
Ex_3_marker <- Ex_3_marker[,1]
table(Ex_3_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_3_marker <- Ex_3_marker[Ex_3_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_3_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_3_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = Ex_3_marker,pt.size = 0,ncol = 3)
dev.off()

#InCGE marker
InCGE_marker <- read.csv(file = './res/step_22_fig_220401/InCGE_marker.csv')
InCGE_marker <- InCGE_marker[,1]
table(InCGE_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
InCGE_marker <- InCGE_marker[InCGE_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/InCGE_marker_vlnplot.pdf',width = 24,height = (round(length(InCGE_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = InCGE_marker,pt.size = 0,ncol = 3)
dev.off()

#InMGE marker
InMGE_marker <- read.csv(file = './res/step_22_fig_220401/InMGE_marker.csv')
InMGE_marker <- InMGE_marker[,1]
table(InMGE_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
InMGE_marker <- InMGE_marker[InMGE_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/InMGE_marker_vlnplot.pdf',width = 24,height = (round(length(InMGE_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',features = InMGE_marker,pt.size = 0,ncol = 3)
dev.off()

# sub cell type marker ----------------------------------------------------
#Ex_1_0 marker
Ex_1_0_marker <- read.csv(file = './res/step_22_fig_220401/Ex_1_0_marker.csv')
Ex_1_0_marker <- Ex_1_0_marker[,1]
table(Ex_1_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_1_0_marker <- Ex_1_0_marker[Ex_1_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_1_0_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_1_0_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = Ex_1_0_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_1_1 marker
Ex_1_1_marker <- read.csv(file = './res/step_22_fig_220401/Ex_1_1_marker.csv')
Ex_1_1_marker <- Ex_1_1_marker[,1]
table(Ex_1_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_1_1_marker <- Ex_1_1_marker[Ex_1_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_1_1_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_1_1_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = Ex_1_1_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_2_0 marker
Ex_2_0_marker <- read.csv(file = './res/step_22_fig_220401/Ex_2_0_marker.csv')
Ex_2_0_marker <- Ex_2_0_marker[,1]
table(Ex_2_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_2_0_marker <- Ex_2_0_marker[Ex_2_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_2_0_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_2_0_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = Ex_2_0_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_2_1 marker
Ex_2_1_marker <- read.csv(file = './res/step_22_fig_220401/Ex_2_1_marker.csv')
Ex_2_1_marker <- Ex_2_1_marker[,1]
table(Ex_2_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_2_1_marker <- Ex_2_1_marker[Ex_2_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_2_1_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_2_1_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = Ex_2_1_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_3_0 marker
Ex_3_0_marker <- read.csv(file = './res/step_22_fig_220401/Ex_3_0_marker.csv')
Ex_3_0_marker <- Ex_3_0_marker[,1]
table(Ex_3_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_3_0_marker <- Ex_3_0_marker[Ex_3_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_3_0_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_3_0_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = Ex_3_0_marker,pt.size = 0,ncol = 3)
dev.off()

#Ex_3_1 marker
Ex_3_1_marker <- read.csv(file = './res/step_22_fig_220401/Ex_3_1_marker.csv')
Ex_3_1_marker <- Ex_3_1_marker[,1]
table(Ex_3_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
Ex_3_1_marker <- Ex_3_1_marker[Ex_3_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/Ex_3_1_marker_vlnplot.pdf',width = 24,height = (round(length(Ex_3_1_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = Ex_3_1_marker,pt.size = 0,ncol = 3)
dev.off()

#InCGE_0 marker
InCGE_0_marker <- read.csv(file = './res/step_22_fig_220401/InCGE_0_marker.csv')
InCGE_0_marker <- InCGE_0_marker[,1]
table(InCGE_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
InCGE_0_marker <- InCGE_0_marker[InCGE_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/InCGE_0_marker_vlnplot.pdf',width = 24,height = (round(length(InCGE_0_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = InCGE_0_marker,pt.size = 0,ncol = 3)
dev.off()

#InCGE_1 marker
InCGE_1_marker <- read.csv(file = './res/step_22_fig_220401/InCGE_1_marker.csv')
InCGE_1_marker <- InCGE_1_marker[,1]
table(InCGE_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
InCGE_1_marker <- InCGE_1_marker[InCGE_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/InCGE_1_marker_vlnplot.pdf',width = 24,height = (round(length(InCGE_1_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = InCGE_1_marker,pt.size = 0,ncol = 3)
dev.off()

#InMGE_0 marker
InMGE_0_marker <- read.csv(file = './res/step_22_fig_220401/InMGE_0_marker.csv')
InMGE_0_marker <- InMGE_0_marker[,1]
table(InMGE_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
InMGE_0_marker <- InMGE_0_marker[InMGE_0_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/InMGE_0_marker_vlnplot.pdf',width = 24,height = (round(length(InMGE_0_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = InMGE_0_marker,pt.size = 0,ncol = 3)
dev.off()

#InMGE_1 marker
InMGE_1_marker <- read.csv(file = './res/step_22_fig_220401/InMGE_1_marker.csv')
InMGE_1_marker <- InMGE_1_marker[,1]
table(InMGE_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
InMGE_1_marker <- InMGE_1_marker[InMGE_1_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/InMGE_1_marker_vlnplot.pdf',width = 24,height = (round(length(InMGE_1_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = InMGE_1_marker,pt.size = 0,ncol = 3)
dev.off()

#InMGE_2 marker
InMGE_2_marker <- read.csv(file = './res/step_22_fig_220401/InMGE_2_marker.csv')
InMGE_2_marker <- InMGE_2_marker[,1]
table(InMGE_2_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data))
InMGE_2_marker <- InMGE_2_marker[InMGE_2_marker %in% rownames(Greenleaf_RNA_Seurat@assays$converted@data)]
pdf(file = './res/step_23_fig_220404/InMGE_2_marker_vlnplot.pdf',width = 24,height = (round(length(InMGE_2_marker)/3)+1)*4)
VlnPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',features = InMGE_2_marker,pt.size = 0,ncol = 3)
dev.off()

# marker overlap ----------------------------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
meta_data <- readRDS(file = './res/step_23_fig_220404/Greenleaf_RNA_Seurat_label_transfer_meta_data_220405.rds')
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,rownames(meta_data)]
Greenleaf_RNA_Seurat$sub_cell_type <- meta_data$predicted.id
Greenleaf_RNA_Seurat$cell_type <- meta_data$cell_type
DimPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP')
Greenleaf_RNA_Seurat <- my_process_seurat(object = Greenleaf_RNA_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']
meta_data <- readRDS(file = './res/step_22_fig_220330/macaque_multiome_Seurat_meta_data.rds')
macaque_multiome_Seurat$sub_cell_type <- meta_data[colnames(macaque_multiome_Seurat),"sub_cell_type"]
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

## cell type marker --------------------------------------------------------
#macaque
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_marker <- FindAllMarkers(object = macaque_multiome_Seurat,assay = 'RNA',test.use = 'bimod',slot = 'data',verbose = TRUE,only.pos = TRUE)
temp <- macaque_marker[macaque_marker$pct.1 > 0.6 & macaque_marker$pct.2 < 0.4,]
table(temp$cluster)

#human
Idents(Greenleaf_RNA_Seurat) <- 'cell_type'
human_marker <- FindAllMarkers(object = Greenleaf_RNA_Seurat,assay = 'converted',test.use = 'bimod',slot = 'data',verbose = TRUE,only.pos = TRUE)
temp <- human_marker[human_marker$pct.1 > 0.4 & human_marker$pct.2 < 0.4,]
table(temp$cluster)

#calculate jaccard index
JD_index <- matrix(nrow = length(unique(Greenleaf_RNA_Seurat$cell_type)),ncol = length(unique(macaque_multiome_Seurat$cell_type)))
JD_index <- as.data.frame(JD_index)
rownames(JD_index) <- unique(Greenleaf_RNA_Seurat$cell_type)
colnames(JD_index) <- unique(macaque_multiome_Seurat$cell_type)

#all marker
for (i in rownames(JD_index)) {
  for (j in colnames(JD_index)) {
    temp_human <- human_marker
    temp_human <- unique(temp_human[temp_human$cluster == i,]$gene)
    temp_macaque <- macaque_marker
    temp_macaque <- unique(temp_macaque[temp_macaque$cluster == j,]$gene)
    JD_index[i,j] <- length(dplyr::intersect(temp_human,temp_macaque))/length(dplyr::union(temp_human,temp_macaque))
  }
}

JD_index <- JD_index[c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE'),c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE')]

pdf(file = './res/step_23_fig_220404/macaque_Greenleaf_cell_type_marker_all_marker_overlap_heatmap.pdf',width = 8,height = 7)
Heatmap(matrix = JD_index,cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = TRUE,show_column_names = TRUE,
        row_title = 'Human marker all',column_title = 'Macaque marker all',
        row_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        column_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        border = TRUE,height = unit(5,'inches'),width = unit(5,'inches'),name = 'Jaccard idx')
dev.off()

#pct2 < 0.6
for (i in rownames(JD_index)) {
  for (j in colnames(JD_index)) {
    temp_human <- human_marker[human_marker$pct.2 < 0.6,]
    temp_human <- unique(temp_human[temp_human$cluster == i,]$gene)
    temp_macaque <- macaque_marker[macaque_marker$pct.2 < 0.6,]
    temp_macaque <- unique(temp_macaque[temp_macaque$cluster == j,]$gene)
    JD_index[i,j] <- length(dplyr::intersect(temp_human,temp_macaque))/length(dplyr::union(temp_human,temp_macaque))
  }
}

JD_index <- JD_index[c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE'),c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE')]

pdf(file = './res/step_23_fig_220404/macaque_Greenleaf_cell_type_marker_pct2_0.6_overlap_heatmap.pdf',width = 8,height = 7)
Heatmap(matrix = JD_index,cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = TRUE,show_column_names = TRUE,
        row_title = 'Human marker pct2 0.6',column_title = 'Macaque marker pct2 0.6',
        row_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        column_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        border = TRUE,height = unit(5,'inches'),width = unit(5,'inches'),name = 'Jaccard idx')
dev.off()

#pct2 < 0.5
for (i in rownames(JD_index)) {
  for (j in colnames(JD_index)) {
    temp_human <- human_marker[human_marker$pct.2 < 0.5,]
    temp_human <- unique(temp_human[temp_human$cluster == i,]$gene)
    temp_macaque <- macaque_marker[macaque_marker$pct.2 < 0.5,]
    temp_macaque <- unique(temp_macaque[temp_macaque$cluster == j,]$gene)
    JD_index[i,j] <- length(dplyr::intersect(temp_human,temp_macaque))/length(dplyr::union(temp_human,temp_macaque))
  }
}

JD_index <- JD_index[c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE'),c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE')]

pdf(file = './res/step_23_fig_220404/macaque_Greenleaf_cell_type_marker_pct2_0.5_overlap_heatmap.pdf',width = 8,height = 7)
Heatmap(matrix = JD_index,cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = TRUE,show_column_names = TRUE,
        row_title = 'Human marker pct2 0.5',column_title = 'Macaque marker pct2 0.5',
        row_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        column_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        border = TRUE,height = unit(5,'inches'),width = unit(5,'inches'),name = 'Jaccard idx')
dev.off()

#pct2 < 0.4
for (i in rownames(JD_index)) {
  for (j in colnames(JD_index)) {
    temp_human <- human_marker[human_marker$pct.2 < 0.4,]
    temp_human <- unique(temp_human[temp_human$cluster == i,]$gene)
    temp_macaque <- macaque_marker[macaque_marker$pct.2 < 0.4,]
    temp_macaque <- unique(temp_macaque[temp_macaque$cluster == j,]$gene)
    JD_index[i,j] <- length(dplyr::intersect(temp_human,temp_macaque))/length(dplyr::union(temp_human,temp_macaque))
  }
}

JD_index <- JD_index[c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE'),c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE')]

pdf(file = './res/step_23_fig_220404/macaque_Greenleaf_cell_type_marker_pct2_0.4_overlap_heatmap.pdf',width = 8,height = 7)
Heatmap(matrix = JD_index,cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = TRUE,show_column_names = TRUE,
        row_title = 'Human marker pct2 0.4',column_title = 'Macaque marker pct2 0.4',
        row_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        column_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        border = TRUE,height = unit(5,'inches'),width = unit(5,'inches'),name = 'Jaccard idx')
dev.off()

#pct2 < 0.4 pct1 > 0.6
for (i in rownames(JD_index)) {
  for (j in colnames(JD_index)) {
    temp_human <- human_marker[human_marker$pct.2 < 0.4 & human_marker$pct.1 > 0.6,]
    temp_human <- unique(temp_human[temp_human$cluster == i,]$gene)
    temp_macaque <- macaque_marker[macaque_marker$pct.2 < 0.4 & macaque_marker$pct.1 > 0.6,]
    temp_macaque <- unique(temp_macaque[temp_macaque$cluster == j,]$gene)
    JD_index[i,j] <- length(dplyr::intersect(temp_human,temp_macaque))/length(dplyr::union(temp_human,temp_macaque))
  }
}

JD_index <- JD_index[c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE'),c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE')]

pdf(file = './res/step_23_fig_220404/macaque_Greenleaf_cell_type_marker_pct2_0.4_pct1_0.6_overlap_heatmap.pdf',width = 8,height = 7)
Heatmap(matrix = JD_index,cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = TRUE,show_column_names = TRUE,
        row_title = 'Human marker pct2 0.4 pct1 0.6',column_title = 'Macaque marker pct2 0.4 pct1 0.6',
        row_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        column_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        border = TRUE,height = unit(5,'inches'),width = unit(5,'inches'),name = 'Jaccard idx')
dev.off()

#pct2 < 0.4 pct1 > 0.5
for (i in rownames(JD_index)) {
  for (j in colnames(JD_index)) {
    temp_human <- human_marker[human_marker$pct.2 < 0.4 & human_marker$pct.1 > 0.5,]
    temp_human <- unique(temp_human[temp_human$cluster == i,]$gene)
    temp_macaque <- macaque_marker[macaque_marker$pct.2 < 0.4 & macaque_marker$pct.1 > 0.5,]
    temp_macaque <- unique(temp_macaque[temp_macaque$cluster == j,]$gene)
    JD_index[i,j] <- length(dplyr::intersect(temp_human,temp_macaque))/length(dplyr::union(temp_human,temp_macaque))
  }
}

JD_index <- JD_index[c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE'),c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE')]

pdf(file = './res/step_23_fig_220404/macaque_Greenleaf_cell_type_marker_pct2_0.4_pct1_0.5_overlap_heatmap.pdf',width = 8,height = 7)
Heatmap(matrix = JD_index,cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = TRUE,show_column_names = TRUE,
        row_title = 'Human marker pct2 0.4 pct1 0.5',column_title = 'Macaque marker pct2 0.4 pct1 0.5',
        row_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        column_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        border = TRUE,height = unit(5,'inches'),width = unit(5,'inches'),name = 'Jaccard idx')
dev.off()

#pct2 < 0.4 pct1 > 0.4
for (i in rownames(JD_index)) {
  for (j in colnames(JD_index)) {
    temp_human <- human_marker[human_marker$pct.2 < 0.4 & human_marker$pct.1 > 0.4,]
    temp_human <- unique(temp_human[temp_human$cluster == i,]$gene)
    temp_macaque <- macaque_marker[macaque_marker$pct.2 < 0.4 & macaque_marker$pct.1 > 0.4,]
    temp_macaque <- unique(temp_macaque[temp_macaque$cluster == j,]$gene)
    JD_index[i,j] <- length(dplyr::intersect(temp_human,temp_macaque))/length(dplyr::union(temp_human,temp_macaque))
  }
}

JD_index <- JD_index[c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE'),c('Ependymal','Per','End','Mic','OPC','RG-3','RG-2','RG-1','Cyc-S','Cyc-G2M','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE')]

pdf(file = './res/step_23_fig_220404/macaque_Greenleaf_cell_type_marker_pct2_0.4_pct1_0.4_overlap_heatmap.pdf',width = 8,height = 7)
Heatmap(matrix = JD_index,cluster_rows = FALSE,cluster_columns = FALSE,
        show_row_names = TRUE,show_column_names = TRUE,
        row_title = 'Human marker pct2 0.4 pct1 0.4',column_title = 'Macaque marker pct2 0.4 pct1 0.4',
        row_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        column_split = factor(c('nn','nn','nn','nn','nn','npc','npc','npc','npc','npc','Ex','Ex','Ex','Ex','In','In'),levels = c('nn','npc','Ex','In')),
        border = TRUE,height = unit(5,'inches'),width = unit(5,'inches'),name = 'Jaccard idx')
dev.off()

# Greenleaf RNA Seurat dimplot --------------------------------------------
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_greenleaf_RNA_harmony_label_transfer_seurat_220309.rds')
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$dataset == 'human']
meta_data <- readRDS(file = './res/step_23_fig_220404/Greenleaf_RNA_Seurat_label_transfer_meta_data_220405.rds')
meta_data <- meta_data[colnames(Greenleaf_RNA_Seurat),]
Greenleaf_RNA_Seurat$sub_cell_type <- meta_data$predicted.id

p1 <- DimPlot(Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(Greenleaf_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = c(as.character(ArchRPalettes$stallion),'#7DD06F','#844081')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_23_fig_220404/Greenleaf_RNA_seurat_dimplot.pdf',width = 18,height = 8)
p1+p2+plot_layout(ncol = 2)
dev.off()