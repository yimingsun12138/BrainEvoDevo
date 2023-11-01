#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: cellchat analysis on macaque multiome data                      ##
## Data: 2023.02.03                                                                ##
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
setwd('/content/data/sunym/project/Brain/')
.libPaths('/content/data/sunym/software/R_lib/R_4.2.2/')
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
library(ggpubr)
library(ggtext)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(BSgenome.Mmusculus.UCSC.mm10)
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(universalmotif)
library(topGO)
library(future.apply)
library(transPlotR)
library(aplot)
library(Biostrings)
library(riverplot)
library(CellChat)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

# create cellchat object --------------------------------------------------

#create cellchat object
express_matrix <- macaque_multiome_Seurat@assays$RNA@data
meta_data <- macaque_multiome_Seurat@meta.data
macaque_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'cell_type')

#set LR database
macaque_CellChat@DB <- CellChatDB.human

# preprocess expression data ----------------------------------------------
macaque_CellChat <- subsetData(object = macaque_CellChat)
future::plan(strategy = 'multisession',workers = 5)
future::plan(strategy = 'sequential')

#identify overexpression genes and interactions
macaque_CellChat <- identifyOverExpressedGenes(object = macaque_CellChat)
macaque_CellChat <- identifyOverExpressedInteractions(object = macaque_CellChat)

#project to human PPI network
macaque_CellChat_PPI <- projectData(object = macaque_CellChat,adjMatrix = PPI.human)

# infer cell cell communication network -----------------------------------
macaque_CellChat <- computeCommunProb(object = macaque_CellChat,raw.use = TRUE)
macaque_CellChat <- filterCommunication(object = macaque_CellChat,min.cells = 10)

macaque_CellChat_PPI <- computeCommunProb(object = macaque_CellChat_PPI,raw.use = FALSE)
macaque_CellChat_PPI <- filterCommunication(object = macaque_CellChat_PPI,min.cells = 10)

#extract inferred cell cell communication network
df.net <- subsetCommunication(object = macaque_CellChat)
df.net.PPI <- subsetCommunication(object = macaque_CellChat_PPI)

#infer cell cell communication at signaling pathway level
macaque_CellChat <- computeCommunProbPathway(object = macaque_CellChat)
macaque_CellChat_PPI <- computeCommunProbPathway(object = macaque_CellChat_PPI)

#calculate aggregated cell cell communication network
macaque_CellChat <- aggregateNet(object = macaque_CellChat)
macaque_CellChat_PPI <- aggregateNet(object = macaque_CellChat_PPI)

#visualization
groupSize <- as.numeric(table(macaque_CellChat@idents))
temp <- names(col_param$celltype)[names(col_param$celltype) %in% macaque_CellChat@idents]
macaque_CellChat@idents <- factor(macaque_CellChat@idents,levels = temp)

netVisual_circle(net = macaque_CellChat@net$count,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Number of interactions")

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_cell_type_communication_strength.pdf',width = 8,height = 8)
netVisual_circle(net = macaque_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength")
dev.off()

netVisual_circle(net = macaque_CellChat_PPI@net$count,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Number of interactions")

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_PPI_cell_type_communication_strength.pdf',width = 8,height = 8)
netVisual_circle(net = macaque_CellChat_PPI@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength")
dev.off()

# PDGFD signal ------------------------------------------------------------
'PDGFD' %in% df.net$ligand
temp <- df.net[df.net$ligand == 'PDGFD',]
temp

#no RG-1 related cell cell communication

#visualize all PDGFD mediated cell cell communication
pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = 'PDGF',geneLR.return = FALSE)

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_PDGFD_signal_network.pdf',width = 5,height = 5)
netVisual_individual(object = macaque_CellChat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle')
dev.off()

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_PPI_PDGFD_signal_network.pdf',width = 5,height = 5)
netVisual_individual(object = macaque_CellChat_PPI,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle')
dev.off()

#PPI projection seems have strong effect on cell cell communication

# cell communication related to RG-1 --------------------------------------

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_RG_1_related_cell_communication.pdf',width = 5,height = 5)
netVisual_circle(net = macaque_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as source",sources.use = 'RG-1')
netVisual_circle(net = macaque_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as target",targets.use = 'RG-1')
dev.off()

netVisual_circle(net = macaque_CellChat_PPI@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as source",sources.use = 'RG-1')
netVisual_circle(net = macaque_CellChat_PPI@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as target",targets.use = 'RG-1')

#no big difference between raw data and PPI projected data.

# overlap with RG-1 marker ------------------------------------------------
RG_1_marker_list <- readRDS(file = './res/step_97_fig_230130/RG_1_marker_list_from_liuyt.rds')

## specific in human -------------------------------------------------------
gene_list <- RG_1_marker_list$`specific-in-H`
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#visualize
pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = macaque_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(macaque_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 3/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 H marker mediated interaction') + 
  scale_x_discrete(labels = levels(macaque_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = macaque_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(macaque_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 2/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(macaque_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_RG_1_interaction_by_H_marker.pdf',width = 8,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

## conserved in primate ----------------------------------------------------
gene_list <- RG_1_marker_list$`conserved-in-HR`
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#visualize
pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = macaque_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(macaque_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 2/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 HR marker mediated interaction') + 
  scale_x_discrete(labels = levels(macaque_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = macaque_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(macaque_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 2/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(macaque_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_RG_1_interaction_by_HR_marker.pdf',width = 8,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

## conserved in all species ------------------------------------------------
gene_list <- RG_1_marker_list$`conserved-in-HRM`
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#visualize
pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = macaque_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(macaque_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 10/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 HRM marker mediated interaction') + 
  scale_x_discrete(labels = levels(macaque_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = macaque_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(macaque_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 5/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(macaque_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_98_fig_230203/macaque_multiome_cellchat_RG_1_interaction_by_HRM_marker.pdf',width = 8,height = 8)
p1+p2+plot_layout(ncol = 1)
dev.off()