#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: cellchat analysis on Greenleaf data                             ##
## Data: 2023.01.30                                                                ##
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
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

# create cellchat object --------------------------------------------------

#create cellchat object
express_matrix <- Greenleaf_RNA_Seurat@assays$RNA@data
meta_data <- Greenleaf_RNA_Seurat@meta.data
Greenleaf_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'ReAnno_celltype')

#set LR database
Greenleaf_CellChat@DB <- CellChatDB.human

# preprocess expression data ----------------------------------------------
Greenleaf_CellChat <- subsetData(object = Greenleaf_CellChat)
future::plan(strategy = 'multisession',workers = 5)
future::plan(strategy = 'sequential')

#identify overexpression genes and interactions
Greenleaf_CellChat <- identifyOverExpressedGenes(object = Greenleaf_CellChat)
Greenleaf_CellChat <- identifyOverExpressedInteractions(object = Greenleaf_CellChat)

#project to human PPI network
Greenleaf_CellChat_PPI <- projectData(object = Greenleaf_CellChat,adjMatrix = PPI.human)

# infer cell cell communication network -----------------------------------
Greenleaf_CellChat <- computeCommunProb(object = Greenleaf_CellChat,raw.use = TRUE)
Greenleaf_CellChat <- filterCommunication(object = Greenleaf_CellChat,min.cells = 10)

Greenleaf_CellChat_PPI <- computeCommunProb(object = Greenleaf_CellChat_PPI,raw.use = FALSE)
Greenleaf_CellChat_PPI <- filterCommunication(object = Greenleaf_CellChat_PPI,min.cells = 10)

#extract inferred cell cell communication network
df.net <- subsetCommunication(object = Greenleaf_CellChat)
df.net.PPI <- subsetCommunication(object = Greenleaf_CellChat_PPI)

#infer cell cell communication at signaling pathway level
Greenleaf_CellChat <- computeCommunProbPathway(object = Greenleaf_CellChat)
Greenleaf_CellChat_PPI <- computeCommunProbPathway(object = Greenleaf_CellChat_PPI)

#calculate aggregated cell cell communication network
Greenleaf_CellChat <- aggregateNet(object = Greenleaf_CellChat)
Greenleaf_CellChat_PPI <- aggregateNet(object = Greenleaf_CellChat_PPI)

#visualization
groupSize <- as.numeric(table(Greenleaf_CellChat@idents))
temp <- names(col_param$celltype)[names(col_param$celltype) %in% Greenleaf_CellChat@idents]
Greenleaf_CellChat@idents <- factor(Greenleaf_CellChat@idents,levels = temp)

netVisual_circle(net = Greenleaf_CellChat@net$count,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Number of interactions")

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_cell_type_communication_strength.pdf',width = 8,height = 8)
netVisual_circle(net = Greenleaf_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength")
dev.off()

netVisual_circle(net = Greenleaf_CellChat_PPI@net$count,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Number of interactions")

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_PPI_cell_type_communication_strength.pdf',width = 8,height = 8)
netVisual_circle(net = Greenleaf_CellChat_PPI@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength")
dev.off()

# PDGFD signal ------------------------------------------------------------
'PDGFD' %in% df.net$ligand
temp <- df.net[df.net$ligand == 'PDGFD',]
temp

#visualize all PDGFD mediated cell cell communication
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = 'PDGF',geneLR.return = FALSE)

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_PDGFD_signal_network.pdf',width = 5,height = 5)
netVisual_individual(object = Greenleaf_CellChat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle')
dev.off()

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_PPI_PDGFD_signal_network.pdf',width = 5,height = 5)
netVisual_individual(object = Greenleaf_CellChat_PPI,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle')
dev.off()

#seems PPI network change the expression of PDGFRB
VlnPlot(object = Greenleaf_RNA_Seurat,features = c('PDGFRB'),pt.size = 0,assay = 'RNA',slot = 'data',group.by = 'ReAnno_celltype')
p <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('PDGFRB'),pt.size = 0.001,cols = ArchRPalettes$whitePurple[-1],order = TRUE) + 
  theme_cowplot() + NoAxes() + 
  theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_97_fig_230130/Greenleaf_RNA_Seurat_raw_data_PDGFRB_expression.pdf',width = 4,height = 4)
p
dev.off()

temp <- Greenleaf_CellChat_PPI@data.project
temp <- CreateSeuratObject(counts = temp,project = 'temp',assay = 'RNA',meta.data = Greenleaf_RNA_Seurat@meta.data[rownames(temp),c("ReAnno_celltype","cell_type")],min.cells = 0,min.features = 0)
temp[['umap']] <- CreateDimReducObject(embeddings = Greenleaf_RNA_Seurat@reductions$umap@cell.embeddings,assay = 'RNA')
p <- FeaturePlot(object = temp,features = c('PDGFRB'),pt.size = 0.001,cols = ArchRPalettes$whitePurple[-1],order = TRUE,slot = 'counts') + 
  theme_cowplot() + NoAxes() + 
  theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(hjust = 0.5))
p
#that is strange

# cell communication related to RG-1 --------------------------------------

#global visualize
pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_1_realated_cell_communication.pdf',width = 5,height = 5)
netVisual_circle(net = Greenleaf_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as source",sources.use = 'RG-1')
netVisual_circle(net = Greenleaf_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as target",targets.use = 'RG-1')
dev.off()

#the strongest cell communication related to RG-1
p <- netVisual_bubble(object = Greenleaf_CellChat,sources.use = 'RG-1',remove.isolate = FALSE)

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_1_sourced_significant_cell_communication.pdf',width = 8,height = 16)
p + theme_bw() + 
  theme(aspect.ratio = 2.5,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 sourced significant cell communication')
dev.off()

p <- netVisual_bubble(object = Greenleaf_CellChat,targets.use = 'RG-1',remove.isolate = FALSE)

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_1_targeted_significant_cell_communication.pdf',width = 8,height = 18)
p + theme_bw() + 
  theme(aspect.ratio = 3,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 targeted significant cell communication')
dev.off()

#the strongest cell communication related to RG-2
p <- netVisual_bubble(object = Greenleaf_CellChat,sources.use = 'RG-2',remove.isolate = FALSE)

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_2_sourced_significant_cell_communication.pdf',width = 8,height = 16)
p + theme_bw() + 
  theme(aspect.ratio = 2.5,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-2 sourced significant cell communication')
dev.off()

p <- netVisual_bubble(object = Greenleaf_CellChat,targets.use = 'RG-2',remove.isolate = FALSE)

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_2_targeted_significant_cell_communication.pdf',width = 8,height = 18)
p + theme_bw() + 
  theme(aspect.ratio = 3,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-2 targeted significant cell communication')
dev.off()

# Identify outgoing signaling contribute to each cell type ----------------
# 
# #compute net centrality score
# Greenleaf_CellChat <- netAnalysis_computeCentrality(object = Greenleaf_CellChat,slot.name = 'netP')
# Greenleaf_CellChat_PPI <- netAnalysis_computeCentrality(object = Greenleaf_CellChat_PPI,slot.name = 'netP')
# 
# #visualize
# netAnalysis_signalingRole_heatmap(object = Greenleaf_CellChat,pattern = 'outgoing')
# 
# #do not know what is wrong here.

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
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = Greenleaf_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(Greenleaf_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 8/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 H marker mediated interaction') + 
  scale_x_discrete(labels = levels(Greenleaf_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = Greenleaf_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(Greenleaf_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 24/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(Greenleaf_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_1_interaction_by_H_marker.pdf',width = 8,height = 12)
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
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = Greenleaf_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(Greenleaf_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 9/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 HR marker mediated interaction') + 
  scale_x_discrete(labels = levels(Greenleaf_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = Greenleaf_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(Greenleaf_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 14/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(Greenleaf_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_1_interaction_by_HR_marker.pdf',width = 8,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

## conserved in all species ----------------------------------------------------
gene_list <- RG_1_marker_list$`conserved-in-HRM`
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#visualize
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = Greenleaf_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(Greenleaf_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 13/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 HRM marker mediated interaction') + 
  scale_x_discrete(labels = levels(Greenleaf_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = Greenleaf_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(Greenleaf_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 18/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(Greenleaf_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_97_fig_230130/Greenleaf_cellchat_RG_1_interaction_by_HRM_marker.pdf',width = 8,height = 12)
p1+p2+plot_layout(ncol = 1)
dev.off()