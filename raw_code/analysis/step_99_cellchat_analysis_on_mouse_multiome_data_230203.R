#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: cellchat analysis on mouse multiome data                        ##
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
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#re-annotate mouse multiome data
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
meta_data <- meta_data[colnames(mouse_multiome_Seurat@assays$RNA@data),]

# create cellchat object --------------------------------------------------

#create cellchat object
express_matrix <- mouse_multiome_Seurat@assays$RNA@data
meta_data <- meta_data
mouse_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'ReAnno_celltype')

#set LR database
mouse_CellChat@DB <- CellChatDB.mouse

CellChatDB.use <- CellChatDB.mouse
temp <- which(CellChatDB.use[['interaction']]$ligand %in% c('H2-BI','H2-Ea-ps'))
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1*temp,]
mouse_CellChat@DB <- CellChatDB.use

# preprocess expression data ----------------------------------------------
mouse_CellChat <- subsetData(object = mouse_CellChat)
future::plan(strategy = 'multisession',workers = 5)
future::plan(strategy = 'sequential')

#identify overexpression genes and interactions
mouse_CellChat <- identifyOverExpressedGenes(object = mouse_CellChat)
mouse_CellChat <- identifyOverExpressedInteractions(object = mouse_CellChat)

#project to human PPI network
mouse_CellChat_PPI <- projectData(object = mouse_CellChat,adjMatrix = PPI.mouse)

# infer cell cell communication network -----------------------------------
mouse_CellChat <- computeCommunProb(object = mouse_CellChat,raw.use = TRUE)
mouse_CellChat <- filterCommunication(object = mouse_CellChat,min.cells = 10)

mouse_CellChat_PPI <- computeCommunProb(object = mouse_CellChat_PPI,raw.use = FALSE)
mouse_CellChat_PPI <- filterCommunication(object = mouse_CellChat_PPI,min.cells = 10)

#extract inferred cell cell communication network
df.net <- subsetCommunication(object = mouse_CellChat)
df.net.PPI <- subsetCommunication(object = mouse_CellChat_PPI)

#infer cell cell communication at signaling pathway level
mouse_CellChat <- computeCommunProbPathway(object = mouse_CellChat)
mouse_CellChat_PPI <- computeCommunProbPathway(object = mouse_CellChat_PPI)

#calculate aggregated cell cell communication network
mouse_CellChat <- aggregateNet(object = mouse_CellChat)
mouse_CellChat_PPI <- aggregateNet(object = mouse_CellChat_PPI)

#visualization
groupSize <- as.numeric(table(mouse_CellChat@idents))
temp <- c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic')
mouse_CellChat@idents <- factor(mouse_CellChat@idents,levels = temp)

netVisual_circle(net = mouse_CellChat@net$count,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Number of interactions")

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_cell_type_communication_strength.pdf',width = 8,height = 8)
netVisual_circle(net = mouse_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength")
dev.off()

netVisual_circle(net = mouse_CellChat_PPI@net$count,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Number of interactions")

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_PPI_cell_type_communication_strength.pdf',width = 8,height = 8)
netVisual_circle(net = mouse_CellChat_PPI@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength")
dev.off()

# pdgfd signal -------------------------------------------------------------
'Pdgfd' %in% df.net$ligand
temp <- df.net[df.net$ligand == 'Pdgfd',]
temp

#no RG-1 related cell cell communication

#visualize all PDGFD mediated cell cell communication
pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = 'PDGF',geneLR.return = FALSE)

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_PDGFD_signal_network.pdf',width = 5,height = 5)
netVisual_individual(object = mouse_CellChat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle')
dev.off()

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_PPI_PDGFD_signal_network.pdf',width = 5,height = 5)
netVisual_individual(object = mouse_CellChat_PPI,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle')
dev.off()

#I do not think PPI network is a very good idea

# cell communication related to RG-1 --------------------------------------

#global visualize
pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_RG_1_realated_cell_communication.pdf',width = 5,height = 5)
netVisual_circle(net = mouse_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as source",sources.use = 'RG-1')
netVisual_circle(net = mouse_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as target",targets.use = 'RG-1')
dev.off()

#PPI global visualize
netVisual_circle(net = mouse_CellChat_PPI@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as source",sources.use = 'RG-1')
netVisual_circle(net = mouse_CellChat_PPI@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,
                 title.name = "RG-1 as target",targets.use = 'RG-1')

#It seems that PPI projection would increase some cell cell communication in LR pairs that detect no interaction in raw data,
# but have very little impact on global cell cell communication pattern.

#the strongest cell communication related to RG-1
p <- netVisual_bubble(object = mouse_CellChat,sources.use = 'RG-1',remove.isolate = FALSE)

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_RG_1_sourced_significant_cell_communication.pdf',width = 8,height = 16)
p + theme_bw() + 
  theme(aspect.ratio = 2.5,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 sourced significant cell communication')
dev.off()

p <- netVisual_bubble(object = mouse_CellChat,targets.use = 'RG-1',remove.isolate = FALSE)

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_RG_1_targeted_significant_cell_communication.pdf',width = 8,height = 18)
p + theme_bw() + 
  theme(aspect.ratio = 3,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 targeted significant cell communication')
dev.off()

# overlap with RG-1 marker ------------------------------------------------
RG_1_marker_list <- readRDS(file = './res/step_97_fig_230130/RG_1_marker_list_from_liuyt.rds')

#convert to mouse gene symbol
gene_list <- unlist(RG_1_marker_list)
temp <- gene_list
gene_list <- stringr::str_to_title(string = gene_list)
names(gene_list) <- temp

gene_list[!(gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts))]

gene_list["SAMD4A"] <- 'Samd4'
gene_list["PEA15"] <- 'Pea15a'
gene_list["DARS1"] <- 'Dars'
gene_list["SYPL1"] <- 'Sypl'
gene_list["SPART"] <- 'Spg20'
gene_list["QKI"] <- 'Qk'

gene_dic <- gene_list

## specific in human -------------------------------------------------------
gene_list <- RG_1_marker_list$`specific-in-H`
gene_list <- gene_dic[gene_list]
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#visualize
pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = mouse_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(mouse_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 4/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 H marker mediated interaction') + 
  scale_x_discrete(labels = levels(mouse_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = mouse_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(mouse_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 3/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(mouse_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_RG_1_interaction_by_H_marker.pdf',width = 8,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

## conserved in primate ----------------------------------------------------
gene_list <- RG_1_marker_list$`conserved-in-HR`
gene_list <- gene_dic[gene_list]
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#visualize
pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = mouse_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(mouse_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 2/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 HR marker mediated interaction') + 
  scale_x_discrete(labels = levels(mouse_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = mouse_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(mouse_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 14/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(mouse_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_RG_1_interaction_by_HR_marker.pdf',width = 8,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

## conserved in species ----------------------------------------------------
gene_list <- RG_1_marker_list$`conserved-in-HRM`
gene_list <- gene_dic[gene_list]
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#visualize
pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

p1 <- netVisual_bubble(object = mouse_CellChat,pairLR.use = pairLR,sources.use = 'RG-1',targets.use = levels(mouse_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 12/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = 'RG-1 HRM marker mediated interaction') + 
  scale_x_discrete(labels = levels(mouse_CellChat@idents)) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab('RG-1 sourced') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

p2 <- netVisual_bubble(object = mouse_CellChat,pairLR.use = pairLR,targets.use = 'RG-1',sources.use = levels(mouse_CellChat@idents)) + 
  theme_bw() + 
  theme(aspect.ratio = 11/14,
        axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold')) + 
  labs(title = '') + 
  scale_x_discrete(labels = levels(mouse_CellChat@idents)) + 
  ylab('RG-1 targeted') + xlab('') + 
  theme(axis.title = element_text(size = 14,face = 'bold'))

pdf(file = './res/step_99_fig_230203/mouse_multiome_cellchat_RG_1_interaction_by_HRM_marker.pdf',width = 8,height = 12)
p1+p2+plot_layout(ncol = 1)
dev.off()

## specific in mouse -------------------------------------------------------
gene_list <- RG_1_marker_list$`specific-in-M`
gene_list <- gene_dic[gene_list]
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net[idx > 0,]

#Basicly mediate no cell cell communication.