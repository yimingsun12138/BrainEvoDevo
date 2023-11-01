#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: species cellchat analysis focused on LGALS3                     ##
## Data: 2023.03.01                                                                ##
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

#RNA-seq data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#parameter
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230211.rds')

# check LGALS3 and LGALS3BP gene expression -------------------------------

#human
p1 <- DimPlot(object = Greenleaf_RNA_Seurat,cols = col_param$celltype,pt.size = 0.001,group.by = 'ReAnno_celltype',label = FALSE) + 
  theme_bw() + theme(aspect.ratio = 1) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p2 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'LGALS3',cols = ArchRPalettes$whitePurple[-1],order = TRUE,pt.size = 0.001) + 
  theme_bw() + theme(aspect.ratio = 1) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p3 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'LGALS3BP',cols = ArchRPalettes$whitePurple[-1],order = TRUE,pt.size = 0.001) + 
  theme_bw() + theme(aspect.ratio = 1) + theme(plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_105_fig_230301/human_LGALS3_gene_expression.pdf',width = 12,height = 5)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#macaque
p1 <- DimPlot(object = macaque_multiome_Seurat,cols = col_param$celltype,pt.size = 0.001,group.by = 'cell_type',label = FALSE) + 
  theme_bw() + theme(aspect.ratio = 1) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'LGALS3',cols = ArchRPalettes$whitePurple[-1],order = TRUE,pt.size = 0.001) + 
  theme_bw() + theme(aspect.ratio = 1) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'LGALS3BP',cols = ArchRPalettes$whitePurple[-1],order = TRUE,pt.size = 0.001) + 
  theme_bw() + theme(aspect.ratio = 1) + theme(plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_105_fig_230301/macaque_LGALS3_gene_expression.pdf',width = 12,height = 5)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#mouse
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
mouse_multiome_Seurat$ReAnno_celltype <- meta_data[colnames(mouse_multiome_Seurat),"ReAnno_celltype"]
p1 <- DimPlot(object = mouse_multiome_Seurat,cols = col_param$celltype,pt.size = 0.001,group.by = 'ReAnno_celltype',label = FALSE) + 
  theme_bw() + theme(aspect.ratio = 1) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p2 <- FeaturePlot(object = mouse_multiome_Seurat,features = 'Lgals3',cols = ArchRPalettes$whitePurple[-1],order = TRUE,pt.size = 0.001) + 
  theme_bw() + theme(aspect.ratio = 1) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p3 <- FeaturePlot(object = mouse_multiome_Seurat,features = 'Lgals3bp',cols = ArchRPalettes$whitePurple[-1],order = TRUE,pt.size = 0.001) + 
  theme_bw() + theme(aspect.ratio = 1) + theme(plot.title = element_text(hjust = 0.5))

pdf(file = './res/step_105_fig_230301/mouse_LGALS3_gene_expression.pdf',width = 12,height = 5)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#seems both macaque and mouse do not express LGALS3 and LGALS3BP

# human cellchat analysis -------------------------------------------------

#create cellchat object
express_matrix <- Greenleaf_RNA_Seurat@assays$RNA@data
meta_data <- Greenleaf_RNA_Seurat@meta.data
Greenleaf_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'ReAnno_celltype')

table(Greenleaf_CellChat@idents)
cell_type_list <- names(col_param$celltype)[names(col_param$celltype) %in% Greenleaf_CellChat@meta$ReAnno_celltype]
Greenleaf_CellChat <- setIdent(object = Greenleaf_CellChat,ident.use = 'ReAnno_celltype',levels = cell_type_list,display.warning = TRUE)

#set LR database
CellChatDB.update <- CellChatDB.human
temp <- data.frame(interaction_name = 'LGALS3_LGALS3BP',pathway_name = 'customized',ligand = 'LGALS3',receptor = 'LGALS3BP',
                   agonist = '',antagonist = '',co_A_receptor = '',co_I_receptor = '',evidence = 'PMID: 34728600',annotation = 'ECM-Receptor',
                   interaction_name_2 = 'LGALS3 - LGALS3BP')
rownames(temp) <- temp$interaction_name
CellChatDB.update$interaction <- rbind(CellChatDB.update$interaction,temp)
Greenleaf_CellChat@DB <- CellChatDB.update

#preprocess expression data
Greenleaf_CellChat <- subsetData(object = Greenleaf_CellChat)
future::plan(strategy = 'multisession',workers = 5)
future::plan(strategy = 'sequential')

#identify overexpression genes and interactions
Greenleaf_CellChat <- identifyOverExpressedGenes(object = Greenleaf_CellChat)
Greenleaf_CellChat <- identifyOverExpressedInteractions(object = Greenleaf_CellChat)

#infer cell cell communication network
Greenleaf_CellChat <- computeCommunProb(object = Greenleaf_CellChat,raw.use = TRUE)
Greenleaf_CellChat <- filterCommunication(object = Greenleaf_CellChat,min.cells = 10)

#extract inferred cell cell communication network
df.net <- subsetCommunication(object = Greenleaf_CellChat)

#infer cell cell communication level at signaling pathway level
Greenleaf_CellChat <- computeCommunProbPathway(object = Greenleaf_CellChat)

#calculate aggregated cell cell communication network
Greenleaf_CellChat <- aggregateNet(object = Greenleaf_CellChat)

#visualization
groupSize <- as.numeric(table(Greenleaf_CellChat@idents))
netVisual_circle(net = Greenleaf_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength",
                 color.use = col_param$celltype[levels(Greenleaf_CellChat@idents)])

#customized LGALS3 LGALS3BP signal
'customized' %in% df.net$pathway_name
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = 'customized',geneLR.return = FALSE)

#save results
pdf(file = './res/step_105_fig_230301/human_LGALS3_mediated_cell_cell_communication.pdf',width = 5,height = 6)
netVisual_individual(object = Greenleaf_CellChat,signaling = 'customized',pairLR.use = 'LGALS3_LGALS3BP',layout = 'circle',
                     color.use = col_param$celltype[levels(Greenleaf_CellChat@idents)])
dev.off()

save(Greenleaf_CellChat,file = '/content/data/sunym/temp/Greenleaf_CellChat_230301.RData')

# macaque cellchat analysis -----------------------------------------------

#create cellchat object
express_matrix <- macaque_multiome_Seurat@assays$RNA@data
meta_data <- macaque_multiome_Seurat@meta.data
macaque_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'cell_type')

table(macaque_CellChat@idents)
cell_type_list <- names(col_param$celltype)[names(col_param$celltype) %in% macaque_CellChat@meta$cell_type]
macaque_CellChat <- setIdent(object = macaque_CellChat,ident.use = 'cell_type',levels = cell_type_list,display.warning = TRUE)

#set LR database
CellChatDB.update <- CellChatDB.human
temp <- data.frame(interaction_name = 'LGALS3_LGALS3BP',pathway_name = 'customized',ligand = 'LGALS3',receptor = 'LGALS3BP',
                   agonist = '',antagonist = '',co_A_receptor = '',co_I_receptor = '',evidence = 'PMID: 34728600',annotation = 'ECM-Receptor',
                   interaction_name_2 = 'LGALS3 - LGALS3BP')
rownames(temp) <- temp$interaction_name
CellChatDB.update$interaction <- rbind(CellChatDB.update$interaction,temp)
macaque_CellChat@DB <- CellChatDB.update

#preprocess expression data
macaque_CellChat <- subsetData(object = macaque_CellChat)
future::plan(strategy = 'multisession',workers = 5)
future::plan(strategy = 'sequential')

#identify overexpression genes and interactions
macaque_CellChat <- identifyOverExpressedGenes(object = macaque_CellChat)
macaque_CellChat <- identifyOverExpressedInteractions(object = macaque_CellChat)

#infer cell cell communication network
macaque_CellChat <- computeCommunProb(object = macaque_CellChat,raw.use = TRUE)
macaque_CellChat <- filterCommunication(object = macaque_CellChat,min.cells = 10)

#extract inferred cell cell communication network
df.net <- subsetCommunication(object = macaque_CellChat)

#infer cell cell communication level at signaling pathway level
macaque_CellChat <- computeCommunProbPathway(object = macaque_CellChat)

#calculate aggregated cell cell communication network
macaque_CellChat <- aggregateNet(object = macaque_CellChat)

#visualization
groupSize <- as.numeric(table(macaque_CellChat@idents))
netVisual_circle(net = macaque_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength",
                 color.use = col_param$celltype[levels(macaque_CellChat@idents)])

#customized LGALS3 LGALS3BP signal
'customized' %in% df.net$pathway_name

#save results
save(macaque_CellChat,file = '/content/data/sunym/temp/macaque_CellChat_230301.RData')

# mouse cellchat analysis -------------------------------------------------

#add annotation
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
mouse_multiome_Seurat$ReAnno_celltype <- meta_data[colnames(mouse_multiome_Seurat),"ReAnno_celltype"]

#create cellchat object
express_matrix <- mouse_multiome_Seurat@assays$RNA@data
meta_data <- mouse_multiome_Seurat@meta.data
mouse_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'ReAnno_celltype')

table(mouse_CellChat@idents)
cell_type_list <- names(col_param$celltype)[names(col_param$celltype) %in% mouse_CellChat@meta$ReAnno_celltype]
mouse_CellChat <- setIdent(object = mouse_CellChat,ident.use = 'ReAnno_celltype',levels = cell_type_list,display.warning = TRUE)

#set LR database
CellChatDB.update <- CellChatDB.mouse
temp <- data.frame(interaction_name = 'LGALS3_LGALS3BP',pathway_name = 'customized',ligand = 'Lgals3',receptor = 'Lgals3bp',
                   agonist = '',antagonist = '',co_A_receptor = '',co_I_receptor = '',evidence = 'PMID: 34728600',annotation = 'ECM-Receptor',
                   interaction_name_2 = 'Lgals3 - Lgals3bp')
rownames(temp) <- temp$interaction_name
CellChatDB.update$interaction <- rbind(CellChatDB.update$interaction,temp)
temp <- which(CellChatDB.update[['interaction']]$ligand %in% c('H2-BI','H2-Ea-ps'))
CellChatDB.update[["interaction"]] <- CellChatDB.update[["interaction"]][-1*temp,]

mouse_CellChat@DB <- CellChatDB.update

#preprocess expression data
mouse_CellChat <- subsetData(object = mouse_CellChat)
future::plan(strategy = 'multisession',workers = 5)
future::plan(strategy = 'sequential')

#identify overexpression genes and interactions
mouse_CellChat <- identifyOverExpressedGenes(object = mouse_CellChat)
mouse_CellChat <- identifyOverExpressedInteractions(object = mouse_CellChat)

#infer cell cell communication network
mouse_CellChat <- computeCommunProb(object = mouse_CellChat,raw.use = TRUE)
mouse_CellChat <- filterCommunication(object = mouse_CellChat,min.cells = 10)

#extract inferred cell cell communication network
df.net <- subsetCommunication(object = mouse_CellChat)

#infer cell cell communication level at signaling pathway level
mouse_CellChat <- computeCommunProbPathway(object = mouse_CellChat)

#calculate aggregated cell cell communication network
mouse_CellChat <- aggregateNet(object = mouse_CellChat)

#visualization
groupSize <- as.numeric(table(mouse_CellChat@idents))
netVisual_circle(net = mouse_CellChat@net$weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge= FALSE,title.name = "Interaction strength",
                 color.use = col_param$celltype[levels(mouse_CellChat@idents)])

#customized LGALS3 LGALS3BP signal
'customized' %in% df.net$pathway_name

#save results
save(mouse_CellChat,file = '/content/data/sunym/temp/mouse_CellChat_230301.RData')
