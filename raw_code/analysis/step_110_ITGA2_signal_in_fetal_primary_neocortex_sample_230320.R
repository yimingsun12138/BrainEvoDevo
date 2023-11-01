#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: ITGA2 signal in fetal primary neocortex sample                  ##
## Data: 2023.03.20                                                                ##
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
library(OpenAI4R)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session()

# load data ---------------------------------------------------------------

#RNA-seq data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#parameter
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#re anno
Greenleaf_RNA_Seurat$cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
mouse_multiome_Seurat$cell_type <- meta_data[colnames(mouse_multiome_Seurat),"ReAnno_celltype"]

#unify annotation name
data_list <- c('Greenleaf_RNA_Seurat','macaque_multiome_Seurat','mouse_multiome_Seurat')
cell_type_list <- c('RG-1','RG-2','Ex-1','Ex-2','Ex-3','Ex-4')
cell_type_list_unify <- c('RG-neu','RG-glia','ExN','ExM','ExUp','ExDp')

for (i in 1:length(data_list)) {
  temp <- get(x = data_list[i])
  temp$unify_cell_type <- temp$cell_type
  
  for (j in 1:length(cell_type_list)) {
    temp@meta.data[temp$cell_type == cell_type_list[j],"unify_cell_type"] <- cell_type_list_unify[j]
  }
  
  assign(x = data_list[i],value = temp)
}

# #modify col_param
# temp <- col_param$celltype
# for (i in 1:length(cell_type_list)) {
#   idx <- which(names(temp) == cell_type_list[i])
#   names(temp)[idx] <- cell_type_list_unify[i]
# }
# 
# col_param$unify_cell_type <- temp
# saveRDS(object = col_param,file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

# ITGA2 related gene expression -------------------------------------------
data_list <- c('Greenleaf_RNA_Seurat','macaque_multiome_Seurat','mouse_multiome_Seurat')
for (i in 1:length(data_list)) {
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'cell_type',label = FALSE,raster = FALSE,cols = col_param$celltype) + 
    theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
    labs(title = data_list[i]) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()) + 
    NoLegend()
  
  assign(x = paste0('p',i),value = p)
  
  #feature plot
  gene_list <- c('ITGA2','COL1A1','COL1A2')
  if(i == 3){
    gene_list <- str_to_title(gene_list)
  }
  
  for (j in 1:length(gene_list)) {
    indi <- TRUE
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],pt.size = 0.001,order = indi,raster = FALSE,cols = ArchRPalettes$whitePurple[-1]) + 
      theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank()) + 
      NoLegend()
    assign(x = paste0('p',i+j*length(data_list)),value = p)
  }
}

char <- paste0('ggarrange(',paste('p',1:12,sep = '',collapse = ','),',nrow = 4,ncol = 3)')

png(filename = './res/step_110_fig_230320/species_primary_COL11A1_COL11A2_ITGA2_expression_featureplot.png',width = 1400,height = 1600)
eval(parse(text = char))
dev.off()

# focus on COL1A1-ITGA2 pair ---------------------------------------------

## human cellchat analysis -------------------------------------------------
#create cellchat object
express_matrix <- Greenleaf_RNA_Seurat@assays$RNA@data
meta_data <- Greenleaf_RNA_Seurat@meta.data
Greenleaf_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'unify_cell_type')

table(Greenleaf_CellChat@idents)
cell_type_list <- names(col_param$unify_cell_type)[names(col_param$unify_cell_type) %in% Greenleaf_CellChat@meta$unify_cell_type]
Greenleaf_CellChat <- setIdent(object = Greenleaf_CellChat,ident.use = 'unify_cell_type',levels = cell_type_list,display.warning = TRUE)

#set LR database
Greenleaf_CellChat@DB <- CellChatDB.human

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
                 color.use = col_param$unify_cell_type[levels(Greenleaf_CellChat@idents)])

#COL11A1-ITGA2 signal
df.net[grep(pattern = 'ITGA2',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = 'COLLAGEN',geneLR.return = FALSE)

#save results
pdf(file = './res/step_110_fig_230320/human_COL1A1_ITGA2_mediated_cell_cell_communication.pdf',width = 5,height = 6)
netVisual_individual(object = Greenleaf_CellChat,signaling = 'COLLAGEN',pairLR.use = 'COL1A1_ITGA2_ITGB1',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(Greenleaf_CellChat@idents)])
dev.off()

## macaque cellchat analysis -----------------------------------------------
#create cellchat object
express_matrix <- macaque_multiome_Seurat@assays$RNA@data
meta_data <- macaque_multiome_Seurat@meta.data
macaque_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'unify_cell_type')

table(macaque_CellChat@idents)
cell_type_list <- names(col_param$unify_cell_type)[names(col_param$unify_cell_type) %in% macaque_CellChat@meta$unify_cell_type]
macaque_CellChat <- setIdent(object = macaque_CellChat,ident.use = 'unify_cell_type',levels = cell_type_list,display.warning = TRUE)

#set LR database
macaque_CellChat@DB <- CellChatDB.human

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
                 color.use = col_param$unify_cell_type[levels(macaque_CellChat@idents)])

#COL11A1-ITGA2 signal
df.net[grep(pattern = 'ITGA2',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = 'COLLAGEN',geneLR.return = FALSE)

#save results
pdf(file = './res/step_110_fig_230320/macaque_COL1A1_ITGA2_mediated_cell_cell_communication.pdf',width = 5,height = 6)
netVisual_individual(object = macaque_CellChat,signaling = 'COLLAGEN',pairLR.use = 'COL1A1_ITGA2_ITGB1',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(macaque_CellChat@idents)])
dev.off()

## mouse cellchat analysis -------------------------------------------------

#create cellchat object
express_matrix <- mouse_multiome_Seurat@assays$RNA@data
meta_data <- mouse_multiome_Seurat@meta.data
mouse_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'unify_cell_type')

table(mouse_CellChat@idents)
cell_type_list <- names(col_param$unify_cell_type)[names(col_param$unify_cell_type) %in% mouse_CellChat@meta$unify_cell_type]
mouse_CellChat <- setIdent(object = mouse_CellChat,ident.use = 'unify_cell_type',levels = cell_type_list,display.warning = TRUE)

#set LR database
CellChatDB.update <- CellChatDB.mouse
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
                 color.use = col_param$unify_cell_type[levels(mouse_CellChat@idents)])

#COL11A1-ITGA2 signal
df.net[grep(pattern = 'Itga2',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = 'COLLAGEN',geneLR.return = FALSE)

#save results
pdf(file = './res/step_110_fig_230320/mouse_COL1A1_ITGA2_mediated_cell_cell_communication.pdf',width = 5,height = 6)
netVisual_individual(object = mouse_CellChat,signaling = 'COLLAGEN',pairLR.use = 'COL1A1_ITGA2_ITGB1',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(mouse_CellChat@idents)])
dev.off()

# generate LGALS3 ITGA2 PDGFD mediated cell cell communication ------------

## human cellchat analysis -------------------------------------------------
#create cellchat object
express_matrix <- Greenleaf_RNA_Seurat@assays$RNA@data
meta_data <- Greenleaf_RNA_Seurat@meta.data
Greenleaf_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'unify_cell_type')

table(Greenleaf_CellChat@idents)
cell_type_list <- names(col_param$unify_cell_type)[names(col_param$unify_cell_type) %in% Greenleaf_CellChat@meta$unify_cell_type]
Greenleaf_CellChat <- setIdent(object = Greenleaf_CellChat,ident.use = 'unify_cell_type',levels = cell_type_list,display.warning = TRUE)

#set LR database
CellChatDB.update <- CellChatDB.human
temp <- data.frame(interaction_name = 'LGALS3_LGALS3BP',pathway_name = 'customized',ligand = 'LGALS3',receptor = 'LGALS3BP',
                   agonist = '',antagonist = '',co_A_receptor = '',co_I_receptor = '',evidence = 'PMID: 34728600',annotation = 'ECM-Receptor',
                   interaction_name_2 = 'LGALS3 - LGALS3BP')
rownames(temp) <- temp$interaction_name
CellChatDB.update$interaction <- rbind(CellChatDB.update$interaction,temp)
CellChatDB.update$interaction[grep(pattern = 'ITGA2_',x = CellChatDB.update$interaction$receptor),"pathway_name"] <- 'ITGA2'

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
                 color.use = col_param$unify_cell_type[levels(Greenleaf_CellChat@idents)])

#save results

#ITGA2
pdf(file = './res/step_110_fig_230320/human_ITGA2_signal_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_aggregate(object = Greenleaf_CellChat,signaling = c('ITGA2'),layout = 'circle',
                    color.use = col_param$unify_cell_type[levels(Greenleaf_CellChat@idents)])
dev.off()

#LGALS3
df.net[grep(pattern = 'LGALS3',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = 'customized',geneLR.return = FALSE)
pdf(file = './res/step_110_fig_230320/human_LGALS3_LGALS3BP_mediated_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_individual(object = Greenleaf_CellChat,signaling = 'customized',pairLR.use = 'LGALS3_LGALS3BP',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(Greenleaf_CellChat@idents)])
dev.off()

#PDGFD
df.net[grep(pattern = 'PDGFD',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = 'PDGF',geneLR.return = FALSE)
pdf(file = './res/step_110_fig_230320/human_PDGFD_PDGFRB_mediated_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_individual(object = Greenleaf_CellChat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(Greenleaf_CellChat@idents)])
dev.off()

#save data
save(Greenleaf_CellChat,file = '/home/sunym/temp/Greenleaf_CellChat.RData')

## macaque cellchat analysis -----------------------------------------------
#create cellchat object
express_matrix <- macaque_multiome_Seurat@assays$RNA@data
meta_data <- macaque_multiome_Seurat@meta.data
macaque_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'unify_cell_type')

table(macaque_CellChat@idents)
cell_type_list <- names(col_param$unify_cell_type)[names(col_param$unify_cell_type) %in% macaque_CellChat@meta$unify_cell_type]
macaque_CellChat <- setIdent(object = macaque_CellChat,ident.use = 'unify_cell_type',levels = cell_type_list,display.warning = TRUE)

#set LR database
CellChatDB.update <- CellChatDB.human
temp <- data.frame(interaction_name = 'LGALS3_LGALS3BP',pathway_name = 'customized',ligand = 'LGALS3',receptor = 'LGALS3BP',
                   agonist = '',antagonist = '',co_A_receptor = '',co_I_receptor = '',evidence = 'PMID: 34728600',annotation = 'ECM-Receptor',
                   interaction_name_2 = 'LGALS3 - LGALS3BP')
rownames(temp) <- temp$interaction_name
CellChatDB.update$interaction <- rbind(CellChatDB.update$interaction,temp)
CellChatDB.update$interaction[grep(pattern = 'ITGA2_',x = CellChatDB.update$interaction$receptor),"pathway_name"] <- 'ITGA2'

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
                 color.use = col_param$unify_cell_type[levels(macaque_CellChat@idents)])

#save results

#ITGA2
pdf(file = './res/step_110_fig_230320/macaque_ITGA2_signal_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_aggregate(object = macaque_CellChat,signaling = c('ITGA2'),layout = 'circle',
                    color.use = col_param$unify_cell_type[levels(macaque_CellChat@idents)])
dev.off()

#LGALS3
df.net[grep(pattern = 'LGALS3',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = 'customized',geneLR.return = FALSE)
pdf(file = './res/step_110_fig_230320/macaque_LGALS3_LGALS3BP_mediated_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_individual(object = macaque_CellChat,signaling = 'customized',pairLR.use = 'LGALS3_LGALS3BP',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(macaque_CellChat@idents)])
dev.off()

#PDGFD
df.net[grep(pattern = 'PDGFD',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = 'PDGF',geneLR.return = FALSE)
pdf(file = './res/step_110_fig_230320/macaque_PDGFD_PDGFRB_mediated_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_individual(object = macaque_CellChat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(macaque_CellChat@idents)])
dev.off()

#save data
save(macaque_CellChat,file = '/home/sunym/temp/macaque_CellChat.RData')

## mouse cellchat analysis -------------------------------------------------

#create cellchat object
express_matrix <- mouse_multiome_Seurat@assays$RNA@data
meta_data <- mouse_multiome_Seurat@meta.data
mouse_CellChat <- createCellChat(object = express_matrix,meta = meta_data,group.by = 'unify_cell_type')

table(mouse_CellChat@idents)
cell_type_list <- names(col_param$unify_cell_type)[names(col_param$unify_cell_type) %in% mouse_CellChat@meta$unify_cell_type]
mouse_CellChat <- setIdent(object = mouse_CellChat,ident.use = 'unify_cell_type',levels = cell_type_list,display.warning = TRUE)

#set LR database
CellChatDB.update <- CellChatDB.mouse
temp <- data.frame(interaction_name = 'LGALS3_LGALS3BP',pathway_name = 'customized',ligand = 'Lgals3',receptor = 'Lgals3bp',
                   agonist = '',antagonist = '',co_A_receptor = '',co_I_receptor = '',evidence = 'PMID: 34728600',annotation = 'ECM-Receptor',
                   interaction_name_2 = 'Lgals3 - Lgals3bp')
rownames(temp) <- temp$interaction_name
CellChatDB.update$interaction <- rbind(CellChatDB.update$interaction,temp)
CellChatDB.update$interaction[grep(pattern = 'ITGA2_',x = CellChatDB.update$interaction$receptor),"pathway_name"] <- 'ITGA2'
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
                 color.use = col_param$unify_cell_type[levels(mouse_CellChat@idents)])

#save results

#ITGA2
pdf(file = './res/step_110_fig_230320/mouse_ITGA2_signal_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_aggregate(object = mouse_CellChat,signaling = c('ITGA2'),layout = 'circle',
                    color.use = col_param$unify_cell_type[levels(mouse_CellChat@idents)])
dev.off()

#LGALS3
df.net[grep(pattern = 'Lgals3',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = 'customized',geneLR.return = FALSE)
pdf(file = './res/step_110_fig_230320/mouse_LGALS3_LGALS3BP_mediated_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_individual(object = mouse_CellChat,signaling = 'customized',pairLR.use = 'LGALS3_LGALS3BP',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(mouse_CellChat@idents)])
dev.off()

#PDGFD
df.net[grep(pattern = 'Pdgf',x = df.net$interaction_name_2,fixed = TRUE),]
pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = 'PDGF',geneLR.return = FALSE)
pdf(file = './res/step_110_fig_230320/mouse_PDGFD_PDGFRB_mediated_cell_cell_communication.pdf',width = 5,height = 10)
netVisual_individual(object = mouse_CellChat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle',
                     color.use = col_param$unify_cell_type[levels(mouse_CellChat@idents)])
dev.off()

#save data
save(mouse_CellChat,file = '/home/sunym/temp/mouse_CellChat.RData')

# generate detailed LR interaction ----------------------------------------

#load data
load(file = '/home/sunym/temp/Greenleaf_CellChat.RData')
load(file = '/home/sunym/temp/macaque_CellChat.RData')
load(file = '/home/sunym/temp/mouse_CellChat.RData')

df.net.human <- subsetCommunication(object = Greenleaf_CellChat)
df.net.macaque <- subsetCommunication(object = macaque_CellChat)
df.net.mouse <- subsetCommunication(object = mouse_CellChat)

temp_fun <- function(dataset,interaction_list){
  interaction_iso <- interaction_list[!(interaction_list %in% dataset$interaction_name_2)]
  if(length(interaction_iso) == 0){
    dataset$interaction_name_2 <- factor(dataset$interaction_name_2,levels = interaction_list)
    return(dataset)
  }else{
    temp <- matrix(data = NA,nrow = length(interaction_iso),ncol = ncol(dataset))
    temp <- as.data.frame(temp)
    colnames(temp) <- colnames(dataset)
    temp$source <- 'RG-neu'
    temp$target <- 'RG-neu'
    temp$pval <- 1
    temp$interaction_name_2 <- interaction_iso
    temp$source.target <- 'RG-neu -> RG-neu'
    dataset <- rbind(dataset,temp)
    dataset$interaction_name_2 <- factor(dataset$interaction_name_2,levels = interaction_list)
    return(dataset)
  }
}

## human specific marker ---------------------------------------------------
gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')

#human
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.human$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.human[idx > 0,]

pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_human <- netVisual_bubble(object = Greenleaf_CellChat,sources.use = 'RG-neu',targets.use = levels(Greenleaf_CellChat@idents),pairLR.use = pairLR,return.data = TRUE)
target_human <- netVisual_bubble(object = Greenleaf_CellChat,sources.use = levels(Greenleaf_CellChat@idents),targets.use = 'RG-neu',pairLR.use = pairLR,return.data = TRUE)

#macaque
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.macaque$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.macaque[idx > 0,]

pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_macaque <- netVisual_bubble(object = macaque_CellChat,sources.use = 'RG-neu',targets.use = levels(macaque_CellChat@idents),pairLR.use = pairLR,return.data = TRUE)
target_macaque <- netVisual_bubble(object = macaque_CellChat,sources.use = levels(macaque_CellChat@idents),targets.use = 'RG-neu',pairLR.use = pairLR,return.data = TRUE)

#mouse
temp <- gene_list
gene_list <- stringr::str_to_title(string = gene_list)
names(gene_list) <- temp
gene_list[!(gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts))]

gene_list["SPART"] <- 'Spg20'
gene_list["SYPL1"] <- 'Sypl'
gene_list["DARS1"] <- 'Dars'
gene_list["PEA15"] <- 'Pea15a'
gene_list[!(gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts))]

idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.mouse$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.mouse[idx > 0,]

pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_mouse <- netVisual_bubble(object = mouse_CellChat,sources.use = 'RG-neu',targets.use = levels(mouse_CellChat@idents),pairLR.use = pairLR,return.data = TRUE)
target_mouse <- netVisual_bubble(object = mouse_CellChat,sources.use = levels(mouse_CellChat@idents),targets.use = 'RG-neu',pairLR.use = pairLR,return.data = TRUE)

#sourced data frame
source_human <- source_human$communication
source_human$interaction_name_2 <- as.character(source_human$interaction_name_2)

source_macaque <- source_macaque$communication
source_macaque$interaction_name_2 <- as.character(source_macaque$interaction_name_2)

source_mouse <- source_mouse$communication
source_mouse$interaction_name_2 <- as.character(source_mouse$interaction_name_2)
source_mouse$interaction_name_2 <- gsub(pattern = '  - ',replacement = ' - ',x = stringr::str_to_upper(string = source_mouse$interaction_name_2),fixed = TRUE)

interaction_name_list <- unique(c(source_human$interaction_name_2,source_macaque$interaction_name_2,source_mouse$interaction_name_2))
interaction_name_list <- na.omit(interaction_name_list)
interaction_name_list[!(interaction_name_list %in% Greenleaf_CellChat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_CellChat@DB$interaction$interaction_name_2[Greenleaf_CellChat@DB$interaction$interaction_name_2 %in% interaction_name_list]

source_human <- temp_fun(dataset = source_human,interaction_list = interaction_name_list)
source_macaque <- temp_fun(dataset = source_macaque,interaction_list = interaction_name_list)
source_mouse <- temp_fun(dataset = source_mouse,interaction_list = interaction_name_list)

values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")

df <- source_human
p1 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = '^RG-neu -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-neu sourced') + 
  labs(title = 'Human') + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 5/14)
if(sum(!is.na(df$prob)) <= 1){
  p1 <- p1+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white')
}else{
  p1 <- p1+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- source_macaque
p2 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = '^RG-neu -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = 'Macaque') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 5/15)
if(sum(!is.na(df$prob)) <= 1){
  p2 <- p2+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white')
}else{
  p2 <- p2+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- source_mouse
p3 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = '^RG-neu -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = 'Mouse') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 5/13)
if(sum(!is.na(df$prob)) <= 1){
  p3 <- p3+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white')
}else{
  p3 <- p3+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

#targeted data frame
target_human <- target_human$communication
target_human$interaction_name_2 <- as.character(target_human$interaction_name_2)

target_macaque <- target_macaque$communication
target_macaque$interaction_name_2 <- as.character(target_macaque$interaction_name_2)

target_mouse <- target_mouse$communication
target_mouse$interaction_name_2 <- as.character(target_mouse$interaction_name_2)
target_mouse$interaction_name_2 <- gsub(pattern = '  - ',replacement = ' - ',x = stringr::str_to_upper(string = target_mouse$interaction_name_2),fixed = TRUE)

interaction_name_list <- unique(c(target_human$interaction_name_2,target_macaque$interaction_name_2,target_mouse$interaction_name_2))
interaction_name_list <- na.omit(interaction_name_list)
interaction_name_list[!(interaction_name_list %in% Greenleaf_CellChat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_CellChat@DB$interaction$interaction_name_2[Greenleaf_CellChat@DB$interaction$interaction_name_2 %in% interaction_name_list]

target_human <- temp_fun(dataset = target_human,interaction_list = interaction_name_list)
target_macaque <- temp_fun(dataset = target_macaque,interaction_list = interaction_name_list)
target_mouse <- temp_fun(dataset = target_mouse,interaction_list = interaction_name_list)

values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")

df <- target_human
p4 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = ' -> RG-neu$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-neu targetd') + 
  labs(title = '') + 
  NoLegend() + 
  theme(aspect.ratio = 21/14)
if(sum(!is.na(df$prob)) <= 1){
  p4 <- p4+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white')
}else{
  p4 <- p4+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- target_macaque
p5 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = ' -> RG-neu$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = '') + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 21/15)
if(sum(!is.na(df$prob)) <= 1){
  p5 <- p5+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white')
}else{
  p5 <- p5+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- target_mouse
p6 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = ' -> RG-neu$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = '') + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 21/13)
if(sum(!is.na(df$prob)) <= 1){
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white')
}else{
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

pdf(file = './res/step_110_fig_230320/RG_1_human_specific_gene_mediated_cell_communication.pdf',width = 12,height = 12)
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
dev.off()

## species conserved marker ---------------------------------------------------
gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')

#human
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.human$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.human[idx > 0,]

pairLR <- extractEnrichedLR(object = Greenleaf_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_human <- netVisual_bubble(object = Greenleaf_CellChat,sources.use = 'RG-neu',targets.use = levels(Greenleaf_CellChat@idents),pairLR.use = pairLR,return.data = TRUE)
target_human <- netVisual_bubble(object = Greenleaf_CellChat,sources.use = levels(Greenleaf_CellChat@idents),targets.use = 'RG-neu',pairLR.use = pairLR,return.data = TRUE)

#macaque
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.macaque$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.macaque[idx > 0,]

pairLR <- extractEnrichedLR(object = macaque_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_macaque <- netVisual_bubble(object = macaque_CellChat,sources.use = 'RG-neu',targets.use = levels(macaque_CellChat@idents),pairLR.use = pairLR,return.data = TRUE)
target_macaque <- netVisual_bubble(object = macaque_CellChat,sources.use = levels(macaque_CellChat@idents),targets.use = 'RG-neu',pairLR.use = pairLR,return.data = TRUE)

#mouse
temp <- gene_list
gene_list <- stringr::str_to_title(string = gene_list)
names(gene_list) <- temp
gene_list[!(gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts))]

idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.mouse$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.mouse[idx > 0,]

pairLR <- extractEnrichedLR(object = mouse_CellChat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_mouse <- netVisual_bubble(object = mouse_CellChat,sources.use = 'RG-neu',targets.use = levels(mouse_CellChat@idents),pairLR.use = pairLR,return.data = TRUE)
target_mouse <- netVisual_bubble(object = mouse_CellChat,sources.use = levels(mouse_CellChat@idents),targets.use = 'RG-neu',pairLR.use = pairLR,return.data = TRUE)

#sourced data frame
source_human <- source_human$communication
source_human$interaction_name_2 <- as.character(source_human$interaction_name_2)

source_macaque <- source_macaque$communication
source_macaque$interaction_name_2 <- as.character(source_macaque$interaction_name_2)

source_mouse <- source_mouse$communication
source_mouse$interaction_name_2 <- as.character(source_mouse$interaction_name_2)
source_mouse$interaction_name_2 <- gsub(pattern = '  - ',replacement = ' - ',x = stringr::str_to_upper(string = source_mouse$interaction_name_2),fixed = TRUE)

interaction_name_list <- unique(c(source_human$interaction_name_2,source_macaque$interaction_name_2,source_mouse$interaction_name_2))
interaction_name_list <- na.omit(interaction_name_list)
interaction_name_list[!(interaction_name_list %in% Greenleaf_CellChat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_CellChat@DB$interaction$interaction_name_2[Greenleaf_CellChat@DB$interaction$interaction_name_2 %in% interaction_name_list]

source_human <- temp_fun(dataset = source_human,interaction_list = interaction_name_list)
source_macaque <- temp_fun(dataset = source_macaque,interaction_list = interaction_name_list)
source_mouse <- temp_fun(dataset = source_mouse,interaction_list = interaction_name_list)

values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")

df <- source_human
p1 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = '^RG-neu -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-neu sourced') + 
  labs(title = 'Human') + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 17/14)
if(sum(!is.na(df$prob)) <= 1){
  p1 <- p1+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white')
}else{
  p1 <- p1+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- source_macaque
p2 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = '^RG-neu -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = 'Macaque') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 17/15)
if(sum(!is.na(df$prob)) <= 1){
  p2 <- p2+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white')
}else{
  p2 <- p2+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- source_mouse
p3 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = '^RG-neu -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = 'Mouse') + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 17/13)
if(sum(!is.na(df$prob)) <= 1){
  p3 <- p3+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white')
}else{
  p3 <- p3+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

#targeted data frame
target_human <- target_human$communication
target_human$interaction_name_2 <- as.character(target_human$interaction_name_2)

target_macaque <- target_macaque$communication
target_macaque$interaction_name_2 <- as.character(target_macaque$interaction_name_2)

target_mouse <- target_mouse$communication
target_mouse$interaction_name_2 <- as.character(target_mouse$interaction_name_2)
target_mouse$interaction_name_2 <- gsub(pattern = '  - ',replacement = ' - ',x = stringr::str_to_upper(string = target_mouse$interaction_name_2),fixed = TRUE)

interaction_name_list <- unique(c(target_human$interaction_name_2,target_macaque$interaction_name_2,target_mouse$interaction_name_2))
interaction_name_list <- na.omit(interaction_name_list)
interaction_name_list[!(interaction_name_list %in% Greenleaf_CellChat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_CellChat@DB$interaction$interaction_name_2[Greenleaf_CellChat@DB$interaction$interaction_name_2 %in% interaction_name_list]

target_human <- temp_fun(dataset = target_human,interaction_list = interaction_name_list)
target_macaque <- temp_fun(dataset = target_macaque,interaction_list = interaction_name_list)
target_mouse <- temp_fun(dataset = target_mouse,interaction_list = interaction_name_list)

values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")

df <- target_human
p4 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = ' -> RG-neu$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-neu targetd') + 
  labs(title = '') + 
  NoLegend() + 
  theme(aspect.ratio = 22/14)
if(sum(!is.na(df$prob)) <= 1){
  p4 <- p4+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white')
}else{
  p4 <- p4+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- target_macaque
p5 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = ' -> RG-neu$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = '') + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 22/15)
if(sum(!is.na(df$prob)) <= 1){
  p5 <- p5+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white')
}else{
  p5 <- p5+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

df <- target_mouse
p6 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
  geom_point(pch = 16) + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  scale_x_discrete(labels = gsub(pattern = ' -> RG-neu$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('') + 
  labs(title = '') + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 22/13)
if(sum(!is.na(df$prob)) <= 1){
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white')
}else{
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

pdf(file = './res/step_110_fig_230320/RG_1_species_conserved_gene_mediated_cell_communication.pdf',width = 12,height = 12)
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
dev.off()