#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: species cellchat analysis focused on RG_1 marker                ##
## Data: 2023.02.11                                                                ##
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
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

Greenleaf_cellchat <- readRDS(file = '/content/data/sunym/temp/Greenleaf_cellchat.rds')
macaque_cellchat <- readRDS(file = '/content/data/sunym/temp/macaque_cellchat.rds')
mouse_cellchat <- readRDS(file = '/content/data/sunym/temp/mouse_cellchat.rds')

col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230211.rds')

# re-set idents -----------------------------------------------------------

#huamn
table(Greenleaf_cellchat@idents)
cell_type_list <- names(col_param$celltype)[names(col_param$celltype) %in% Greenleaf_cellchat@meta$ReAnno_celltype]
Greenleaf_cellchat <- setIdent(object = Greenleaf_cellchat,ident.use = 'ReAnno_celltype',levels = cell_type_list,display.warning = TRUE)

Greenleaf_cellchat <- filterCommunication(object = Greenleaf_cellchat,min.cells = 10)
Greenleaf_cellchat <- computeCommunProbPathway(object = Greenleaf_cellchat)
Greenleaf_cellchat <- aggregateNet(object = Greenleaf_cellchat)

#macaque
table(macaque_cellchat@idents)
cell_type_list <- names(col_param$celltype)[names(col_param$celltype) %in% macaque_cellchat@meta$cell_type]
macaque_cellchat <- setIdent(object = macaque_cellchat,ident.use = 'cell_type',levels = cell_type_list,display.warning = TRUE)

macaque_cellchat <- filterCommunication(object = macaque_cellchat,min.cells = 10)
macaque_cellchat <- computeCommunProbPathway(object = macaque_cellchat)
macaque_cellchat <- aggregateNet(object = macaque_cellchat)

#mouse
table(mouse_cellchat@idents)
cell_type_list <- levels(mouse_cellchat@idents)
mouse_cellchat <- setIdent(object = mouse_cellchat,ident.use = 'ReAnno_celltype',levels = cell_type_list,display.warning = TRUE)

mouse_cellchat <- filterCommunication(object = mouse_cellchat,min.cells = 10)
mouse_cellchat <- computeCommunProbPathway(object = mouse_cellchat)
mouse_cellchat <- aggregateNet(object = mouse_cellchat)

# modify mouse cell type color --------------------------------------------
# col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')
# 
# col_fun <- colorRamp2(breaks = c(0,1),colors = c(col_param$celltype['Cycling'],col_param$celltype['IP']))
# col_Cyc_IP <- col_fun(0.5)
# 
# col_param$celltype
# col_param$celltype <- c(col_param$celltype[c("RG","RG-1","RG-2","Cycling")],
#                         'Cyc-RG' = as.character(col_param$celltype["Cycling"]),
#                         'Cyc-IP' = as.character(col_Cyc_IP),
#                         col_param$celltype[c("IP","Ex-1","Ex-2","Ex-3","Ex-4")],
#                         'In' = as.character(col_param$celltype["InMGE"]),
#                         col_param$celltype[c("InMGE","InCGE","OPC","End","Per","VLMC","Mic")])
# 
# col_param$celltype.human <- NULL
# col_param$celltype.mouse <- NULL
# 
# saveRDS(object = col_param,file = './data/parameter/shared_param/MetaValue_color_param_230211.rds')

# PDGFD mediated cell cell communication ----------------------------------

#human
pairLR <- extractEnrichedLR(object = Greenleaf_cellchat,signaling = 'PDGF',geneLR.return = FALSE)
pdf(file = './res/step_102_fig_230211/Greenleaf_cellchat_PDGFD_mediated_cell_cell_communication.pdf',width = 5,height = 6)
netVisual_individual(object = Greenleaf_cellchat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle',
                     color.use = col_param$celltype[levels(Greenleaf_cellchat@idents)])
dev.off()

#macaque
pdf(file = './res/step_102_fig_230211/macaque_cellchat_PDGFD_mediated_cell_cell_communication.pdf',width = 5,height = 6)
netVisual_individual(object = macaque_cellchat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle',
                     color.use = col_param$celltype[levels(macaque_cellchat@idents)])
dev.off()

#mouse
pdf(file = './res/step_102_fig_230211/mouse_cellchat_PDGFD_mediated_cell_cell_communication.pdf',width = 5,height = 6)
netVisual_individual(object = mouse_cellchat,signaling = 'PDGF',pairLR.use = 'PDGFD_PDGFRB',layout = 'circle',
                     color.use = col_param$celltype[levels(mouse_cellchat@idents)])
dev.off()

# re-analysis RG-1 marker mediated cell communication ---------------------
df.net.human <- subsetCommunication(object = Greenleaf_cellchat)
df.net.macaque <- subsetCommunication(object = macaque_cellchat)
df.net.mouse <- subsetCommunication(object = mouse_cellchat)

temp_fun <- function(dataset,interaction_list){
  interaction_iso <- interaction_list[!(interaction_list %in% dataset$interaction_name_2)]
  if(length(interaction_iso) == 0){
    dataset$interaction_name_2 <- factor(dataset$interaction_name_2,levels = interaction_list)
    return(dataset)
  }else{
    temp <- matrix(data = NA,nrow = length(interaction_iso),ncol = ncol(dataset))
    temp <- as.data.frame(temp)
    colnames(temp) <- colnames(dataset)
    temp$source <- 'RG-1'
    temp$target <- 'RG-1'
    temp$pval <- 1
    temp$interaction_name_2 <- interaction_iso
    temp$source.target <- 'RG-1 -> RG-1'
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

pairLR <- extractEnrichedLR(object = Greenleaf_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_human <- netVisual_bubble(object = Greenleaf_cellchat,sources.use = 'RG-1',targets.use = levels(Greenleaf_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_human <- netVisual_bubble(object = Greenleaf_cellchat,sources.use = levels(Greenleaf_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

#macaque
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.macaque$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.macaque[idx > 0,]

pairLR <- extractEnrichedLR(object = macaque_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_macaque <- netVisual_bubble(object = macaque_cellchat,sources.use = 'RG-1',targets.use = levels(macaque_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_macaque <- netVisual_bubble(object = macaque_cellchat,sources.use = levels(macaque_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

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

pairLR <- extractEnrichedLR(object = mouse_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_mouse <- netVisual_bubble(object = mouse_cellchat,sources.use = 'RG-1',targets.use = levels(mouse_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_mouse <- netVisual_bubble(object = mouse_cellchat,sources.use = levels(mouse_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

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
interaction_name_list[!(interaction_name_list %in% Greenleaf_cellchat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_cellchat@DB$interaction$interaction_name_2[Greenleaf_cellchat@DB$interaction$interaction_name_2 %in% interaction_name_list]

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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-1 sourced') + 
  labs(title = 'Human') + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 4/14)
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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 4/15)
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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 4/13)
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
interaction_name_list[!(interaction_name_list %in% Greenleaf_cellchat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_cellchat@DB$interaction$interaction_name_2[Greenleaf_cellchat@DB$interaction$interaction_name_2 %in% interaction_name_list]

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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-1 targetd') + 
  labs(title = '') + 
  NoLegend() + 
  theme(aspect.ratio = 20/14)
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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 20/15)
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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 20/13)
if(sum(!is.na(df$prob)) <= 1){
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white')
}else{
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

pdf(file = './res/step_102_fig_230211/RG_1_human_specific_gene_mediated_cell_communication.pdf',width = 12,height = 12)
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
dev.off()

## primate conserved marker ---------------------------------------------------
gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')

#human
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.human$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.human[idx > 0,]

pairLR <- extractEnrichedLR(object = Greenleaf_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_human <- netVisual_bubble(object = Greenleaf_cellchat,sources.use = 'RG-1',targets.use = levels(Greenleaf_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_human <- netVisual_bubble(object = Greenleaf_cellchat,sources.use = levels(Greenleaf_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

#macaque
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.macaque$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.macaque[idx > 0,]

pairLR <- extractEnrichedLR(object = macaque_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_macaque <- netVisual_bubble(object = macaque_cellchat,sources.use = 'RG-1',targets.use = levels(macaque_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_macaque <- netVisual_bubble(object = macaque_cellchat,sources.use = levels(macaque_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

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

pairLR <- extractEnrichedLR(object = mouse_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_mouse <- netVisual_bubble(object = mouse_cellchat,sources.use = 'RG-1',targets.use = levels(mouse_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_mouse <- netVisual_bubble(object = mouse_cellchat,sources.use = levels(mouse_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

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
interaction_name_list[!(interaction_name_list %in% Greenleaf_cellchat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_cellchat@DB$interaction$interaction_name_2[Greenleaf_cellchat@DB$interaction$interaction_name_2 %in% interaction_name_list]

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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-1 sourced') + 
  labs(title = 'Human') + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  NoLegend() + 
  theme(aspect.ratio = 7/14)
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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 7/15)
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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 7/13)
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
interaction_name_list[!(interaction_name_list %in% Greenleaf_cellchat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_cellchat@DB$interaction$interaction_name_2[Greenleaf_cellchat@DB$interaction$interaction_name_2 %in% interaction_name_list]

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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-1 targetd') + 
  labs(title = '') + 
  NoLegend() + 
  theme(aspect.ratio = 8/14)
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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 8/15)
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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  theme(aspect.ratio = 8/13)
if(sum(!is.na(df$prob)) <= 1){
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white')
}else{
  p6 <- p6+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white',
                                 limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
                                 labels = c('min','max'))
}

pdf(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_mediated_cell_communication.pdf',width = 12,height = 12)
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

pairLR <- extractEnrichedLR(object = Greenleaf_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_human <- netVisual_bubble(object = Greenleaf_cellchat,sources.use = 'RG-1',targets.use = levels(Greenleaf_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_human <- netVisual_bubble(object = Greenleaf_cellchat,sources.use = levels(Greenleaf_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

#macaque
idx <- base::do.call(what = cbind,args = base::lapply(X = gene_list,FUN = function(x){
  temp <- grepl(pattern = x,x = df.net.macaque$interaction_name_2,fixed = TRUE)
  return(temp)
}))
idx <- rowSums(idx)
temp <- df.net.macaque[idx > 0,]

pairLR <- extractEnrichedLR(object = macaque_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_macaque <- netVisual_bubble(object = macaque_cellchat,sources.use = 'RG-1',targets.use = levels(macaque_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_macaque <- netVisual_bubble(object = macaque_cellchat,sources.use = levels(macaque_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

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

pairLR <- extractEnrichedLR(object = mouse_cellchat,signaling = unique(temp$pathway_name),geneLR.return = FALSE)
pairLR <- pairLR[pairLR$interaction_name %in% temp$interaction_name,,drop = FALSE]

source_mouse <- netVisual_bubble(object = mouse_cellchat,sources.use = 'RG-1',targets.use = levels(mouse_cellchat@idents),pairLR.use = pairLR,return.data = TRUE)
target_mouse <- netVisual_bubble(object = mouse_cellchat,sources.use = levels(mouse_cellchat@idents),targets.use = 'RG-1',pairLR.use = pairLR,return.data = TRUE)

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
interaction_name_list[!(interaction_name_list %in% Greenleaf_cellchat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_cellchat@DB$interaction$interaction_name_2[Greenleaf_cellchat@DB$interaction$interaction_name_2 %in% interaction_name_list]

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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-1 sourced') + 
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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
interaction_name_list[!(interaction_name_list %in% Greenleaf_cellchat@DB$interaction$interaction_name_2)]
interaction_name_list <- Greenleaf_cellchat@DB$interaction$interaction_name_2[Greenleaf_cellchat@DB$interaction$interaction_name_2 %in% interaction_name_list]

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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
  scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
               labels = names(values)[values %in% df$pval],name = 'p-value') + 
  guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
  geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
  theme(axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
  xlab('') + ylab('RG-1 targetd') + 
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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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
  scale_x_discrete(labels = gsub(pattern = ' -> RG-1$',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
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

pdf(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_mediated_cell_communication.pdf',width = 12,height = 12)
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
dev.off()

# #legend
# source_human <- source_human$communication
# source_human$interaction_name_2 <- as.character(source_human$interaction_name_2)
# 
# source_macaque <- source_macaque$communication
# source_macaque$interaction_name_2 <- as.character(source_macaque$interaction_name_2)
# 
# source_mouse <- source_mouse$communication
# source_mouse$interaction_name_2 <- as.character(source_mouse$interaction_name_2)
# source_mouse$interaction_name_2 <- gsub(pattern = '  - ',replacement = ' - ',x = stringr::str_to_upper(string = source_mouse$interaction_name_2),fixed = TRUE)
# 
# interaction_name_list <- unique(c(source_human$interaction_name_2,source_macaque$interaction_name_2,source_mouse$interaction_name_2))
# interaction_name_list <- na.omit(interaction_name_list)
# interaction_name_list[!(interaction_name_list %in% Greenleaf_cellchat@DB$interaction$interaction_name_2)]
# interaction_name_list <- Greenleaf_cellchat@DB$interaction$interaction_name_2[Greenleaf_cellchat@DB$interaction$interaction_name_2 %in% interaction_name_list]
# 
# source_human <- temp_fun(dataset = source_human,interaction_list = interaction_name_list)
# source_macaque <- temp_fun(dataset = source_macaque,interaction_list = interaction_name_list)
# source_mouse <- temp_fun(dataset = source_mouse,interaction_list = interaction_name_list)
# 
# values <- c(1, 2, 3)
# names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
# 
# df <- source_human
# p1 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
#   geom_point(pch = 16) + theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
#   scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
#   scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
#                labels = names(values)[values %in% df$pval],name = 'p-value') + 
#   guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
#   geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
#   geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
#   theme(axis.title = element_text(face = 'bold',size = 14),
#         plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
#   xlab('') + ylab('RG-1 sourced') + 
#   labs(title = 'Human') + 
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank()) + 
#   theme(aspect.ratio = 17/14)
# if(sum(!is.na(df$prob)) <= 1){
#   p1 <- p1+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white')
# }else{
#   p1 <- p1+scale_color_gradientn(colours = c('lightgrey',col_param$species['human']),na.value = 'white',
#                                  limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
#                                  breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
#                                  labels = c('min','max'))
# }
# 
# df <- source_macaque
# p2 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
#   geom_point(pch = 16) + theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
#   scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
#   scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
#                labels = names(values)[values %in% df$pval],name = 'p-value') + 
#   guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
#   geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
#   geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
#   theme(axis.title = element_text(face = 'bold',size = 14),
#         plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
#   xlab('') + ylab('') + 
#   labs(title = 'Macaque') + 
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank()) + 
#   theme(aspect.ratio = 17/15)
# if(sum(!is.na(df$prob)) <= 1){
#   p2 <- p2+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white')
# }else{
#   p2 <- p2+scale_color_gradientn(colours = c('lightgrey',col_param$species['macaque']),na.value = 'white',
#                                  limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
#                                  breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
#                                  labels = c('min','max'))
# }
# 
# df <- source_mouse
# p3 <- ggplot(data = df,mapping = aes(x = source.target,y = interaction_name_2,color = prob,size = pval)) + 
#   geom_point(pch = 16) + theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
#   scale_x_discrete(labels = gsub(pattern = '^RG-1 -> ',replacement = '',x = levels(df$source.target),fixed = FALSE),position = 'bottom') + 
#   scale_radius(range = c(min(df$pval),max(df$pval)),breaks = sort(unique(df$pval)),
#                labels = names(values)[values %in% df$pval],name = 'p-value') + 
#   guides(color = guide_colourbar(barwidth = 0.5,title = "Commun. Prob.")) + 
#   geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
#   geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 0.5, 1), lwd = 0.1, colour = 'lightgrey') + 
#   theme(axis.title = element_text(face = 'bold',size = 14),
#         plot.title = element_text(face = 'bold',size = 16,hjust = 0.5)) + 
#   xlab('') + ylab('') + 
#   labs(title = 'Mouse') + 
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank()) + 
#   theme(aspect.ratio = 17/13)
# if(sum(!is.na(df$prob)) <= 1){
#   p3 <- p3+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white')
# }else{
#   p3 <- p3+scale_color_gradientn(colours = c('lightgrey',col_param$species['mouse']),na.value = 'white',
#                                  limits = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
#                                  breaks = c(quantile(x = df$prob,0,na.rm = TRUE),quantile(x = df$prob,1,na.rm = TRUE)),
#                                  labels = c('min','max'))
# }
# 
# pdf(file = './res/step_102_fig_230211/RG_1_marker_netvisual_bubble_legend.pdf',width = 14,height = 5)
# p1+p2+p3+plot_layout(ncol = 3)
# dev.off()

## review certain genes ----------------------------------------------------

#BMP7 mediated cell cell communication seems odd
gene_list <- c('BMP7','BMPR1A','BMPR2','ACVR2B','ACVR1','BMPR2','ACVR2A')
for (i in 1:length(gene_list)) {
  p <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_list[i],pt.size = 0.001) + theme_bw() + theme(aspect.ratio = 1) + NoAxes() + NoLegend()
  assign(x = paste('p',i,sep = '_'),value = p)
}
p_1+p_2+p_3+p_4+p_5+p_6+p_7+plot_layout(ncol = 3)

for (i in 1:length(gene_list)) {
  p <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_list[i],pt.size = 0.001) + theme_bw() + theme(aspect.ratio = 1) + NoAxes() + NoLegend()
  assign(x = paste('p',i,sep = '_'),value = p)
}
p_1+p_2+p_3+p_4+p_5+p_6+p_7+plot_layout(ncol = 3)
#semms that BMP7 is not a very good marker in macaque (quite sparse)

#SPP1 − (ITGAV+ITGB5)
gene_list <- c('SPP1','ITGAV','ITGB5')
for (i in 1:length(gene_list)) {
  p <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_list[i],pt.size = 0.001) + theme_bw() + theme(aspect.ratio = 1) + NoAxes() + NoLegend()
  assign(x = paste('p',i,sep = '_'),value = p)
}
p_1+p_2+p_3+plot_layout(ncol = 3)

for (i in 1:length(gene_list)) {
  p <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_list[i],pt.size = 0.001) + theme_bw() + theme(aspect.ratio = 1) + NoAxes() + NoLegend()
  assign(x = paste('p',i,sep = '_'),value = p)
}
p_1+p_2+p_3+plot_layout(ncol = 3)

#ITGA5 seems not a very good RG-1 primate conserved marker