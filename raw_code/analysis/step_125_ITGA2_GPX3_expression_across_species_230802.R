#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: ITGA2_GPX3_expression_across_species                            ##
## Data: 2023.08.02                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.1/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(dplyr)
library(Seurat)
library(ggplot2)
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
library(paletteer)
library(ggpattern)
library(ggrepel)
library(openxlsx)
library(readxl)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0613')

# load data ---------------------------------------------------------------

#human data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PDhuman_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')

#chimp
Kanton_chimp_RNA_Seurat <- readRDS(file = './data/public/Organoid_single_cell_genomic_atlas_uncovers_human_specific_features_of_brain_development/processed_data/Chimp_filted_Seurat_object.rds')
Pollen_chimp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Co_RNA_Seurat_230302.rds')

#macaque
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_human_symbol_220802.rds')

#tree shrew
TreeShrew_RNA_Seurat <- readRDS(file = './res/step_125_fig_230802/TreeShrew_data/Treeshrew_SN1SN2_merged_Seurat_Preprocess_PC30_ScaleCountSample_AddAnno.rds')

#mouse
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_human_symbol_221009.rds')


# unify anno --------------------------------------------------------------
table(Greenleaf_RNA_Seurat$ReAnno_celltype)
Greenleaf_RNA_Seurat$unify_anno <- as.character(Greenleaf_RNA_Seurat$ReAnno_celltype)

table(PDhuman_RNA_Seurat$Cluster)
PDhuman_RNA_Seurat$unify_anno <- as.character(PDhuman_RNA_Seurat$Cluster)

table(Kanton_chimp_RNA_Seurat$PredCellType)
Kanton_chimp_RNA_Seurat$unify_anno <- as.character(Kanton_chimp_RNA_Seurat$PredCellType)

table(Pollen_chimp_RNA_Seurat$cell_type)
Pollen_chimp_RNA_Seurat$unify_anno <- as.character(Pollen_chimp_RNA_Seurat$cell_type)


table(macaque_multiome_Seurat$cell_type)
macaque_multiome_Seurat$unify_anno <- as.character(macaque_multiome_Seurat$cell_type)

table(TreeShrew_RNA_Seurat$predicted.celltype)
TreeShrew_RNA_Seurat$unify_anno <- as.character(TreeShrew_RNA_Seurat$predicted.celltype)

meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
mouse_multiome_Seurat$unify_anno <- as.character(meta_data[colnames(mouse_multiome_Seurat),"ReAnno_celltype"])
table(mouse_multiome_Seurat$unify_anno)

# normalize data ----------------------------------------------------------
label_list <- c('Greenleaf_RNA_Seurat' = 'Trevino et al.','PDhuman_RNA_Seurat' = 'Polioudakis et al.',
                'Kanton_chimp_RNA_Seurat' = 'Kanton et al.','Pollen_chimp_RNA_Seurat' = 'Pollen et al.',
                'macaque_multiome_Seurat' = 'macaque','TreeShrew_RNA_Seurat' = 'Tree Shrew',
                'mouse_multiome_Seurat' = 'mouse')
species_list <- c('Greenleaf_RNA_Seurat' = 'Human','PDhuman_RNA_Seurat' = 'Human',
                  'Kanton_chimp_RNA_Seurat' = 'Chimp','Pollen_chimp_RNA_Seurat' = 'Chimp',
                  'macaque_multiome_Seurat' = 'Macaque','TreeShrew_RNA_Seurat' = 'Tree Shrew',
                  'mouse_multiome_Seurat' = 'Mouse')

expression_matrix <- base::do.call(what = rbind,args = base::lapply(X = names(label_list),FUN = function(x){
  temp_data <- get(x = x)
  temp_data <- NormalizeData(object = temp_data,assay = 'RNA')
  if(!all(c('ITGA2','GPX3') %in% rownames(temp_data@assays$RNA@data))){
    stop('genes not found!')
  }
  temp_matrix <- temp_data@assays$RNA@data[c('ITGA2','GPX3'),temp_data$unify_anno %in% c('RG-1','RG','vRG','oRG')]
  colnames(temp_matrix) <- paste(x,colnames(temp_matrix),sep = '_')
  temp_matrix <- t(temp_matrix)
  temp_matrix <- as.data.frame(temp_matrix)
  temp_matrix$label <- label_list[x]
  temp_matrix$species <- species_list[x]
  return(temp_matrix)
  }))

expression_matrix$species <- factor(expression_matrix$species,levels = c('Human','Chimp','Macaque','Tree Shrew','Mouse'))
expression_matrix$label <- factor(expression_matrix$label,levels = as.character(label_list))

# plot --------------------------------------------------------------------
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

p1 <- ggplot(data = expression_matrix,mapping = aes(x = label,y = ITGA2,fill = species,color = species)) + 
  geom_violin() + 
  scale_fill_manual(values = c("#3FA4D9","#5EB6A5","#B2C224","#FF7F50","#8966A9")) + 
  scale_color_manual(values = c("#3FA4D9","#5EB6A5","#B2C224","#FF7F50","#8966A9")) + 
  theme_bw() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank()) + 
  xlab('Data source')

p2 <- ggplot(data = expression_matrix,mapping = aes(x = label,y = GPX3,fill = species,color = species)) + 
  geom_violin() + 
  scale_fill_manual(values = c("#3FA4D9","#5EB6A5","#B2C224","#FF7F50","#8966A9")) + 
  scale_color_manual(values = c("#3FA4D9","#5EB6A5","#B2C224","#FF7F50","#8966A9")) + 
  theme_bw() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid = element_blank()) + 
  xlab('Data source')

pdf(file = './res/step_125_fig_230802/ITGA2_GPX3_vlnplot.pdf',width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

#barplot
bar_table_ITGA2 <- base::do.call(what = rbind,args = base::lapply(X = as.character(label_list),FUN = function(x){
  mean_value <- mean(expression_matrix[expression_matrix$label == x,"ITGA2"])
  if(mean_value == 0){
    sd_value <- 0
  }else{
    sd_value <- sd(expression_matrix[expression_matrix$label == x,"ITGA2"])
  }
  temp <- data.frame(label = x,mean = mean_value,sd = sd_value)
}))
bar_table_ITGA2$species <- as.character(species_list)
bar_table_ITGA2$label <- factor(bar_table_ITGA2$label,levels = as.character(label_list))
bar_table_ITGA2$species <- factor(bar_table_ITGA2$species,levels = c('Human','Chimp','Macaque','Tree Shrew','Mouse'))

p1 <- ggplot(data = bar_table_ITGA2,mapping = aes(x = label,y = mean,fill = species)) + 
  geom_bar(stat = 'identity',position = position_dodge(width = 0.7),width = 0.6) + 
  geom_errorbar(mapping = aes(ymin = mean - 0.125*sd,ymax = mean + 0.125*sd),position = position_dodge(width = 0.7),width = 0.3,linewidth = 0.2) + 
  theme_bw() + 
  scale_fill_manual(values = c("#3FA4D9","#5EB6A5","#B2C224","#FF7F50","#8966A9")) + 
  theme(aspect.ratio = 1.5,
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = 'right') + 
  ylab('ITGA2')

bar_table_GPX3 <- base::do.call(what = rbind,args = base::lapply(X = as.character(label_list),FUN = function(x){
  mean_value <- mean(expression_matrix[expression_matrix$label == x,"GPX3"])
  if(mean_value == 0){
    sd_value <- 0
  }else{
    sd_value <- sd(expression_matrix[expression_matrix$label == x,"GPX3"])
  }
  temp <- data.frame(label = x,mean = mean_value,sd = sd_value)
}))
bar_table_GPX3$species <- as.character(species_list)
bar_table_GPX3$label <- factor(bar_table_GPX3$label,levels = as.character(label_list))
bar_table_GPX3$species <- factor(bar_table_GPX3$species,levels = c('Human','Chimp','Macaque','Tree Shrew','Mouse'))

p2 <- ggplot(data = bar_table_GPX3,mapping = aes(x = label,y = mean,fill = species)) + 
  geom_bar(stat = 'identity',position = position_dodge(width = 0.7),width = 0.6) + 
  geom_errorbar(mapping = aes(ymin = mean - 0.125*sd,ymax = mean + 0.125*sd),position = position_dodge(width = 0.7),width = 0.3,linewidth = 0.2) + 
  theme_bw() + 
  scale_fill_manual(values = c("#3FA4D9","#5EB6A5","#B2C224","#FF7F50","#8966A9")) + 
  theme(aspect.ratio = 1.5,
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = 'right') + 
  ylab('GPX3')

pdf(file = './res/step_125_fig_230802/ITGA2_GPX3_barplot_0.125sd.pdf',width = 8,height = 5)
p1+p2+plot_layout(ncol = 2)
dev.off()