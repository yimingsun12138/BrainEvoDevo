#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Fig_4_add_and_modify_figs                                       ##
## Data: 2023.04.05                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.2.3/')
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session()

# fig a -------------------------------------------------------------------

#load col_param
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#re anno mouse multiome Seurat
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
mouse_multiome_Seurat$cell_type <- meta_data[colnames(mouse_multiome_Seurat),"ReAnno_celltype"]

#get gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')

#' #GO extracellular matrix / extracellular matrix binding / extracellular matrix assembly / extracellular matrix organization / extracellular matrix-cell signaling
#' GO_term <- c('GO:0031012' = 'CC',
#'              'GO:0050840' = 'MF',
#'              #'GO:0085029' = 'BP',
#'              #'GO:0030198' = 'BP',
#'              'GO:0035426' = 'BP')
#' gene_list <- base::lapply(X = names(GO_term),FUN = function(x){
#'   temp <- getGOgeneSet(x = x,OrgDb = 'org.Hs.eg.db',ont = GO_term[x],keytype = 'SYMBOL')
#'   return(temp)
#' })
#' gene_list <- unique(unlist(gene_list))
#' 
#' #used gene list
#' gene_list <- c(HG_gene_list,PC_gene_list,SC_gene_list)[c(HG_gene_list,PC_gene_list,SC_gene_list) %in% gene_list]

gene_list <- c('ANXA2','ITGA2','LGALS3','LRIG1','ITGB5','BMP7','COL11A1','F3','COL4A5','SFRP1','MFGE8','GPC4','TNC')
gene_list <- c(gene_list[gene_list %in% SC_gene_list],gene_list[gene_list %in% PC_gene_list],gene_list[gene_list %in% HG_gene_list])

#plot
cell_type_list <- c('RG-1','IP','Ex-1','Ex-2','Ex-3','Ex-4')
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$ReAnno_celltype %in% cell_type_list]
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% cell_type_list]
mouse_multiome_Seurat <- mouse_multiome_Seurat[,mouse_multiome_Seurat$cell_type %in% cell_type_list]

p1 <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = gene_list,cols = c('lightgrey',col_param$species['human']),group.by = 'ReAnno_celltype')
p2 <- DotPlot(object = macaque_multiome_Seurat,assay = 'RNA',features = gene_list,cols = c('lightgrey',col_param$species['macaque']),group.by = 'cell_type')
p3 <- DotPlot(object = mouse_multiome_Seurat,assay = 'RNA',features = str_to_title(gene_list),cols = c('lightgrey',col_param$species['mouse']),group.by = 'cell_type')

scale_max <- max(c(p1$data$pct.exp,p2$data$pct.exp,p3$data$pct.exp))

p <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = gene_list,cols = c('lightgrey',col_param$species['human']),group.by = 'ReAnno_celltype',scale.min = 0,scale.max = scale_max)
p$data <- p$data[p$data$id %in% cell_type_list,]
p$data$id <- factor(p$data$id,levels = cell_type_list)

p1 <- p + coord_flip() + RotatedAxis() + 
  guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = 3,
        text = element_text(size = 10,family = 'sans')) + 
  labs(title = 'human')

pdf(file = './res/step_115_fig_230405/human_ECM_signature_dotplot.pdf',width = 5,height = 5)
p1
dev.off()

p <- DotPlot(object = macaque_multiome_Seurat,assay = 'RNA',features = gene_list,cols = c('lightgrey',col_param$species['macaque']),group.by = 'cell_type',scale.min = 0,scale.max = scale_max)
p$data <- p$data[p$data$id %in% cell_type_list,]
p$data$id <- factor(p$data$id,levels = cell_type_list)

p2 <- p + coord_flip() + RotatedAxis() + 
  guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = 3,
        text = element_text(size = 10,family = 'sans')) + 
  labs(title = 'macaque')

pdf(file = './res/step_115_fig_230405/macaque_ECM_signature_dotplot.pdf',width = 5,height = 5)
p2
dev.off()

p <- DotPlot(object = mouse_multiome_Seurat,assay = 'RNA',features = str_to_title(gene_list),cols = c('lightgrey',col_param$species['mouse']),group.by = 'cell_type',scale.min = 0,scale.max = scale_max)
p$data <- p$data[p$data$id %in% cell_type_list,]
p$data$id <- factor(p$data$id,levels = cell_type_list)

p3 <- p + coord_flip() + RotatedAxis() + 
  guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = 3,
        text = element_text(size = 10,family = 'sans')) + 
  labs(title = 'mouse')

pdf(file = './res/step_115_fig_230405/mouse_ECM_signature_dotplot.pdf',width = 5,height = 5)
p3
dev.off()

# lfc using Pollen data ---------------------------------------------------

#load data
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')

#normalize data
Hp_RNA_Seurat <- NormalizeData(Hp_RNA_Seurat)
Ho_RNA_Seurat <- NormalizeData(Ho_RNA_Seurat)

#re-name
Hp_RNA_Seurat <- RenameCells(object = Hp_RNA_Seurat,new.names = paste('Hp',colnames(Hp_RNA_Seurat),sep = '_'))
Ho_RNA_Seurat <- RenameCells(object = Ho_RNA_Seurat,new.names = paste('Ho',colnames(Ho_RNA_Seurat),sep = '_'))
table(rownames(Hp_RNA_Seurat@assays$RNA@counts) == rownames(Ho_RNA_Seurat@assays$RNA@counts))

#merge RG
RG_RNA_Seurat <- cbind(expm1(Hp_RNA_Seurat[,Hp_RNA_Seurat$cell_type == 'RG']@assays$RNA@data),
                       expm1(Ho_RNA_Seurat[,Ho_RNA_Seurat$cell_type == 'RG']@assays$RNA@data))
meta_data <- data.frame(cell_name = colnames(RG_RNA_Seurat),
                        dataset = c(rep('Primary',times = sum(Hp_RNA_Seurat$cell_type == 'RG')),rep('Organoid',times = sum(Ho_RNA_Seurat$cell_type == 'RG'))),
                        cell_type = 'RG-1')
rownames(meta_data) <- meta_data$cell_name
RG_RNA_Seurat <- CreateSeuratObject(counts = RG_RNA_Seurat,project = 'Pollen. et al.',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#normalize
RG_RNA_Seurat <- NormalizeData(RG_RNA_Seurat)

#calculate df
df_matrix <- my_DF_wilcox_test(mat1 = expm1(RG_RNA_Seurat[,RG_RNA_Seurat$dataset == 'Primary']@assays$RNA@data),
                               mat2 = expm1(RG_RNA_Seurat[,RG_RNA_Seurat$dataset == 'Organoid']@assays$RNA@data),
                               alternative = 'two.sided',paired = FALSE,workers = 4,future.globals.maxSize = 200*(1024^4))
df_matrix$sig <- -log10(df_matrix$fdr)

ggplot(data = df_matrix,mapping = aes(x = log2FC,y = sig)) + geom_point()
gene_list <- c('ANXA2','ITGA2','LGALS3','LRIG1','ITGB5','BMP7','COL11A1','F3','COL4A5','SFRP1','MFGE8','GPC4','TNC')
df_matrix[gene_list,]

#plot
df_matrix$color <- 'not_sig'
df_matrix[df_matrix$sig >= 3 & df_matrix$log2FC >= 1,"color"] <- 'pos'
df_matrix[df_matrix$sig >= 3 & df_matrix$log2FC <= -1,"color"] <- 'neg'
df_matrix$label <- ''
df_matrix[df_matrix$color == 'pos' & rownames(df_matrix) %in% gene_list,"label"] <- rownames(df_matrix[df_matrix$color == 'pos' & rownames(df_matrix) %in% gene_list,])
df_matrix[rownames(df_matrix) %in% df_matrix$label,"color"] <- 'labeled'
df_matrix$alpha = 'default'
df_matrix[df_matrix$color == 'labeled',"alpha"] <- 'labeled'

p <- ggplot(data = df_matrix,mapping = aes(x = log2FC,y = sig,color = color,alpha = alpha)) + 
  geom_point(data = df_matrix[df_matrix$alpha == 'default',],size = 2) + 
  geom_point(data = df_matrix[df_matrix$alpha == 'labeled',],size = 2,position = 'identity') + 
  scale_color_manual(values = c('not_sig' = 'lightgrey','neg' = '#0576C4','pos' = '#F18D1A','labeled' = 'red'),
                     limits = c('labeled','pos','neg','not_sig')) + 
  scale_alpha_manual(values = c('default' = 0.4,'labeled' = 1)) + 
  geom_hline(yintercept = 3,linewidth = 0.5,linetype = 'dashed') + 
  geom_vline(xintercept = -1,linewidth = 0.5,linetype = 'dashed') + 
  geom_vline(xintercept = 1,linewidth = 0.5,linetype = 'dashed') + 
  geom_label_repel(mapping = aes(label = label),size = 3,
                   box.padding = unit(0.5,'lines'),
                   point.padding = unit(0.8,'lines'),
                   max.overlaps = 10^10,
                   segment.color = 'black',
                   show.legend = FALSE) + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.ticks = element_blank())

pdf(file = './res/step_115_fig_230405/Pollen_RG_df_gene_plot.pdf',width = 6,height = 6)
p
dev.off()

# fig.e legend ------------------------------------------------------------

#load col param
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
PFC_RNA_Seurat <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_PFC_brain_preprocess_PC30_by_lyt_230326.rds')
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')
Paola_3_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_3m_Seurat.rds')
Paola_6_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_6m_Seurat.rds')
gc()

#preprocess PD human
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,npcs = 50,preprocess = TRUE)
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'Cluster',label = TRUE)

#normalize data
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in 1:length(data_list)) {
  temp <- get(data_list[i])
  temp <- NormalizeData(temp)
  assign(x = data_list[i],value = temp)
  rm(temp)
  gc()
}

meta_data <- readRDS(file = './res/step_113_fig_230326/PFC_RG_meta_data.rds')
#PFC_RNA_Seurat <- PFC_RNA_Seurat[,rownames(meta_data)]
PFC_RNA_Seurat@meta.data[rownames(meta_data),"cell.type"] <- meta_data$cell_type

#ori_anno
Greenleaf_RNA_Seurat$ori_anno <- Greenleaf_RNA_Seurat$ReAnno_celltype
PD_RNA_Seurat$ori_anno <- PD_RNA_Seurat$Cluster
PFC_RNA_Seurat$ori_anno <- PFC_RNA_Seurat$cell.type
Hp_RNA_Seurat$ori_anno <- Hp_RNA_Seurat$cell_type
Ho_RNA_Seurat$ori_anno <- Ho_RNA_Seurat$cell_type
Ho_1.5_Seurat$ori_anno <- Ho_1.5_Seurat$FinalName
Ho_1.5_Seurat$ori_anno <- as.character(Ho_1.5_Seurat$ori_anno)
Ho_2_Seurat$ori_anno <- Ho_2_Seurat$FinalName
Ho_2_Seurat$ori_anno <- as.character(Ho_2_Seurat$ori_anno)
Ho_3_Seurat$ori_anno <- Ho_3_Seurat$FinalName
Ho_3_Seurat$ori_anno <- as.character(Ho_3_Seurat$ori_anno)
Ho_4_Seurat$ori_anno <- Ho_4_Seurat$FinalName
Ho_4_Seurat$ori_anno <- as.character(Ho_4_Seurat$ori_anno)
Paola_3_Seurat$ori_anno <- Paola_3_Seurat$CellType
Paola_6_Seurat$ori_anno <- Paola_6_Seurat$CellType

#colnames
Greenleaf_RNA_Seurat <- RenameCells(object = Greenleaf_RNA_Seurat,new.names = paste('Greenleaf',colnames(Greenleaf_RNA_Seurat),sep = '_'))
PD_RNA_Seurat <- RenameCells(object = PD_RNA_Seurat,new.names = paste('PD',colnames(PD_RNA_Seurat),sep = '_'))
PFC_RNA_Seurat <- RenameCells(object = PFC_RNA_Seurat,new.names = paste('PFC',colnames(PFC_RNA_Seurat),sep = '_'))
Hp_RNA_Seurat <- RenameCells(object = Hp_RNA_Seurat,new.names = paste('Hp',colnames(Hp_RNA_Seurat),sep = '_'))
Ho_RNA_Seurat <- RenameCells(object = Ho_RNA_Seurat,new.names = paste('Ho',colnames(Ho_RNA_Seurat),sep = '_'))
Ho_1.5_Seurat <- RenameCells(object = Ho_1.5_Seurat,new.names = paste('Ho_1.5',colnames(Ho_1.5_Seurat),sep = '_'))
Ho_2_Seurat <- RenameCells(object = Ho_2_Seurat,new.names = paste('Ho_2',colnames(Ho_2_Seurat),sep = '_'))
Ho_3_Seurat <- RenameCells(object = Ho_3_Seurat,new.names = paste('Ho_3',colnames(Ho_3_Seurat),sep = '_'))
Ho_4_Seurat <- RenameCells(object = Ho_4_Seurat,new.names = paste('Ho_4',colnames(Ho_4_Seurat),sep = '_'))
Paola_3_Seurat <- RenameCells(object = Paola_3_Seurat,new.names = paste('Paola_3',colnames(Paola_3_Seurat),sep = '_'))
Paola_6_Seurat <- RenameCells(object = Paola_6_Seurat,new.names = paste('Paola_6',colnames(Paola_6_Seurat),sep = '_'))

Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1')]
PD_RNA_Seurat <- PD_RNA_Seurat[,PD_RNA_Seurat$ori_anno %in% c('oRG','vRG')]
PFC_RNA_Seurat <- PFC_RNA_Seurat[,PFC_RNA_Seurat$ori_anno %in% c('RG-neu')]
Hp_RNA_Seurat <- Hp_RNA_Seurat[,Hp_RNA_Seurat$cell_type == 'RG']
Ho_RNA_Seurat <- Ho_RNA_Seurat[,Ho_RNA_Seurat$cell_type == 'RG']
Ho_1.5_Seurat <- Ho_1.5_Seurat[,Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG')]
Ho_2_Seurat <- Ho_2_Seurat[,Ho_2_Seurat$ori_anno %in% c('aRG','oRG')]
Ho_3_Seurat <- Ho_3_Seurat[,Ho_3_Seurat$ori_anno %in% c('aRG','oRG')]
Ho_4_Seurat <- Ho_4_Seurat[,Ho_4_Seurat$ori_anno %in% c('aRG','oRG')]
Paola_3_Seurat <- Paola_3_Seurat[,Paola_3_Seurat$ori_anno %in% c('oRG','RG')]
Paola_6_Seurat <- Paola_6_Seurat[,Paola_6_Seurat$ori_anno %in% c('oRG','RG')]
gc()

gene_list <- Reduce(f = intersect,x = list(rownames(Greenleaf_RNA_Seurat),
                                           rownames(PD_RNA_Seurat),
                                           rownames(PFC_RNA_Seurat),
                                           rownames(Hp_RNA_Seurat),
                                           rownames(Ho_RNA_Seurat),
                                           rownames(Ho_1.5_Seurat),
                                           rownames(Paola_3_Seurat)))

RG_RNA_Seurat <- cbind(expm1(Greenleaf_RNA_Seurat@assays$RNA@data[gene_list,Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1')]),
                       expm1(PD_RNA_Seurat@assays$RNA@data[gene_list,PD_RNA_Seurat$ori_anno %in% c('oRG','vRG')]),
                       expm1(PFC_RNA_Seurat@assays$RNA@data[gene_list,PFC_RNA_Seurat$ori_anno %in% c('RG-neu')]),
                       expm1(Hp_RNA_Seurat@assays$RNA@data[gene_list,Hp_RNA_Seurat$cell_type == 'RG']),
                       expm1(Ho_RNA_Seurat@assays$RNA@data[gene_list,Ho_RNA_Seurat$cell_type == 'RG']),
                       expm1(Ho_1.5_Seurat@assays$RNA@data[gene_list,Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Ho_2_Seurat@assays$RNA@data[gene_list,Ho_2_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Ho_3_Seurat@assays$RNA@data[gene_list,Ho_3_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Ho_4_Seurat@assays$RNA@data[gene_list,Ho_4_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Paola_3_Seurat@assays$RNA@data[gene_list,Paola_3_Seurat$ori_anno %in% c('oRG','RG')]),
                       expm1(Paola_6_Seurat@assays$RNA@data[gene_list,Paola_6_Seurat$ori_anno %in% c('oRG','RG')]))

meta_data <- data.frame(cell_name = colnames(RG_RNA_Seurat),
                        dataset = c(rep('Trevino. et al.',times = sum(Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1'))),
                                    rep('Polioudakis. et al.',times = sum(PD_RNA_Seurat$ori_anno %in% c('oRG','vRG'))),
                                    rep('Bhaduri. et al.',times = sum(PFC_RNA_Seurat$ori_anno %in% c('RG-neu'))),
                                    rep('Pollen. et al. Primary',times = sum(Hp_RNA_Seurat$cell_type == 'RG')),
                                    rep('Pollen. et al. Organoid',times = sum(Ho_RNA_Seurat$cell_type == 'RG')),
                                    rep('Uzquiano. et al.',times = sum(c(Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG'),Ho_2_Seurat$ori_anno %in% c('aRG','oRG'),Ho_3_Seurat$ori_anno %in% c('aRG','oRG'),Ho_4_Seurat$ori_anno %in% c('aRG','oRG')))),
                                    rep('Velasco. et al.',times = sum(c(Paola_3_Seurat$ori_anno %in% c('oRG','RG'),Paola_6_Seurat$ori_anno %in% c('oRG','RG'))))))
rownames(meta_data) <- meta_data$cell_name

RG_RNA_Seurat <- CreateSeuratObject(counts = RG_RNA_Seurat,project = 'RG',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#add meta data
RG_RNA_Seurat$Age <- ''

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Greenleaf_RNA_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- Greenleaf_RNA_Seurat@meta.data[cell_list,"Age"]

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(PD_RNA_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- PD_RNA_Seurat@meta.data[cell_list,"Gestation_week"]

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(PFC_RNA_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- PFC_RNA_Seurat@meta.data[cell_list,"age"]

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_1.5_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '1.5m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_2_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '2m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_3_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '3m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_4_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '4m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Paola_3_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '3m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Paola_6_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '6m'

RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '17')] <- '17GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '18')] <- '18GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '19')] <- '19GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '20')] <- '20GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '22')] <- '22GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '25')] <- '25GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == 'pcw16')] <- '16GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == 'pcw20')] <- '20GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == 'pcw21')] <- '21GW'

RG_RNA_Seurat$sample <- paste(RG_RNA_Seurat$dataset,RG_RNA_Seurat$Age,sep = ' ')
RG_RNA_Seurat@meta.data[RG_RNA_Seurat$sample == 'Pollen. et al. Primary ',"sample"] <- 'Pollen. et al. Primary'
RG_RNA_Seurat@meta.data[RG_RNA_Seurat$sample == 'Pollen. et al. Organoid ',"sample"] <- 'Pollen. et al. Organoid'

RG_RNA_Seurat$group <- 'Primary'
RG_RNA_Seurat@meta.data[RG_RNA_Seurat$dataset %in% c('Uzquiano. et al.','Velasco. et al.'),"group"] <- 'Organoid'
RG_RNA_Seurat@meta.data[colnames(RG_RNA_Seurat) %in% colnames(Ho_RNA_Seurat),"group"] <- 'Organoid'

#normalize
RG_RNA_Seurat <- NormalizeData(RG_RNA_Seurat)
gc()

#factor
RG_RNA_Seurat$dataset <- factor(RG_RNA_Seurat$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al. Primary','Pollen. et al. Organoid','Uzquiano. et al.','Velasco. et al.'))
RG_RNA_Seurat$sample <- factor(RG_RNA_Seurat$sample,levels = c('Trevino. et al. 16GW','Trevino. et al. 20GW','Trevino. et al. 21GW',
                                                               'Polioudakis. et al. 17GW','Polioudakis. et al. 18GW',
                                                               'Bhaduri. et al. 17GW','Bhaduri. et al. 18GW','Bhaduri. et al. 19GW','Bhaduri. et al. 20GW','Bhaduri. et al. 22GW','Bhaduri. et al. 25GW',
                                                               'Pollen. et al. Primary','Pollen. et al. Organoid',
                                                               'Uzquiano. et al. 1.5m','Uzquiano. et al. 2m','Uzquiano. et al. 3m','Uzquiano. et al. 4m',
                                                               'Velasco. et al. 3m','Velasco. et al. 6m'))
RG_RNA_Seurat$group <- factor(RG_RNA_Seurat$group,levels = c('Primary','Organoid'))

#add author
RG_RNA_Seurat$Author <- as.character(RG_RNA_Seurat$dataset)
RG_RNA_Seurat@meta.data[RG_RNA_Seurat$dataset %in% c('Pollen. et al. Primary','Pollen. et al. Organoid'),"Author"] <- 'Pollen. et al.'
RG_RNA_Seurat$Author <- factor(RG_RNA_Seurat$Author,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
table(RG_RNA_Seurat$Author)

#plot legend
pdf(file = './res/step_115_fig_230405/public_dataset_color_legend.pdf',width = 8,height = 6)
VlnPlot(object = RG_RNA_Seurat,features = 'nCount_RNA',
        cols = c('#FF4136','#FF851B','#FFDC00','#39CCCC','#0074D9','#2ECC40'),
        group.by = 'Author') + theme(aspect.ratio = 1)
dev.off()

# ECM signature -----------------------------------------------------------

#load col param
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
PFC_RNA_Seurat <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_PFC_brain_preprocess_PC30_by_lyt_230326.rds')
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')
Paola_3_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_3m_Seurat.rds')
Paola_6_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_6m_Seurat.rds')
gc()

#preprocess PD human
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,npcs = 50,preprocess = TRUE)
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'Cluster',label = TRUE)

#normalize data
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in 1:length(data_list)) {
  temp <- get(data_list[i])
  temp <- NormalizeData(temp)
  assign(x = data_list[i],value = temp)
  rm(temp)
  gc()
}

meta_data <- readRDS(file = './res/step_113_fig_230326/PFC_RG_meta_data.rds')
#PFC_RNA_Seurat <- PFC_RNA_Seurat[,rownames(meta_data)]
PFC_RNA_Seurat@meta.data[rownames(meta_data),"cell.type"] <- meta_data$cell_type

#ori_anno
Greenleaf_RNA_Seurat$ori_anno <- Greenleaf_RNA_Seurat$ReAnno_celltype
PD_RNA_Seurat$ori_anno <- PD_RNA_Seurat$Cluster
PFC_RNA_Seurat$ori_anno <- PFC_RNA_Seurat$cell.type
Hp_RNA_Seurat$ori_anno <- Hp_RNA_Seurat$cell_type
Ho_RNA_Seurat$ori_anno <- Ho_RNA_Seurat$cell_type
Ho_1.5_Seurat$ori_anno <- Ho_1.5_Seurat$FinalName
Ho_1.5_Seurat$ori_anno <- as.character(Ho_1.5_Seurat$ori_anno)
Ho_2_Seurat$ori_anno <- Ho_2_Seurat$FinalName
Ho_2_Seurat$ori_anno <- as.character(Ho_2_Seurat$ori_anno)
Ho_3_Seurat$ori_anno <- Ho_3_Seurat$FinalName
Ho_3_Seurat$ori_anno <- as.character(Ho_3_Seurat$ori_anno)
Ho_4_Seurat$ori_anno <- Ho_4_Seurat$FinalName
Ho_4_Seurat$ori_anno <- as.character(Ho_4_Seurat$ori_anno)
Paola_3_Seurat$ori_anno <- Paola_3_Seurat$CellType
Paola_6_Seurat$ori_anno <- Paola_6_Seurat$CellType

#colnames
Greenleaf_RNA_Seurat <- RenameCells(object = Greenleaf_RNA_Seurat,new.names = paste('Greenleaf',colnames(Greenleaf_RNA_Seurat),sep = '_'))
PD_RNA_Seurat <- RenameCells(object = PD_RNA_Seurat,new.names = paste('PD',colnames(PD_RNA_Seurat),sep = '_'))
PFC_RNA_Seurat <- RenameCells(object = PFC_RNA_Seurat,new.names = paste('PFC',colnames(PFC_RNA_Seurat),sep = '_'))
Hp_RNA_Seurat <- RenameCells(object = Hp_RNA_Seurat,new.names = paste('Hp',colnames(Hp_RNA_Seurat),sep = '_'))
Ho_RNA_Seurat <- RenameCells(object = Ho_RNA_Seurat,new.names = paste('Ho',colnames(Ho_RNA_Seurat),sep = '_'))
Ho_1.5_Seurat <- RenameCells(object = Ho_1.5_Seurat,new.names = paste('Ho_1.5',colnames(Ho_1.5_Seurat),sep = '_'))
Ho_2_Seurat <- RenameCells(object = Ho_2_Seurat,new.names = paste('Ho_2',colnames(Ho_2_Seurat),sep = '_'))
Ho_3_Seurat <- RenameCells(object = Ho_3_Seurat,new.names = paste('Ho_3',colnames(Ho_3_Seurat),sep = '_'))
Ho_4_Seurat <- RenameCells(object = Ho_4_Seurat,new.names = paste('Ho_4',colnames(Ho_4_Seurat),sep = '_'))
Paola_3_Seurat <- RenameCells(object = Paola_3_Seurat,new.names = paste('Paola_3',colnames(Paola_3_Seurat),sep = '_'))
Paola_6_Seurat <- RenameCells(object = Paola_6_Seurat,new.names = paste('Paola_6',colnames(Paola_6_Seurat),sep = '_'))

#re-anno
Greenleaf_RNA_Seurat$re_anno <- Greenleaf_RNA_Seurat$ori_anno
PD_RNA_Seurat$re_anno <- PD_RNA_Seurat$ori_anno
PFC_RNA_Seurat$re_anno <- PFC_RNA_Seurat$ori_anno
Hp_RNA_Seurat$re_anno <- Hp_RNA_Seurat$ori_anno
Ho_RNA_Seurat$re_anno <- Ho_RNA_Seurat$ori_anno
Ho_1.5_Seurat$re_anno <- Ho_1.5_Seurat$ori_anno
Ho_2_Seurat$re_anno <- Ho_2_Seurat$ori_anno
Ho_3_Seurat$re_anno <- Ho_3_Seurat$ori_anno
Ho_4_Seurat$re_anno <- Ho_4_Seurat$ori_anno
Paola_3_Seurat$re_anno <- Paola_3_Seurat$ori_anno
Paola_6_Seurat$re_anno <- Paola_6_Seurat$ori_anno

Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1'),"re_anno"] <- 'RG-neu'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$ori_anno %in% c('oRG','vRG'),"re_anno"] <- 'RG-neu'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('RG-neu'),"re_anno"] <- 'RG-neu'
Hp_RNA_Seurat@meta.data[Hp_RNA_Seurat$cell_type == 'RG',"re_anno"] <- 'RG-neu'
Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type == 'RG',"re_anno"] <- 'RG-neu'
Ho_1.5_Seurat@meta.data[Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Ho_2_Seurat@meta.data[Ho_2_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Ho_3_Seurat@meta.data[Ho_3_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Ho_4_Seurat@meta.data[Ho_4_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Paola_3_Seurat@meta.data[Paola_3_Seurat$ori_anno %in% c('oRG','RG'),"re_anno"] <- 'RG-neu'
Paola_6_Seurat@meta.data[Paola_6_Seurat$ori_anno %in% c('oRG','RG'),"re_anno"] <- 'RG-neu'

#get gene list
gene_list <- getGOgeneSet(x = 'GO:0031012',OrgDb = org.Hs.eg.db,ont = 'CC',keytype = 'SYMBOL')
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat) & gene_list %in% rownames(PD_RNA_Seurat) & gene_list %in% rownames(PFC_RNA_Seurat) & gene_list %in% rownames(Hp_RNA_Seurat) & gene_list %in% rownames(Ho_RNA_Seurat) & gene_list %in% rownames(Ho_1.5_Seurat) & gene_list %in% rownames(Ho_2_Seurat) & gene_list %in% rownames(Ho_3_Seurat) & gene_list %in% rownames(Ho_4_Seurat) & gene_list %in% rownames(Paola_3_Seurat) & gene_list %in% rownames(Paola_6_Seurat)]

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in data_list) {
  assign(x = i,value = AddModuleScore(object = get(i),features = list(gene_list),name = 'ECM'))
  print(paste(i,'done!',sep = ' '))
  gc()
}

#plot
temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$re_anno %in% 'RG-neu',]
  meta_data <- meta_data[,c("re_anno","ECM1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$group <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat'),"group"] <- 'organoid'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('PFC_RNA_Seurat'),"author"] <- 'Bhaduri. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'
temp[temp$dataset %in% c('Paola_3_Seurat','Paola_6_Seurat'),"author"] <- 'Velasco. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = data_list)
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Bhaduri. et al','Uzquiano. et al','Velasco. et al'))

p1 <- ggplot(data = temp[temp$group == 'primary',],mapping = aes(x = author,y = ECM1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.1,0.15)) + 
  #scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00')) + 
  theme(aspect.ratio = 7/3.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('ECM module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p2 <- ggplot(data = temp[temp$group == 'organoid',],mapping = aes(x = author,y = ECM1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.1,0.15)) + 
  #scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#0074D9','#2ECC40')) + 
  theme(aspect.ratio = 7/2.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('ECM module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_115_fig_230405/primary_organoid_ECM_signature_boxplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

t.test(x = temp[temp$group == 'primary',"ECM1"],
       y = temp[temp$group == 'organoid',"ECM1"],
       alternative = 'greater')

# ECM binding signature -----------------------------------------------------------

#load col param
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
PFC_RNA_Seurat <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_PFC_brain_preprocess_PC30_by_lyt_230326.rds')
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')
Paola_3_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_3m_Seurat.rds')
Paola_6_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_6m_Seurat.rds')
gc()

#preprocess PD human
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,npcs = 50,preprocess = TRUE)
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'Cluster',label = TRUE)

#normalize data
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in 1:length(data_list)) {
  temp <- get(data_list[i])
  temp <- NormalizeData(temp)
  assign(x = data_list[i],value = temp)
  rm(temp)
  gc()
}

meta_data <- readRDS(file = './res/step_113_fig_230326/PFC_RG_meta_data.rds')
#PFC_RNA_Seurat <- PFC_RNA_Seurat[,rownames(meta_data)]
PFC_RNA_Seurat@meta.data[rownames(meta_data),"cell.type"] <- meta_data$cell_type

#ori_anno
Greenleaf_RNA_Seurat$ori_anno <- Greenleaf_RNA_Seurat$ReAnno_celltype
PD_RNA_Seurat$ori_anno <- PD_RNA_Seurat$Cluster
PFC_RNA_Seurat$ori_anno <- PFC_RNA_Seurat$cell.type
Hp_RNA_Seurat$ori_anno <- Hp_RNA_Seurat$cell_type
Ho_RNA_Seurat$ori_anno <- Ho_RNA_Seurat$cell_type
Ho_1.5_Seurat$ori_anno <- Ho_1.5_Seurat$FinalName
Ho_1.5_Seurat$ori_anno <- as.character(Ho_1.5_Seurat$ori_anno)
Ho_2_Seurat$ori_anno <- Ho_2_Seurat$FinalName
Ho_2_Seurat$ori_anno <- as.character(Ho_2_Seurat$ori_anno)
Ho_3_Seurat$ori_anno <- Ho_3_Seurat$FinalName
Ho_3_Seurat$ori_anno <- as.character(Ho_3_Seurat$ori_anno)
Ho_4_Seurat$ori_anno <- Ho_4_Seurat$FinalName
Ho_4_Seurat$ori_anno <- as.character(Ho_4_Seurat$ori_anno)
Paola_3_Seurat$ori_anno <- Paola_3_Seurat$CellType
Paola_6_Seurat$ori_anno <- Paola_6_Seurat$CellType

#colnames
Greenleaf_RNA_Seurat <- RenameCells(object = Greenleaf_RNA_Seurat,new.names = paste('Greenleaf',colnames(Greenleaf_RNA_Seurat),sep = '_'))
PD_RNA_Seurat <- RenameCells(object = PD_RNA_Seurat,new.names = paste('PD',colnames(PD_RNA_Seurat),sep = '_'))
PFC_RNA_Seurat <- RenameCells(object = PFC_RNA_Seurat,new.names = paste('PFC',colnames(PFC_RNA_Seurat),sep = '_'))
Hp_RNA_Seurat <- RenameCells(object = Hp_RNA_Seurat,new.names = paste('Hp',colnames(Hp_RNA_Seurat),sep = '_'))
Ho_RNA_Seurat <- RenameCells(object = Ho_RNA_Seurat,new.names = paste('Ho',colnames(Ho_RNA_Seurat),sep = '_'))
Ho_1.5_Seurat <- RenameCells(object = Ho_1.5_Seurat,new.names = paste('Ho_1.5',colnames(Ho_1.5_Seurat),sep = '_'))
Ho_2_Seurat <- RenameCells(object = Ho_2_Seurat,new.names = paste('Ho_2',colnames(Ho_2_Seurat),sep = '_'))
Ho_3_Seurat <- RenameCells(object = Ho_3_Seurat,new.names = paste('Ho_3',colnames(Ho_3_Seurat),sep = '_'))
Ho_4_Seurat <- RenameCells(object = Ho_4_Seurat,new.names = paste('Ho_4',colnames(Ho_4_Seurat),sep = '_'))
Paola_3_Seurat <- RenameCells(object = Paola_3_Seurat,new.names = paste('Paola_3',colnames(Paola_3_Seurat),sep = '_'))
Paola_6_Seurat <- RenameCells(object = Paola_6_Seurat,new.names = paste('Paola_6',colnames(Paola_6_Seurat),sep = '_'))

#re-anno
Greenleaf_RNA_Seurat$re_anno <- Greenleaf_RNA_Seurat$ori_anno
PD_RNA_Seurat$re_anno <- PD_RNA_Seurat$ori_anno
PFC_RNA_Seurat$re_anno <- PFC_RNA_Seurat$ori_anno
Hp_RNA_Seurat$re_anno <- Hp_RNA_Seurat$ori_anno
Ho_RNA_Seurat$re_anno <- Ho_RNA_Seurat$ori_anno
Ho_1.5_Seurat$re_anno <- Ho_1.5_Seurat$ori_anno
Ho_2_Seurat$re_anno <- Ho_2_Seurat$ori_anno
Ho_3_Seurat$re_anno <- Ho_3_Seurat$ori_anno
Ho_4_Seurat$re_anno <- Ho_4_Seurat$ori_anno
Paola_3_Seurat$re_anno <- Paola_3_Seurat$ori_anno
Paola_6_Seurat$re_anno <- Paola_6_Seurat$ori_anno

Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1'),"re_anno"] <- 'RG-neu'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$ori_anno %in% c('oRG','vRG'),"re_anno"] <- 'RG-neu'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('RG-neu'),"re_anno"] <- 'RG-neu'
Hp_RNA_Seurat@meta.data[Hp_RNA_Seurat$cell_type == 'RG',"re_anno"] <- 'RG-neu'
Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type == 'RG',"re_anno"] <- 'RG-neu'
Ho_1.5_Seurat@meta.data[Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Ho_2_Seurat@meta.data[Ho_2_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Ho_3_Seurat@meta.data[Ho_3_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Ho_4_Seurat@meta.data[Ho_4_Seurat$ori_anno %in% c('aRG','oRG'),"re_anno"] <- 'RG-neu'
Paola_3_Seurat@meta.data[Paola_3_Seurat$ori_anno %in% c('oRG','RG'),"re_anno"] <- 'RG-neu'
Paola_6_Seurat@meta.data[Paola_6_Seurat$ori_anno %in% c('oRG','RG'),"re_anno"] <- 'RG-neu'

#get gene list
gene_list <- getGOgeneSet(x = 'GO:0050840',OrgDb = org.Hs.eg.db,ont = 'MF',keytype = 'SYMBOL')
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat) & gene_list %in% rownames(PD_RNA_Seurat) & gene_list %in% rownames(PFC_RNA_Seurat) & gene_list %in% rownames(Hp_RNA_Seurat) & gene_list %in% rownames(Ho_RNA_Seurat) & gene_list %in% rownames(Ho_1.5_Seurat) & gene_list %in% rownames(Ho_2_Seurat) & gene_list %in% rownames(Ho_3_Seurat) & gene_list %in% rownames(Ho_4_Seurat) & gene_list %in% rownames(Paola_3_Seurat) & gene_list %in% rownames(Paola_6_Seurat)]

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in data_list) {
  assign(x = i,value = AddModuleScore(object = get(i),features = list(gene_list),name = 'ECM'))
  print(paste(i,'done!',sep = ' '))
  gc()
}

#plot
temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$re_anno %in% 'RG-neu',]
  meta_data <- meta_data[,c("re_anno","ECM1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$group <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat'),"group"] <- 'organoid'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('PFC_RNA_Seurat'),"author"] <- 'Bhaduri. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'
temp[temp$dataset %in% c('Paola_3_Seurat','Paola_6_Seurat'),"author"] <- 'Velasco. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = data_list)
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Bhaduri. et al','Uzquiano. et al','Velasco. et al'))

p1 <- ggplot(data = temp[temp$group == 'primary',],mapping = aes(x = author,y = ECM1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.1,0.3)) + 
  #scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00')) + 
  theme(aspect.ratio = 7/3.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('ECM binding module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p2 <- ggplot(data = temp[temp$group == 'organoid',],mapping = aes(x = author,y = ECM1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.1,0.3)) + 
  #scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#0074D9','#2ECC40')) + 
  theme(aspect.ratio = 7/2.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('ECM binding module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_115_fig_230405/primary_organoid_ECM_binding_signature_boxplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

t.test(x = temp[temp$group == 'primary',"ECM1"],
       y = temp[temp$group == 'organoid',"ECM1"],
       alternative = 'greater')

# output bar plot in supplementary -----------------------------------------

#load col param
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
PFC_RNA_Seurat <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_PFC_brain_preprocess_PC30_by_lyt_230326.rds')
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')
Paola_3_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_3m_Seurat.rds')
Paola_6_Seurat <- readRDS(file = './data/public/Individual_brain_organoids_reproducibly_form_cell_diversity_of_the_human_cerebral_cortex/processed_data/Organoid_6m_Seurat.rds')
gc()

#preprocess PD human
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,npcs = 50,preprocess = TRUE)
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'Cluster',label = TRUE)

#normalize data
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in 1:length(data_list)) {
  temp <- get(data_list[i])
  temp <- NormalizeData(temp)
  assign(x = data_list[i],value = temp)
  rm(temp)
  gc()
}

meta_data <- readRDS(file = './res/step_113_fig_230326/PFC_RG_meta_data.rds')
#PFC_RNA_Seurat <- PFC_RNA_Seurat[,rownames(meta_data)]
PFC_RNA_Seurat@meta.data[rownames(meta_data),"cell.type"] <- meta_data$cell_type

#ori_anno
Greenleaf_RNA_Seurat$ori_anno <- Greenleaf_RNA_Seurat$ReAnno_celltype
PD_RNA_Seurat$ori_anno <- PD_RNA_Seurat$Cluster
PFC_RNA_Seurat$ori_anno <- PFC_RNA_Seurat$cell.type
Hp_RNA_Seurat$ori_anno <- Hp_RNA_Seurat$cell_type
Ho_RNA_Seurat$ori_anno <- Ho_RNA_Seurat$cell_type
Ho_1.5_Seurat$ori_anno <- Ho_1.5_Seurat$FinalName
Ho_1.5_Seurat$ori_anno <- as.character(Ho_1.5_Seurat$ori_anno)
Ho_2_Seurat$ori_anno <- Ho_2_Seurat$FinalName
Ho_2_Seurat$ori_anno <- as.character(Ho_2_Seurat$ori_anno)
Ho_3_Seurat$ori_anno <- Ho_3_Seurat$FinalName
Ho_3_Seurat$ori_anno <- as.character(Ho_3_Seurat$ori_anno)
Ho_4_Seurat$ori_anno <- Ho_4_Seurat$FinalName
Ho_4_Seurat$ori_anno <- as.character(Ho_4_Seurat$ori_anno)
Paola_3_Seurat$ori_anno <- Paola_3_Seurat$CellType
Paola_6_Seurat$ori_anno <- Paola_6_Seurat$CellType

#colnames
Greenleaf_RNA_Seurat <- RenameCells(object = Greenleaf_RNA_Seurat,new.names = paste('Greenleaf',colnames(Greenleaf_RNA_Seurat),sep = '_'))
PD_RNA_Seurat <- RenameCells(object = PD_RNA_Seurat,new.names = paste('PD',colnames(PD_RNA_Seurat),sep = '_'))
PFC_RNA_Seurat <- RenameCells(object = PFC_RNA_Seurat,new.names = paste('PFC',colnames(PFC_RNA_Seurat),sep = '_'))
Hp_RNA_Seurat <- RenameCells(object = Hp_RNA_Seurat,new.names = paste('Hp',colnames(Hp_RNA_Seurat),sep = '_'))
Ho_RNA_Seurat <- RenameCells(object = Ho_RNA_Seurat,new.names = paste('Ho',colnames(Ho_RNA_Seurat),sep = '_'))
Ho_1.5_Seurat <- RenameCells(object = Ho_1.5_Seurat,new.names = paste('Ho_1.5',colnames(Ho_1.5_Seurat),sep = '_'))
Ho_2_Seurat <- RenameCells(object = Ho_2_Seurat,new.names = paste('Ho_2',colnames(Ho_2_Seurat),sep = '_'))
Ho_3_Seurat <- RenameCells(object = Ho_3_Seurat,new.names = paste('Ho_3',colnames(Ho_3_Seurat),sep = '_'))
Ho_4_Seurat <- RenameCells(object = Ho_4_Seurat,new.names = paste('Ho_4',colnames(Ho_4_Seurat),sep = '_'))
Paola_3_Seurat <- RenameCells(object = Paola_3_Seurat,new.names = paste('Paola_3',colnames(Paola_3_Seurat),sep = '_'))
Paola_6_Seurat <- RenameCells(object = Paola_6_Seurat,new.names = paste('Paola_6',colnames(Paola_6_Seurat),sep = '_'))

#filter
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1')]
PD_RNA_Seurat <- PD_RNA_Seurat[,PD_RNA_Seurat$ori_anno %in% c('oRG','vRG')]
PFC_RNA_Seurat <- PFC_RNA_Seurat[,PFC_RNA_Seurat$ori_anno %in% c('RG-neu')]
Hp_RNA_Seurat <- Hp_RNA_Seurat[,Hp_RNA_Seurat$cell_type == 'RG']
Ho_RNA_Seurat <- Ho_RNA_Seurat[,Ho_RNA_Seurat$cell_type == 'RG']
Ho_1.5_Seurat <- Ho_1.5_Seurat[,Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG')]
Ho_2_Seurat <- Ho_2_Seurat[,Ho_2_Seurat$ori_anno %in% c('aRG','oRG')]
Ho_3_Seurat <- Ho_3_Seurat[,Ho_3_Seurat$ori_anno %in% c('aRG','oRG')]
Ho_4_Seurat <- Ho_4_Seurat[,Ho_4_Seurat$ori_anno %in% c('aRG','oRG')]
Paola_3_Seurat <- Paola_3_Seurat[,Paola_3_Seurat$ori_anno %in% c('oRG','RG')]
Paola_6_Seurat <- Paola_6_Seurat[,Paola_6_Seurat$ori_anno %in% c('oRG','RG')]
gc()

gene_list <- Reduce(f = intersect,x = list(rownames(Greenleaf_RNA_Seurat),
                                           rownames(PD_RNA_Seurat),
                                           rownames(PFC_RNA_Seurat),
                                           rownames(Hp_RNA_Seurat),
                                           rownames(Ho_RNA_Seurat),
                                           rownames(Ho_1.5_Seurat),
                                           rownames(Paola_3_Seurat)))

RG_RNA_Seurat <- cbind(expm1(Greenleaf_RNA_Seurat@assays$RNA@data[gene_list,Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1')]),
                       expm1(PD_RNA_Seurat@assays$RNA@data[gene_list,PD_RNA_Seurat$ori_anno %in% c('oRG','vRG')]),
                       expm1(PFC_RNA_Seurat@assays$RNA@data[gene_list,PFC_RNA_Seurat$ori_anno %in% c('RG-neu')]),
                       expm1(Hp_RNA_Seurat@assays$RNA@data[gene_list,Hp_RNA_Seurat$cell_type == 'RG']),
                       expm1(Ho_RNA_Seurat@assays$RNA@data[gene_list,Ho_RNA_Seurat$cell_type == 'RG']),
                       expm1(Ho_1.5_Seurat@assays$RNA@data[gene_list,Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Ho_2_Seurat@assays$RNA@data[gene_list,Ho_2_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Ho_3_Seurat@assays$RNA@data[gene_list,Ho_3_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Ho_4_Seurat@assays$RNA@data[gene_list,Ho_4_Seurat$ori_anno %in% c('aRG','oRG')]),
                       expm1(Paola_3_Seurat@assays$RNA@data[gene_list,Paola_3_Seurat$ori_anno %in% c('oRG','RG')]),
                       expm1(Paola_6_Seurat@assays$RNA@data[gene_list,Paola_6_Seurat$ori_anno %in% c('oRG','RG')]))

meta_data <- data.frame(cell_name = colnames(RG_RNA_Seurat),
                        dataset = c(rep('Trevino. et al.',times = sum(Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1'))),
                                    rep('Polioudakis. et al.',times = sum(PD_RNA_Seurat$ori_anno %in% c('oRG','vRG'))),
                                    rep('Bhaduri. et al.',times = sum(PFC_RNA_Seurat$ori_anno %in% c('RG-neu'))),
                                    rep('Pollen. et al. Primary',times = sum(Hp_RNA_Seurat$cell_type == 'RG')),
                                    rep('Pollen. et al. Organoid',times = sum(Ho_RNA_Seurat$cell_type == 'RG')),
                                    rep('Uzquiano. et al.',times = sum(c(Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG'),Ho_2_Seurat$ori_anno %in% c('aRG','oRG'),Ho_3_Seurat$ori_anno %in% c('aRG','oRG'),Ho_4_Seurat$ori_anno %in% c('aRG','oRG')))),
                                    rep('Velasco. et al.',times = sum(c(Paola_3_Seurat$ori_anno %in% c('oRG','RG'),Paola_6_Seurat$ori_anno %in% c('oRG','RG'))))))
rownames(meta_data) <- meta_data$cell_name

RG_RNA_Seurat <- CreateSeuratObject(counts = RG_RNA_Seurat,project = 'RG',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

#add meta data
RG_RNA_Seurat$Age <- ''

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Greenleaf_RNA_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- Greenleaf_RNA_Seurat@meta.data[cell_list,"Age"]

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(PD_RNA_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- PD_RNA_Seurat@meta.data[cell_list,"Gestation_week"]

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(PFC_RNA_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- PFC_RNA_Seurat@meta.data[cell_list,"age"]

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_1.5_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '1.5m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_2_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '2m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_3_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '3m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Ho_4_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '4m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Paola_3_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '3m'

cell_list <- colnames(RG_RNA_Seurat)[colnames(RG_RNA_Seurat) %in% colnames(Paola_6_Seurat)]
RG_RNA_Seurat@meta.data[cell_list,"Age"] <- '6m'

RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '17')] <- '17GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '18')] <- '18GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '19')] <- '19GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '20')] <- '20GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '22')] <- '22GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == '25')] <- '25GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == 'pcw16')] <- '16GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == 'pcw20')] <- '20GW'
RG_RNA_Seurat$Age[which(RG_RNA_Seurat$Age == 'pcw21')] <- '21GW'

RG_RNA_Seurat$sample <- paste(RG_RNA_Seurat$dataset,RG_RNA_Seurat$Age,sep = ' ')
RG_RNA_Seurat@meta.data[RG_RNA_Seurat$sample == 'Pollen. et al. Primary ',"sample"] <- 'Pollen. et al. Primary'
RG_RNA_Seurat@meta.data[RG_RNA_Seurat$sample == 'Pollen. et al. Organoid ',"sample"] <- 'Pollen. et al. Organoid'

RG_RNA_Seurat$group <- 'Primary'
RG_RNA_Seurat@meta.data[RG_RNA_Seurat$dataset %in% c('Uzquiano. et al.','Velasco. et al.'),"group"] <- 'Organoid'
RG_RNA_Seurat@meta.data[colnames(RG_RNA_Seurat) %in% colnames(Ho_RNA_Seurat),"group"] <- 'Organoid'

#normalize
RG_RNA_Seurat <- NormalizeData(RG_RNA_Seurat)
gc()

#factor
RG_RNA_Seurat$dataset <- factor(RG_RNA_Seurat$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al. Primary','Pollen. et al. Organoid','Uzquiano. et al.','Velasco. et al.'))
RG_RNA_Seurat$sample <- factor(RG_RNA_Seurat$sample,levels = c('Trevino. et al. 16GW','Trevino. et al. 20GW','Trevino. et al. 21GW',
                                                               'Polioudakis. et al. 17GW','Polioudakis. et al. 18GW',
                                                               'Bhaduri. et al. 17GW','Bhaduri. et al. 18GW','Bhaduri. et al. 19GW','Bhaduri. et al. 20GW','Bhaduri. et al. 22GW','Bhaduri. et al. 25GW',
                                                               'Pollen. et al. Primary','Pollen. et al. Organoid',
                                                               'Uzquiano. et al. 1.5m','Uzquiano. et al. 2m','Uzquiano. et al. 3m','Uzquiano. et al. 4m',
                                                               'Velasco. et al. 3m','Velasco. et al. 6m'))
RG_RNA_Seurat$group <- factor(RG_RNA_Seurat$group,levels = c('Primary','Organoid'))

#add meta data
RG_RNA_Seurat$donor <- ''
RG_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"donor"] <- Greenleaf_RNA_Seurat$Sample.ID
RG_RNA_Seurat@meta.data[colnames(PD_RNA_Seurat),"donor"] <- PD_RNA_Seurat$Donor
RG_RNA_Seurat@meta.data[colnames(PFC_RNA_Seurat),"donor"] <- PFC_RNA_Seurat$individual
RG_RNA_Seurat@meta.data[colnames(Hp_RNA_Seurat),"donor"] <- Hp_RNA_Seurat$Age
RG_RNA_Seurat@meta.data[colnames(Ho_RNA_Seurat),"donor"] <- Ho_RNA_Seurat$Age
RG_RNA_Seurat@meta.data[colnames(Ho_1.5_Seurat),"donor"] <- Ho_1.5_Seurat$dataset
RG_RNA_Seurat@meta.data[colnames(Ho_2_Seurat),"donor"] <- Ho_2_Seurat$dataset
RG_RNA_Seurat@meta.data[colnames(Ho_3_Seurat),"donor"] <- Ho_3_Seurat$dataset
RG_RNA_Seurat@meta.data[colnames(Ho_4_Seurat),"donor"] <- Ho_4_Seurat$dataset
RG_RNA_Seurat@meta.data[colnames(Paola_3_Seurat),"donor"] <- Paola_3_Seurat$Organoid
RG_RNA_Seurat@meta.data[colnames(Paola_6_Seurat),"donor"] <- Paola_6_Seurat$Organoid

RG_RNA_Seurat$sample_2 <- paste(RG_RNA_Seurat$sample,RG_RNA_Seurat$donor,sep = '_')

#filter data
filter_list <- names(table(RG_RNA_Seurat$sample_2))[table(RG_RNA_Seurat$sample_2) < 10 & !grepl(pattern = '^Pollen',x = names(table(RG_RNA_Seurat$sample_2)),fixed = FALSE)]
RG_RNA_Seurat <- RG_RNA_Seurat[,!(RG_RNA_Seurat$sample_2 %in% filter_list)]

#normalize data
RG_RNA_Seurat <- NormalizeData(RG_RNA_Seurat)

#calculate percentage
gene_list <- c('ANXA2','ITGA2','LGALS3','LRIG1','ITGB5','BMP7','COL11A1','F3','COL4A5','SFRP1','MFGE8','GPC4','TNC')
plot_data <- DotPlot(object = RG_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'sample_2')$data
plot_data$dataset <- gsub(pattern = ' et al.*',replacement = '',x = plot_data$id,fixed = FALSE)
plot_data$dataset <- paste(plot_data$dataset,'et al.',sep = ' ')
plot_data$dataset <- factor(plot_data$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
plot_data$group <- 'Organoid'
plot_data[grep(pattern = 'Trevino|Polioudakis|Bhaduri|Pollen. et al. Primary',x = plot_data$id,fixed = FALSE),"group"] <- 'Primary'

for (i in gene_list) {
  pct_matrix <- plot_data[plot_data$features.plot == i & plot_data$dataset != 'Pollen. et al.' & !grepl(pattern = '^Velasco. et al. 6m',x = plot_data$id,fixed = FALSE),]
  pct_matrix %>% dplyr::summarize(mean_pct = mean(avg.exp),sd_pct = sd(avg.exp),.by = dataset) -> df_bar
  df_bar$dataset <- factor(df_bar$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
  p_value <- t.test(x = pct_matrix$avg.exp[pct_matrix$group == 'Primary'],y = pct_matrix$avg.exp[pct_matrix$group == 'Organoid'])$p.value
  p_value <- format(p_value,scientific = TRUE,digits = 3)
  df_bar$group = 'Organoid'
  df_bar[grep(pattern = 'Trevino|Polioudakis|Bhaduri',x = df_bar$dataset,fixed = FALSE),'group'] <- 'Primary'
  df_bar$group <- factor(df_bar$group,levels = c('Primary','Organoid'))
  
  p <- ggplot(data = df_bar,mapping = aes(x = dataset,y = mean_pct,fill = dataset,pattern = group)) + 
    geom_bar_pattern(stat = 'identity',width = 0.6,pattern_color = 'white',pattern_fill = 'white',pattern_density = 0.1) + 
    geom_text(x = 3.5,y = max(df_bar$mean_pct),label = paste('p =',p_value,sep = ' '),size = 4,family = 'sans',fontface = 'italic') + 
    scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00','#0074D9','#2ECC40')) + 
    scale_pattern_discrete(choices = c('none','stripe')) + 
    geom_errorbar(mapping = aes(ymin = mean_pct - 0.5*sd_pct,ymax = mean_pct + 0.5*sd_pct),width = 0.2,color = 'black',linewidth = 0.2) + 
    theme_ArchR() + theme(aspect.ratio = 2) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 8)) + 
    labs(title = i) + CenterTitle() + 
    ylab('Average expression')
  
  char <- paste0('./res/step_115_fig_230405/bar_plot_output/primary_organoid_',i,'_expression_level_barplot.pdf')
  pdf(file = char,width = 5,height = 5)
  print(p)
  dev.off()
}
