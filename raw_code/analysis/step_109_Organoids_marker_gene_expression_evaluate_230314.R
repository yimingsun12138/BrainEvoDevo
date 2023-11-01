#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Organoids marker gene expression evaluate                       ##
## Data: 2023.03.14                                                                ##
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
#Greenleaf data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

#human primary
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Hp_RNA_Seurat <- Hp_RNA_Seurat[,Hp_RNA_Seurat$Age >= 17]

#macaque multiome data
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

#mouse multiome data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#macaque primary data
Mp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Mp_RNA_Seurat_230302.rds')

#human organoid
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')

Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')

meta_data <- do.call(what = rbind,args = base::lapply(X = c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),FUN = function(x){
  temp <- get(x)@meta.data[,c('dataset','percent.mito','percent.ribo','S.Score','G2M.Score','Phase','CC.Difference','FinalName')]
  rownames(temp) <- paste(x,rownames(temp),sep = '_')
  return(temp)
}))
express_matrix <- do.call(what = cbind,args = base::lapply(X = c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),FUN = function(x){
  temp <- get(x)@assays$RNA@counts
  colnames(temp) <- paste(x,colnames(temp),sep = '_')
  return(temp)
}))
Ho_Paola_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'Organoid',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

rm(express_matrix)
gc()

Ho_Paola_Seurat <- NormalizeData(Ho_Paola_Seurat)

#Chimp organoid
Co_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Co_RNA_Seurat_230302.rds')

#human primary RNA Seurat from PD
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')

#parameter
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230211.rds')

# re anno data and col parameter ------------------------------------------
#Greenleaf data
Greenleaf_RNA_Seurat$cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype
Greenleaf_RNA_Seurat$ori_cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('RG-1'),"cell_type"] <- 'RG-neu'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"cell_type"] <- 'Ex'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('InCGE','InMGE'),"cell_type"] <- 'In'
Greenleaf_RNA_Seurat@meta.data[!(Greenleaf_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

#macaque multiome
macaque_multiome_Seurat$ori_cell_type <- macaque_multiome_Seurat$cell_type
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$ori_cell_type %in% c('RG-1'),"cell_type"] <- 'RG-neu'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$ori_cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"cell_type"] <- 'Ex'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$ori_cell_type %in% c('InCGE','InMGE'),"cell_type"] <- 'In'
macaque_multiome_Seurat@meta.data[!(macaque_multiome_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'
table(macaque_multiome_Seurat$cell_type)

#mouse multiome
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
mouse_multiome_Seurat$ori_cell_type <- meta_data[colnames(mouse_multiome_Seurat),"ReAnno_celltype"]
mouse_multiome_Seurat$cell_type <- mouse_multiome_Seurat$ori_cell_type

mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$ori_cell_type %in% c('Cyc-IP','Cyc-RG'),"cell_type"] <- 'Cycling'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$ori_cell_type %in% c('RG-1'),"cell_type"] <- 'RG-neu'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$ori_cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"cell_type"] <- 'Ex'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$ori_cell_type %in% c('InMGE','InCGE'),"cell_type"] <- 'In'
mouse_multiome_Seurat@meta.data[!(mouse_multiome_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'
table(mouse_multiome_Seurat$cell_type)

#Human organoid
Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type %in% c('RG'),"cell_type"] <- 'RG-neu'
Ho_RNA_Seurat@meta.data[!(Ho_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

#human  Paola organoid
Ho_Paola_Seurat$cell_type <- Ho_Paola_Seurat$FinalName
table(Ho_Paola_Seurat$cell_type)

Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('aRG','oRG'),"cell_type"] <- 'RG-neu'
Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('CFuPN','CPN','Newborn CFuPN','Newborn CPN','Newborn DL PN','Newborn PN','PN','Subcortical neurons','Subcortical progenitors'),"cell_type"] <- 'Ex'
Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('Immature IN'),"cell_type"] <- 'In'
Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('IP'),"cell_type"] <- 'IP'
Ho_Paola_Seurat@meta.data[!(Ho_Paola_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

table(Ho_Paola_Seurat$cell_type)

Ho_1.5_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_1.5_Seurat',colnames(Ho_1.5_Seurat),sep = '_'),"cell_type"]
Ho_2_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_2_Seurat',colnames(Ho_2_Seurat),sep = '_'),"cell_type"]
Ho_3_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_3_Seurat',colnames(Ho_3_Seurat),sep = '_'),"cell_type"]
Ho_4_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_4_Seurat',colnames(Ho_4_Seurat),sep = '_'),"cell_type"]

#human primary
table(Hp_RNA_Seurat$cell_type)

Hp_RNA_Seurat@meta.data[Hp_RNA_Seurat$cell_type %in% c('RG'),"cell_type"] <- 'RG-neu'
Hp_RNA_Seurat@meta.data[!(Hp_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

#PD human
PD_RNA_Seurat$cell_type <- 'others'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('oRG','vRG'),'cell_type'] <- 'RG-neu'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('vRG'),'cell_type'] <- 'vRG'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('PgG2M','PgS'),'cell_type'] <- 'Cycling'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('IP'),'cell_type'] <- 'IP'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('ExDp1','ExDp2','ExM','ExM-U','ExN'),'cell_type'] <- 'Ex'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('InCGE','InMGE'),'cell_type'] <- 'In'

# subset data -------------------------------------------------------------
cell_type_list <- c('RG-neu','vRG','IP','Ex','In','others')
for (i in c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Ho_Paola_Seurat')) {
  temp <- get(x = i)
  if(i != 'Ho_Paola_Seurat'){
    temp <- NormalizeData(temp)
  }
  temp <- temp[,temp$cell_type %in% cell_type_list]
  temp$cell_type <- factor(temp$cell_type,levels = cell_type_list)
  Idents(temp) <- 'cell_type'
  assign(x = i,value = temp)
  gc()
}

# human specific RG-neu signature gene  ---------------------------------------------------------
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
gene_list <- HG_gene_list
scale.min <- 0
scale.max <- base::lapply(X = list(Greenleaf_RNA_Seurat,PD_RNA_Seurat,Hp_RNA_Seurat,Ho_RNA_Seurat,Ho_1.5_Seurat,Ho_2_Seurat,Ho_3_Seurat,Ho_4_Seurat),FUN = function(x){
  temp <- DotPlot(object = x,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE)
  return(max(temp$data$pct.exp))
})
scale.max <- max(unlist(scale.max))

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  p <- DotPlot(object = get(data_list[i]),assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
               cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max,
               idents = c('RG-neu','vRG','IP','Ex','In')) + 
    coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    theme(aspect.ratio = length(gene_list)/length(unique(get(data_list[i])$cell_type)),
          text = element_text(size = 10,family = 'sans')) + 
    labs(title = data_list[i])
  if(i != 1){
    p <- p + theme(axis.text.y = element_blank())
  }
  if(i != length(data_list)){
    p <- p + NoLegend()
  }
  assign(x = paste0('p',i),value = p)
}

pdf(file = './res/step_109_fig_230314/human_primary_organoid_HG_RG_neu_signature_gene_dotplot.pdf',width = 18,height = 16)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 8)
dev.off()

# primate specific RG-neu signature gene  ---------------------------------------------------------
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]
gene_list <- PC_gene_list
scale.min <- 0
scale.max <- base::lapply(X = list(Greenleaf_RNA_Seurat,PD_RNA_Seurat,Hp_RNA_Seurat,Ho_RNA_Seurat,Ho_1.5_Seurat,Ho_2_Seurat,Ho_3_Seurat,Ho_4_Seurat),FUN = function(x){
  temp <- DotPlot(object = x,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE)
  return(max(temp$data$pct.exp))
})
scale.max <- max(unlist(scale.max))

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  p <- DotPlot(object = get(data_list[i]),assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
               cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max,
               idents = c('RG-neu','vRG','IP','Ex','In')) + 
    coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    theme(aspect.ratio = length(gene_list)/length(unique(get(data_list[i])$cell_type)),
          text = element_text(size = 10,family = 'sans')) + 
    labs(title = data_list[i])
  if(i != 1){
    p <- p + theme(axis.text.y = element_blank())
  }
  if(i != length(data_list)){
    p <- p + NoLegend()
  }
  assign(x = paste0('p',i),value = p)
}

pdf(file = './res/step_109_fig_230314/human_primary_organoid_PC_RG_neu_signature_gene_dotplot.pdf',width = 18,height = 16)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 8)
dev.off()

# species conserved RG-neu signature gene  ---------------------------------------------------------
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
gene_list <- SC_gene_list
scale.min <- 0
scale.max <- base::lapply(X = list(Greenleaf_RNA_Seurat,PD_RNA_Seurat,Hp_RNA_Seurat,Ho_RNA_Seurat,Ho_1.5_Seurat,Ho_2_Seurat,Ho_3_Seurat,Ho_4_Seurat),FUN = function(x){
  temp <- DotPlot(object = x,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE)
  return(max(temp$data$pct.exp))
})
scale.max <- max(unlist(scale.max))

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  p <- DotPlot(object = get(data_list[i]),assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
               cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max,
               idents = c('RG-neu','vRG','IP','Ex','In')) + 
    coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    theme(aspect.ratio = length(gene_list)/length(unique(get(data_list[i])$cell_type)),
          text = element_text(size = 10,family = 'sans')) + 
    labs(title = data_list[i])
  if(i != 1){
    p <- p + theme(axis.text.y = element_blank())
  }
  if(i != length(data_list)){
    p <- p + NoLegend()
  }
  assign(x = paste0('p',i),value = p)
}

pdf(file = './res/step_109_fig_230314/human_primary_organoid_SC_RG_neu_signature_gene_dotplot.pdf',width = 18,height = 16)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 8)
dev.off()

# # featureplot -------------------------------------------------------------
# Ho_1.5_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_1.5_Seurat',colnames(Ho_1.5_Seurat),sep = '_'),"cell_type"]
# Ho_2_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_2_Seurat',colnames(Ho_2_Seurat),sep = '_'),"cell_type"]
# Ho_3_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_3_Seurat',colnames(Ho_3_Seurat),sep = '_'),"cell_type"]
# Ho_4_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_4_Seurat',colnames(Ho_4_Seurat),sep = '_'),"cell_type"]
# 
# col_list <- col_param$celltype[c("RG","Cycling",'IP',"Ex-1","InMGE","Mic")]
# names(col_list) <- c('RG-neu','Cycling','IP','Ex','In','others')
# 
# # most important ECM marker -----------------------------------------------
# data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
# for (i in 1:length(data_list)) {
#   p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'cell_type',label = TRUE,repel = TRUE,raster = FALSE,cols = col_list) + 
#     theme_bw() + theme(aspect.ratio = 1) + NoLegend() + CenterTitle() + 
#     labs(title = data_list[i]) + 
#     theme(axis.text = element_blank(),
#           axis.title = element_blank(),
#           axis.ticks = element_blank())
#   assign(x = paste0('p',i),value = p)
#   
#   #feature plot
#   gene_list <- c('ANXA2','ITGA2','LGALS3','ITGB5','COL11A1')
#   for (j in 1:length(gene_list)) {
#     
#     indi <- TRUE
#     p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],pt.size = 0.001,order = indi,raster = FALSE,cols = ArchRPalettes$whitePurple[-1]) + 
#       theme_bw() + theme(aspect.ratio = 1) + NoLegend() + CenterTitle() + 
#       theme(axis.text = element_blank(),
#             axis.title = element_blank(),
#             axis.ticks = element_blank(),
#             plot.title = element_blank())
#     assign(x = paste0('p',i+j*length(data_list)),value = p)
#   }
# }
# 
# char <- paste0('p',1:(8*6)) %>% paste(collapse = '+') %>% paste0('+plot_layout(ncol = 8)')
# 
# png(filename = './res/step_109_fig_230314/human_primary_organoid_ECM_gene_featureplot.png',width = 2400,height = 1800)
# eval(parse(text = char))
# dev.off()

# human primary organoid expression pct change ----------------------------
#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

#calculate percentage
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_1.5 <- DotPlot(object = Ho_1.5_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_2 <- DotPlot(object = Ho_2_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_3 <- DotPlot(object = Ho_3_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_4 <- DotPlot(object = Ho_4_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"pct.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"pct.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"pct.exp"],p_Paola_1.5$data[p_Paola_1.5$data$id == 'RG-neu' & p_Paola_1.5$data$features.plot == y,"pct.exp"],p_Paola_2$data[p_Paola_2$data$id == 'RG-neu' & p_Paola_2$data$features.plot == y,"pct.exp"],p_Paola_3$data[p_Paola_3$data$id == 'RG-neu' & p_Paola_3$data$features.plot == y,"pct.exp"],p_Paola_4$data[p_Paola_4$data$id == 'RG-neu' & p_Paola_4$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = rep(x = get(x),times = 8),gene_group = x,pct = c(primary,organoid),
                     group = c(rep(x = 'primary',times = length(primary)),rep(x = 'organoid',times = length(organoid))),
                     sample = c(rep('Trevino, et al',times = length(get(x))),rep('Polioudakis, et al',times = length(get(x))),rep('Pollen, et al',times = length(get(x))),
                                rep('Pollen, et al',times = length(get(x))),rep('Uzquiano, et al, 1.5m',times = length(get(x))),rep('Uzquiano, et al, 2m',times = length(get(x))),
                                rep('Uzquiano, et al, 3m',times = length(get(x))),rep('Uzquiano, et al, 4m',times = length(get(x)))))
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))
pct_matrix$sample <- paste(pct_matrix$group,pct_matrix$sample,sep = '#')
pct_matrix$sample <- factor(pct_matrix$sample,levels = c('primary#Trevino, et al','primary#Polioudakis, et al','primary#Pollen, et al',
                                                         'organoid#Pollen, et al','organoid#Uzquiano, et al, 1.5m','organoid#Uzquiano, et al, 2m',
                                                         'organoid#Uzquiano, et al, 3m','organoid#Uzquiano, et al, 4m'))

#plot
pdf(file = './res/step_113_fig_230326/why_give_up/expression_percent_boxplot_split_dataset.pdf',width = 14,height = 6)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = pct,fill = sample)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  theme_ArchR() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 10)) + 
  xlab('expression percent')
dev.off()
# human primary organoid normalized expression change ---------------------

## split Paola data --------------------------------------------------------
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  temp <- get(data_list[i])
  temp <- NormalizeData(temp)
  assign(x = data_list[i],value = temp)
  gc()
}

#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type')
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type')
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type')
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type')
  p_Paola_1.5 <- DotPlot(object = Ho_1.5_Seurat,features = get(x),group.by = 'cell_type')
  p_Paola_2 <- DotPlot(object = Ho_2_Seurat,features = get(x),group.by = 'cell_type')
  p_Paola_3 <- DotPlot(object = Ho_3_Seurat,features = get(x),group.by = 'cell_type')
  p_Paola_4 <- DotPlot(object = Ho_4_Seurat,features = get(x),group.by = 'cell_type')
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp"],p_Paola_1.5$data[p_Paola_1.5$data$id == 'RG-neu' & p_Paola_1.5$data$features.plot == y,"avg.exp"],p_Paola_2$data[p_Paola_2$data$id == 'RG-neu' & p_Paola_2$data$features.plot == y,"avg.exp"],p_Paola_3$data[p_Paola_3$data$id == 'RG-neu' & p_Paola_3$data$features.plot == y,"avg.exp"],p_Paola_4$data[p_Paola_4$data$id == 'RG-neu' & p_Paola_4$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = rep(x = get(x),times = 8),
                     gene_group = x,express = c(primary,organoid),
                     group = c(rep(x = 'primary',times = length(primary)),
                               rep(x = 'organoid',times = length(organoid))),
                     sample = c(rep(x = 'Trevino, et al',times = length(get(x))),rep(x = 'Polioudakis, et al',times = length(get(x))),rep(x = 'Pollen, et al',times = length(get(x))),rep(x = 'Pollen, et al',times = length(get(x))),
                                rep(x = 'Uzquiano, et al, 1.5m',times = length(get(x))),rep(x = 'Uzquiano, et al, 2m',times = length(get(x))),rep(x = 'Uzquiano, et al, 3m',times = length(get(x))),rep(x = 'Uzquiano, et al, 4m',times = length(get(x)))))
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))
pct_matrix$sample <- factor(pct_matrix$sample,
                            levels = c('Trevino, et al','Polioudakis, et al',
                                       'Pollen, et al',
                                       'Uzquiano, et al, 1.5m','Uzquiano, et al, 2m',
                                       'Uzquiano, et al, 3m','Uzquiano, et al, 4m'))

#plot
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = express,fill = sample,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  theme_ArchR() + 
  ylim(c(0,10)) + 
  scale_fill_manual(values = c('#F7DC6F','#3498DB','#2ECC71','#E74C3C','#9B59B6','#F39C12','#1ABC9C')) + 
  scale_linetype_manual(values = c('solid','dashed')) + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test')

## do not split ------------------------------------------------------------
#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = rep(x = get(x),times = 5),
                     gene_group = x,express = c(primary,organoid),
                     group = c(rep(x = 'primary',times = length(primary)),
                               rep(x = 'organoid',times = length(organoid))),
                     sample = c(rep(x = 'Trevino, et al',times = length(get(x))),rep(x = 'Polioudakis, et al',times = length(get(x))),rep(x = 'Pollen, et al',times = length(get(x))),
                                rep(x = 'Pollen, et al',times = length(get(x))),rep(x = 'Uzquiano, et al',times = length(get(x)))))
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))
pct_matrix$sample <- factor(pct_matrix$sample,
                            levels = c('Trevino, et al','Polioudakis, et al',
                                       'Pollen, et al','Uzquiano, et al'))

#plot
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = log(express),fill = sample,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  theme_ArchR() + 
  ylim(c(-3,3)) + 
  scale_fill_manual(values = c('#F7DC6F','#3498DB','#2ECC71','#E74C3C','#9B59B6','#F39C12','#1ABC9C')) + 
  scale_linetype_manual(values = c('solid','dashed')) + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test')

## do not use pollen primary ------------------------------------------------------------
#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = rep(x = get(x),times = 4),
                     gene_group = x,express = c(primary,organoid),
                     group = c(rep(x = 'primary',times = length(primary)),
                               rep(x = 'organoid',times = length(organoid))),
                     sample = c(rep(x = 'Trevino, et al',times = length(get(x))),rep(x = 'Polioudakis, et al',times = length(get(x))),
                                rep(x = 'Pollen, et al',times = length(get(x))),rep(x = 'Uzquiano, et al',times = length(get(x)))))
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))
pct_matrix$sample <- factor(pct_matrix$sample,
                            levels = c('Trevino, et al','Polioudakis, et al',
                                       'Pollen, et al','Uzquiano, et al'))

#plot
pdf(file = './res/step_109_fig_230314/human_primary_organoid_RG_neu_signature_expression_boxplot.pdf',width = 8,height = 4)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = express,fill = sample,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  theme_ArchR() + 
  ylim(c(0,4)) + 
  scale_fill_manual(values = c('#F7DC6F','#3498DB','#2ECC71','#E74C3C','#9B59B6','#F39C12','#1ABC9C')) + 
  scale_linetype_manual(values = c('solid','dashed')) + 
  scale_x_discrete(labels = c('Human-Specific','Primate-Specific','Species-Conserved')) + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test') + 
  theme(aspect.ratio = 0.6) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('expression') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')
dev.off()

## split Paola data and do not include Ho --------------------------------------------------------
#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_1.5 <- DotPlot(object = Ho_1.5_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_2 <- DotPlot(object = Ho_2_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_3 <- DotPlot(object = Ho_3_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola_4 <- DotPlot(object = Ho_4_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp"],p_Paola_1.5$data[p_Paola_1.5$data$id == 'RG-neu' & p_Paola_1.5$data$features.plot == y,"avg.exp"],p_Paola_2$data[p_Paola_2$data$id == 'RG-neu' & p_Paola_2$data$features.plot == y,"avg.exp"],p_Paola_3$data[p_Paola_3$data$id == 'RG-neu' & p_Paola_3$data$features.plot == y,"avg.exp"],p_Paola_4$data[p_Paola_4$data$id == 'RG-neu' & p_Paola_4$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = rep(x = get(x),times = 7),
                     gene_group = x,express = c(primary,organoid),
                     group = c(rep(x = 'primary',times = length(primary)),
                               rep(x = 'organoid',times = length(organoid))),
                     sample = c(rep(x = 'Trevino, et al',times = length(get(x))),rep(x = 'Polioudakis, et al',times = length(get(x))),rep(x = 'Pollen, et al',times = length(get(x))),
                                rep(x = 'Uzquiano, et al, 1.5m',times = length(get(x))),rep(x = 'Uzquiano, et al, 2m',times = length(get(x))),rep(x = 'Uzquiano, et al, 3m',times = length(get(x))),rep(x = 'Uzquiano, et al, 4m',times = length(get(x)))))
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))
pct_matrix$sample <- factor(pct_matrix$sample,
                            levels = c('Trevino, et al','Polioudakis, et al',
                                       'Pollen, et al',
                                       'Uzquiano, et al, 1.5m','Uzquiano, et al, 2m',
                                       'Uzquiano, et al, 3m','Uzquiano, et al, 4m'))

#plot
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = log(express),fill = sample,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  theme_ArchR() + 
  ylim(c(-3,3)) + 
  scale_fill_manual(values = c('#F7DC6F','#3498DB','#2ECC71','#E74C3C','#9B59B6','#F39C12','#1ABC9C')) + 
  scale_linetype_manual(values = c('solid','dashed')) + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test')

# featureplot of some L/R -------------------------------------------------
# col_list <- col_param$celltype[c("RG","Cycling",'IP',"Ex-1","InMGE","Mic")]
# names(col_list) <- c('RG-neu','Cycling','IP','Ex','In','others')

#LGALS3
data_list <- c('Greenleaf_RNA_Seurat','macaque_multiome_Seurat','mouse_multiome_Seurat')
for (i in 1:length(data_list)) {
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'ori_cell_type',label = FALSE,raster = FALSE,cols = col_param$celltype) + 
    theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
    labs(title = data_list[i]) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  if(i != length(data_list)){
    p <- p + NoLegend()
  }
  assign(x = paste0('p',i),value = p)

  #feature plot
  gene_list <- c('LGALS3','LGALS3BP')
  if(i == 3){
    gene_list <- str_to_title(gene_list)
  }
  
  for (j in 1:length(gene_list)) {
    indi <- TRUE
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],pt.size = 0.001,order = indi,raster = FALSE,cols = ArchRPalettes$whitePurple[-1]) + 
      theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
    if(i != length(data_list)){
      p <- p + NoLegend()
    }
    assign(x = paste0('p',i+j*length(data_list)),value = p)
  }
}

char <- paste0('ggarrange(',paste('p',1:9,sep = '',collapse = ','),',nrow = 3,ncol = 3)')

png(filename = './res/step_109_fig_230314/species_primary_LGALS3_LGALS3BP_expression_featureplot.png',width = 1400,height = 1200)
eval(parse(text = char))
dev.off()

#PDGF
data_list <- c('Greenleaf_RNA_Seurat','macaque_multiome_Seurat','mouse_multiome_Seurat')
for (i in 1:length(data_list)) {
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'ori_cell_type',label = FALSE,raster = FALSE,cols = col_param$celltype) + 
    theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
    labs(title = data_list[i]) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  if(i != length(data_list)){
    p <- p + NoLegend()
  }
  assign(x = paste0('p',i),value = p)
  
  #feature plot
  gene_list <- c('PDGFD','PDGFRB')
  if(i == 3){
    gene_list <- str_to_title(gene_list)
  }
  
  for (j in 1:length(gene_list)) {
    indi <- TRUE
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],pt.size = 0.001,order = indi,raster = FALSE,cols = ArchRPalettes$whitePurple[-1]) + 
      theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
    if(i != length(data_list)){
      p <- p + NoLegend()
    }
    assign(x = paste0('p',i+j*length(data_list)),value = p)
  }
}

char <- paste0('ggarrange(',paste('p',1:9,sep = '',collapse = ','),',nrow = 3,ncol = 3)')

png(filename = './res/step_109_fig_230314/species_primary_PDGFD_PDGFRB_expression_featureplot.png',width = 1400,height = 1200)
eval(parse(text = char))
dev.off()

#COL11A1 ITGA2
data_list <- c('Greenleaf_RNA_Seurat','macaque_multiome_Seurat','mouse_multiome_Seurat')
for (i in 1:length(data_list)) {
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'ori_cell_type',label = FALSE,raster = FALSE,cols = col_param$celltype) + 
    theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
    labs(title = data_list[i]) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  if(i != length(data_list)){
    p <- p + NoLegend()
  }
  assign(x = paste0('p',i),value = p)
  
  #feature plot
  gene_list <- c('COL11A1','ITGA2')
  if(i == 3){
    gene_list <- str_to_title(gene_list)
  }
  
  for (j in 1:length(gene_list)) {
    indi <- TRUE
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],pt.size = 0.001,order = indi,raster = FALSE,cols = ArchRPalettes$whitePurple[-1]) + 
      theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + 
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
    if(i != length(data_list)){
      p <- p + NoLegend()
    }
    assign(x = paste0('p',i+j*length(data_list)),value = p)
  }
}

char <- paste0('ggarrange(',paste('p',1:9,sep = '',collapse = ','),',nrow = 3,ncol = 3)')

png(filename = './res/step_109_fig_230314/species_primary_COL11A1_ITGA2_expression_featureplot.png',width = 1400,height = 1200)
eval(parse(text = char))
dev.off()

# seems hard to prove that Per is important