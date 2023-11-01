#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Fig4_check_public_data                                          ##
## Data: 2023.03.26                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session()

# load data ---------------------------------------------------------------

#load col param
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')
PFC_RNA_Seurat <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_PFC_brain_preprocess_PC30_by_lyt_230326.rds')

gc()

#preprocess PD human
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,npcs = 50,preprocess = TRUE)
PD_RNA_Seurat <- my_process_seurat(object = PD_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'Cluster',label = TRUE)

#normalize data
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  temp <- get(data_list[i])
  temp <- NormalizeData(temp)
  assign(x = data_list[i],value = temp)
  rm(temp)
  gc()
}

# check cell type ---------------------------------------------------------

#ori_anno
Greenleaf_RNA_Seurat$ori_anno <- Greenleaf_RNA_Seurat$ReAnno_celltype
PD_RNA_Seurat$ori_anno <- PD_RNA_Seurat$Cluster
PFC_RNA_Seurat$ori_anno <- PFC_RNA_Seurat$cell.type
Hp_RNA_Seurat$ori_anno <- Hp_RNA_Seurat$cell_type
Ho_RNA_Seurat$ori_anno <- Ho_RNA_Seurat$cell_type
Ho_1.5_Seurat$ori_anno <- Ho_1.5_Seurat$FinalName
Ho_2_Seurat$ori_anno <- Ho_2_Seurat$FinalName
Ho_3_Seurat$ori_anno <- Ho_3_Seurat$FinalName
Ho_4_Seurat$ori_anno <- Ho_4_Seurat$FinalName

#featureplot
gene_list <- c('SOX9','VIM','PAX6','CRYAB','HOPX','EGFR','EOMES','PPP1R17','NEUROD2','NEUROD6','GAD2','DLX5','OLIG2')
for (i in 1:length(data_list)) {
  #dimplot
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'ori_anno',label = TRUE,repel = TRUE,raster = FALSE) + 
    theme_classic() + theme(aspect.ratio = 1) + 
    labs(title = data_list[i]) + CenterTitle() + NoLegend() + NoAxes() + 
    theme(panel.background = element_rect(fill = NA,color = 'black',linewidth = 0.5))
  
  assign(x = paste('p',i,sep = '_'),value = p)
  
  #feature plot
  for (j in 1:length(gene_list)) {
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],cols = ArchRPalettes$whitePurple[-1],pt.size = 0.001,order = TRUE,raster = FALSE) + 
      theme_classic() + theme(aspect.ratio = 1) + 
      labs(title = gene_list[j]) + CenterTitle() + NoLegend() + NoAxes() + 
      theme(panel.background = element_rect(fill = NA,color = 'black',linewidth = 0.5))
    
    assign(x = paste('p',i+j*length(data_list),sep = '_'),value = p)
  }
}

#plot
png(filename = './res/step_112_fig_230326/check_cell_type_classic_marker_featureplot.png',width = 9*400,height = 14*420)
eval(parse(text = paste('p',1:126,sep = '_') %>% paste(collapse = '+') %>% paste('plot_layout(ncol = 9)',sep = '+')))
dev.off()

#featureplot
gene_list <- c('SOX9','VIM','PAX6','CRYAB','HOPX','EGFR','EOMES','PPP1R17','NEUROD2','NEUROD6','GAD2','DLX5','OLIG2')
for (i in 1:length(data_list)) {
  #dimplot
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'ori_anno',label = TRUE,repel = TRUE,raster = FALSE) + 
    theme_classic() + theme(aspect.ratio = 1) + 
    labs(title = data_list[i]) + CenterTitle() + NoLegend() + NoAxes() + 
    theme(panel.background = element_rect(fill = NA,color = 'black',linewidth = 0.5))
  
  assign(x = paste('p',i,sep = '_'),value = p)
  
  #feature plot
  for (j in 1:length(gene_list)) {
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],cols = ArchRPalettes$whitePurple[-1],pt.size = 0.001,order = FALSE,raster = FALSE) + 
      theme_classic() + theme(aspect.ratio = 1) + 
      labs(title = gene_list[j]) + CenterTitle() + NoLegend() + NoAxes() + 
      theme(panel.background = element_rect(fill = NA,color = 'black',linewidth = 0.5))
    
    assign(x = paste('p',i+j*length(data_list),sep = '_'),value = p)
  }
}

#plot
png(filename = './res/step_112_fig_230326/check_cell_type_classic_marker_featureplot_no_order.png',width = 9*400,height = 14*420)
eval(parse(text = paste('p',1:126,sep = '_') %>% paste(collapse = '+') %>% paste('plot_layout(ncol = 9)',sep = '+')))
dev.off()

# RG-1 marker dotplot -----------------------------------------------------

#load gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_1.5_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_1.5_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_1.5_Seurat)]

#add temp_anno
Greenleaf_RNA_Seurat$temp_anno <- Greenleaf_RNA_Seurat$ori_anno
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ori_anno %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"temp_anno"] <- 'Ex'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ori_anno %in% c('InCGE','InMGE'),"temp_anno"] <- 'In'

PD_RNA_Seurat$temp_anno <- PD_RNA_Seurat$ori_anno
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$ori_anno %in% c('ExDp1','ExDp2','ExM','ExM-U','ExN'),"temp_anno"] <- 'Ex'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$ori_anno %in% c('InCGE','InMGE'),"temp_anno"] <- 'In'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$ori_anno %in% c('PgG2M','PgS'),"temp_anno"] <- 'Cycling'

PFC_RNA_Seurat$temp_anno <- PFC_RNA_Seurat$ori_anno
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('Dividing'),"temp_anno"] <- 'Cycling'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('Endo'),"temp_anno"] <- 'End'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('Interneuron'),"temp_anno"] <- 'In'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('IPC'),"temp_anno"] <- 'IP'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('Microglia'),"temp_anno"] <- 'Mic'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('Neuron'),"temp_anno"] <- 'Ex'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$ori_anno %in% c('Oligo'),"temp_anno"] <- 'OPC'

Hp_RNA_Seurat$temp_anno <- Hp_RNA_Seurat$ori_anno
Ho_RNA_Seurat$temp_anno <- Ho_RNA_Seurat$ori_anno

data_list <- c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  temp <- get(data_list[i])
  temp$temp_anno <- as.character(temp$ori_anno)
  temp@meta.data[temp$ori_anno %in% c('CFuPN','CPN','Newborn CFuPN','Newborn CPN','Newborn DL PN','Newborn PN','PN','Subcortical neurons','Subcortical progenitors'),"temp_anno"] <- 'Ex'
  temp@meta.data[temp$ori_anno %in% c('Immature IN'),"temp_anno"] <- 'In'
  assign(x = data_list[i],value = temp)
  rm(temp)
  gc()
}

#cell type list
table(Greenleaf_RNA_Seurat$temp_anno)
table(PD_RNA_Seurat$temp_anno)
table(Hp_RNA_Seurat$temp_anno)
table(Ho_1.5_Seurat$temp_anno)
table(Ho_4_Seurat$temp_anno)
cell_type_list <- c('RG-1','RG-2','RG','oRG','vRG','aRG','IP','Ex','In')

#dotplot
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (g in c('HG_gene_list','PC_gene_list','SC_gene_list')) {
  for (i in 1:length(data_list)) {
    #set color
    if(i <= 4){
      temp_col <- '#FF5733'
    }else{
      temp_col <- '#00BFFF'
    }
    #get p
    p <- DotPlot(object = get(data_list[i]),assay = 'RNA',features = get(g),cols = c('lightgrey',temp_col),
                 group.by = 'temp_anno')
    p$data <- p$data[p$data$id %in% cell_type_list,]
    p$data$id <- factor(p$data$id,levels = cell_type_list)
    #beauty p
    p <- p + coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
      guides(size = guide_legend(title = "Percent Expressed")) + 
      theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) + 
      theme(aspect.ratio = (length(get(g)) + 0.2)/(length(unique(p$data$id)) + 0.2),
            text = element_text(size = 10,family = 'sans')) + 
      labs(title = data_list[i])
    #whether show legand
    if(i != length(data_list)){
      p <- p + NoLegend()
    }
    if(i != 1){
      p <- p + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
    }
    #assign
    assign(x = paste('p',i,sep = '_'),value = p)
  }
  #plot
  char <- paste0('./res/step_112_fig_230326/primary_organoid_',g,'_dotplot.pdf')
  pdf(file = char,width = 15,height = length(get(g)) * 0.4)
  print(eval(parse(text = paste('p',1:9,sep = '_') %>% paste(collapse = '+') %>% paste('plot_layout(ncol = 9)',sep = '+'))))
  dev.off()
}