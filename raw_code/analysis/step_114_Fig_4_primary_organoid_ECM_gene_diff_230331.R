#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Fig_4_primary_organoid_ECM_gene_diff                            ##
## Data: 2023.03.31                                                                ##
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

# re-anno cell type -------------------------------------------------------
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

# merge RG-neu ------------------------------------------------------------

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

# expression percentage diff ----------------------------------------------

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
gene_list <- c('SFRP1','TNC','COL4A5','GPC4','COL11A1','BMP7','F3','ITGB5','LGALS3','ITGA2','ANXA2')
plot_data <- DotPlot(object = RG_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'sample_2')$data
plot_data$dataset <- gsub(pattern = ' et al.*',replacement = '',x = plot_data$id,fixed = FALSE)
plot_data$dataset <- paste(plot_data$dataset,'et al.',sep = ' ')
plot_data$dataset <- factor(plot_data$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
plot_data$group <- 'Organoid'
plot_data[grep(pattern = 'Trevino|Polioudakis|Bhaduri|Pollen. et al. Primary',x = plot_data$id,fixed = FALSE),"group"] <- 'Primary'

## no Pollen data ----------------------------------------------------------

for (i in gene_list) {
  pct_matrix <- plot_data[plot_data$features.plot == i & plot_data$dataset != 'Pollen. et al.' & !grepl(pattern = '^Velasco. et al. 6m',x = plot_data$id,fixed = FALSE),]
  pct_matrix %>% dplyr::summarize(mean_pct = mean(pct.exp),sd_pct = sd(pct.exp),.by = dataset) -> df_bar
  df_bar$dataset <- factor(df_bar$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
  p_value <- t.test(x = pct_matrix$pct.exp[pct_matrix$group == 'Primary'],y = pct_matrix$pct.exp[pct_matrix$group == 'Organoid'])$p.value
  p_value <- round(p_value,digits = 3)
  
  p <- ggplot(data = df_bar,mapping = aes(x = dataset,y = mean_pct,fill = dataset)) + 
    geom_bar(stat = 'identity',width = 0.6) + 
    geom_text(x = 3.5,y = max(df_bar$mean_pct),label = paste('p =',p_value,sep = ' '),size = 4,family = 'sans',fontface = 'italic') + 
    scale_fill_manual(values = c('#32CD32','#87CEEB','#B57EDC','#FFF44F','#FFC0CB')) + 
    geom_errorbar(mapping = aes(ymin = mean_pct - 0.5*sd_pct,ymax = mean_pct + 0.5*sd_pct),width = 0.2,color = 'black',linewidth = 0.2) + 
    theme_ArchR() + theme(aspect.ratio = 2) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 8)) + 
    labs(title = i) + CenterTitle()
  
  char <- paste0('./res/step_114_fig_230331/express_pct/primary_organoid_',i,'_expression_pct_barplot.pdf')
  pdf(file = char,width = 5,height = 5)
  print(p)
  dev.off()
}

#does not seems good

## Pollen data -------------------------------------------------------------
for (i in gene_list) {
  pct_matrix <- plot_data[plot_data$features.plot == i & plot_data$dataset == 'Pollen. et al.',]
  pct_matrix %>% dplyr::summarize(mean_pct = mean(pct.exp),sd_pct = sd(pct.exp),.by = group) -> df_bar
  df_bar$group <- factor(df_bar$group,levels = c('Primary','Organoid'))
  p_value <- t.test(x = pct_matrix$pct.exp[pct_matrix$group == 'Primary'],y = pct_matrix$pct.exp[pct_matrix$group == 'Organoid'])$p.value
  p_value <- format(p_value,scientific = TRUE,digits = 3)
  
  p <- ggplot(data = df_bar,mapping = aes(x = group,y = mean_pct,fill = group)) + 
    geom_bar(stat = 'identity',width = 0.6) + 
    geom_text(x = 1.5,y = max(df_bar$mean_pct),label = paste('p =',p_value,sep = ' '),size = 4,family = 'sans',fontface = 'italic') + 
    scale_fill_manual(values = c('#FF5733','#00BFFF')) + 
    geom_errorbar(mapping = aes(ymin = mean_pct - 0.5*sd_pct,ymax = mean_pct + 0.5*sd_pct),width = 0.2,color = 'black',linewidth = 0.2) + 
    theme_ArchR() + theme(aspect.ratio = 4) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 8)) + 
    labs(title = i) + CenterTitle()
  
  char <- paste0('./res/step_114_fig_230331/express_pct/Pollen_primary_organoid_',i,'_expression_pct_barplot.pdf')
  pdf(file = char,width = 5,height = 5)
  print(p)
  dev.off()
}

# expression level diff ----------------------------------------------

#calculate percentage
gene_list <- c('SFRP1','TNC','COL4A5','GPC4','COL11A1','BMP7','F3','ITGB5','LGALS3','ITGA2','ANXA2')
plot_data <- DotPlot(object = RG_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'sample_2')$data
plot_data$dataset <- gsub(pattern = ' et al.*',replacement = '',x = plot_data$id,fixed = FALSE)
plot_data$dataset <- paste(plot_data$dataset,'et al.',sep = ' ')
plot_data$dataset <- factor(plot_data$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
plot_data$group <- 'Organoid'
plot_data[grep(pattern = 'Trevino|Polioudakis|Bhaduri|Pollen. et al. Primary',x = plot_data$id,fixed = FALSE),"group"] <- 'Primary'

## no Pollen data ----------------------------------------------------------

for (i in gene_list) {
  pct_matrix <- plot_data[plot_data$features.plot == i & plot_data$dataset != 'Pollen. et al.' & !grepl(pattern = '^Velasco. et al. 6m',x = plot_data$id,fixed = FALSE),]
  pct_matrix %>% dplyr::summarize(mean_pct = mean(avg.exp),sd_pct = sd(avg.exp),.by = dataset) -> df_bar
  df_bar$dataset <- factor(df_bar$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
  p_value <- t.test(x = pct_matrix$avg.exp[pct_matrix$group == 'Primary'],y = pct_matrix$avg.exp[pct_matrix$group == 'Organoid'])$p.value
  p_value <- round(p_value,digits = 3)
  
  p <- ggplot(data = df_bar,mapping = aes(x = dataset,y = mean_pct,fill = dataset)) + 
    geom_bar(stat = 'identity',width = 0.6) + 
    geom_text(x = 3.5,y = max(df_bar$mean_pct),label = paste('p =',p_value,sep = ' '),size = 4,family = 'sans',fontface = 'italic') + 
    scale_fill_manual(values = c('#32CD32','#87CEEB','#B57EDC','#FFF44F','#FFC0CB')) + 
    geom_errorbar(mapping = aes(ymin = mean_pct - 0.5*sd_pct,ymax = mean_pct + 0.5*sd_pct),width = 0.2,color = 'black',linewidth = 0.2) + 
    theme_ArchR() + theme(aspect.ratio = 2) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 8)) + 
    labs(title = i) + CenterTitle()
  
  char <- paste0('./res/step_114_fig_230331/express_level/primary_organoid_',i,'_expression_level_barplot.pdf')
  pdf(file = char,width = 5,height = 5)
  print(p)
  dev.off()
}

## Pollen data -------------------------------------------------------------
for (i in gene_list) {
  pct_matrix <- plot_data[plot_data$features.plot == i & plot_data$dataset == 'Pollen. et al.',]
  pct_matrix %>% dplyr::summarize(mean_pct = mean(avg.exp),sd_pct = sd(avg.exp),.by = group) -> df_bar
  df_bar$group <- factor(df_bar$group,levels = c('Primary','Organoid'))
  p_value <- t.test(x = pct_matrix$avg.exp[pct_matrix$group == 'Primary'],y = pct_matrix$avg.exp[pct_matrix$group == 'Organoid'])$p.value
  p_value <- round(p_value,digits = 3)
  
  p <- ggplot(data = df_bar,mapping = aes(x = group,y = mean_pct,fill = group)) + 
    geom_bar(stat = 'identity',width = 0.6) + 
    geom_text(x = 1.5,y = max(df_bar$mean_pct),label = paste('p =',p_value,sep = ' '),size = 4,family = 'sans',fontface = 'italic') + 
    scale_fill_manual(values = c('#FF5733','#00BFFF')) + 
    geom_errorbar(mapping = aes(ymin = mean_pct - 0.5*sd_pct,ymax = mean_pct + 0.5*sd_pct),width = 0.2,color = 'black',linewidth = 0.2) + 
    theme_ArchR() + theme(aspect.ratio = 4) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 8)) + 
    labs(title = i) + CenterTitle()
  
  char <- paste0('./res/step_114_fig_230331/express_level/Pollen_primary_organoid_',i,'_expression_level_barplot.pdf')
  pdf(file = char,width = 5,height = 5)
  print(p)
  dev.off()
}

# output expression level diff --------------------------------------------
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
gene_list <- c('ITGA2','ANXA2','BMP7','COL4A5','LGALS3','TNC')
plot_data <- DotPlot(object = RG_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'sample_2')$data
plot_data$dataset <- gsub(pattern = ' et al.*',replacement = '',x = plot_data$id,fixed = FALSE)
plot_data$dataset <- paste(plot_data$dataset,'et al.',sep = ' ')
plot_data$dataset <- factor(plot_data$dataset,levels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al.','Uzquiano. et al.','Velasco. et al.'))
plot_data$group <- 'Organoid'
plot_data[grep(pattern = 'Trevino|Polioudakis|Bhaduri|Pollen. et al. Primary',x = plot_data$id,fixed = FALSE),"group"] <- 'Primary'

## no Pollen ---------------------------------------------------------------
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
  
  char <- paste0('./res/step_114_fig_230331/output/primary_organoid_',i,'_expression_level_barplot.pdf')
  pdf(file = char,width = 5,height = 5)
  print(p)
  dev.off()
}
