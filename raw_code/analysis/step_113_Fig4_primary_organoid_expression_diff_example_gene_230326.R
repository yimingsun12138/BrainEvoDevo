#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Fig4_primary_organoid_expression_diff_example_gene              ##
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

# pre-process PFC data ----------------------------------------------------
#load data
PFC_RNA_Seurat <- readRDS(file = './data/public/An_atlas_of_cortical_arealization_identifies_dynamic_molecular_signatures/human_PFC_brain_preprocess_PC30_by_lyt_230326.rds')

#RG only needed
cell_list <- colnames(PFC_RNA_Seurat)[PFC_RNA_Seurat$cell.type == 'RG']
meta_data <- PFC_RNA_Seurat@meta.data[cell_list,c("cell.name","age","individual","structure","area","cell.type","ngene","numi")]
PFC_RNA_Seurat <- PFC_RNA_Seurat@assays$RNA@counts[,cell_list]
PFC_RNA_Seurat <- CreateSeuratObject(counts = PFC_RNA_Seurat,project = 'PFC',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

PFC_RNA_Seurat <- my_process_seurat(object = PFC_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',vars.to.regress = c('age','individual','nCount_RNA'),nfeatures = 2000,npcs = 50,preprocess = TRUE)
PFC_RNA_Seurat <- my_process_seurat(object = PFC_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 20,resolution = 0.5,group.by = 'individual',label = TRUE)

#check marker
FeaturePlot(object = PFC_RNA_Seurat,features = c('PAX6','EGFR','OLIG2','KCNE5','CHST7','TRIB2','DUSP6','PDGFRA','MAP3K1'),order = TRUE)
DimPlot(object = PFC_RNA_Seurat,group.by = 'individual',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)
DimPlot(object = PFC_RNA_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

#add anno
PFC_RNA_Seurat$cell_type <- 'RG-neu'
PFC_RNA_Seurat@meta.data[PFC_RNA_Seurat$seurat_clusters %in% c('5','6','8'),"cell_type"] <- 'RG-glia'

#save meta data
saveRDS(object = PFC_RNA_Seurat@meta.data,file = './res/step_113_fig_230326/PFC_RG_meta_data.rds')

#dimplot
DimPlot(object = PFC_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

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

# #filter
# Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$ori_anno %in% c('RG-1')]
# PD_RNA_Seurat <- PD_RNA_Seurat[,PD_RNA_Seurat$ori_anno %in% c('oRG','vRG')]
# PFC_RNA_Seurat <- PFC_RNA_Seurat[,PFC_RNA_Seurat$ori_anno %in% c('RG-neu')]
# Hp_RNA_Seurat <- Hp_RNA_Seurat[,Hp_RNA_Seurat$cell_type == 'RG']
# Ho_RNA_Seurat <- Ho_RNA_Seurat[,Ho_RNA_Seurat$cell_type == 'RG']
# Ho_1.5_Seurat <- Ho_1.5_Seurat[,Ho_1.5_Seurat$ori_anno %in% c('aRG','oRG')]
# Ho_2_Seurat <- Ho_2_Seurat[,Ho_2_Seurat$ori_anno %in% c('aRG','oRG')]
# Ho_3_Seurat <- Ho_3_Seurat[,Ho_3_Seurat$ori_anno %in% c('aRG','oRG')]
# Ho_4_Seurat <- Ho_4_Seurat[,Ho_4_Seurat$ori_anno %in% c('aRG','oRG')]
# Paola_3_Seurat <- Paola_3_Seurat[,Paola_3_Seurat$ori_anno %in% c('oRG','RG')]
# Paola_6_Seurat <- Paola_6_Seurat[,Paola_6_Seurat$ori_anno %in% c('oRG','RG')]
# gc()

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

# merge none 10x dataset --------------------------------------------------
RG_Pollen_Seurat <- cbind(Hp_RNA_Seurat@assays$RNA@data[gene_list,Hp_RNA_Seurat$cell_type == 'RG'],
                          Ho_RNA_Seurat@assays$RNA@data[gene_list,Ho_RNA_Seurat$cell_type == 'RG'])
RG_Pollen_Seurat <- expm1(RG_Pollen_Seurat)
meta_data <- data.frame(cell_name = colnames(RG_Pollen_Seurat),
                        dataset = 'Pollen. et al.',
                        group = c(rep('Primary',times = sum(Hp_RNA_Seurat$cell_type == 'RG')),
                                  rep('Organoid',times = sum(Ho_RNA_Seurat$cell_type == 'RG'))))
rownames(meta_data) <- meta_data$cell_name
RG_Pollen_Seurat <- CreateSeuratObject(counts = RG_Pollen_Seurat,project = 'Pollen',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
RG_Pollen_Seurat$sample <- paste(RG_Pollen_Seurat$dataset,RG_Pollen_Seurat$group,sep = ' ')

RG_Pollen_Seurat$group <- factor(RG_Pollen_Seurat$group,levels = c('Primary','Organoid'))
RG_Pollen_Seurat$sample <- factor(RG_Pollen_Seurat$sample,levels = c('Pollen. et al. Primary','Pollen. et al. Organoid'))

RG_Pollen_Seurat <- NormalizeData(RG_Pollen_Seurat)

# RG marker dotplot -------------------------------------------------------
pal.GEX <- rev(paletteer_c("grDevices::Blues 3", 10))
pal.GEX[1] <- "#e0ecf4"

HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(RG_RNA_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(RG_RNA_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(RG_RNA_Seurat)]

for (i in c('HG_gene_list','PC_gene_list','SC_gene_list')) {
  gene_list <- get(i)
  p <- DotPlot(object = RG_RNA_Seurat[,!(RG_RNA_Seurat$dataset %in% c('Pollen. et al. Primary','Pollen. et al. Organoid'))],assay = 'RNA',features = gene_list,group.by = 'dataset',cols = c('white','#00366CFF'))
  p <- p + coord_flip() + RotatedAxis() + 
    guides(color = guide_colorbar(title = "Average Expression")) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    theme(aspect.ratio = (length(gene_list) + 0.2)/(5*1.6 + 0.2),
          text = element_text(size = 10,family = 'sans')) + 
    labs(title = '')
  
  char <- paste0('./res/step_113_fig_230326/',i,'_dotplot_for_gene_pick.pdf')
  pdf(file = char,width = 8,height = length(gene_list) * 0.25)
  print(p)
  dev.off()
}

#set gene list
HG_gene_list <- c('PEA15','GATM','LGALS3','ATP1B2','ANXA5','KLF6','GPX3','ITGA2','TPM2',
                  'PDGFD','TFPI','LDLR','LRRC17','ANXA6','CAST','SOAT1','ANXA2','STOX1')
PC_gene_list <- c('CLU','HOPX','FAM107A','COL11A1','BMP7','ITGB5','GFAP','SLC7A11','SLC4A4')
SC_gene_list <- c('TNC','PTN','VIM','FGFR1','SOX6','SOX3','CD9','ID4','YAP1')

for (i in c('HG_gene_list','PC_gene_list','SC_gene_list')) {
  gene_list <- get(i)
  p <- DotPlot(object = RG_RNA_Seurat[,!(RG_RNA_Seurat$dataset %in% c('Pollen. et al. Primary','Pollen. et al. Organoid'))],assay = 'RNA',features = gene_list,group.by = 'dataset',cols = c('white','#00366CFF'))
  p <- p + coord_flip() + RotatedAxis() + 
    guides(color = guide_colorbar(title = "Average Expression")) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    theme(aspect.ratio = (length(gene_list) + 0.2)/(5*1.6 + 0.2),
          text = element_text(size = 10,family = 'sans')) + 
    labs(title = '') + NoLegend()
  
  char <- paste0('./res/step_113_fig_230326/',i,'_dotplot_for_show.pdf')
  pdf(file = char,width = 3,height = 30 * 0.25)
  print(p)
  dev.off()
}

# statistic ---------------------------------------------------------------

## average express ---------------------------------------------------------
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(RG_RNA_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(RG_RNA_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(RG_RNA_Seurat)]

#expression
exp_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(g){
  gene_list <- get(g)
  p <- DotPlot(object = RG_RNA_Seurat[,!(RG_RNA_Seurat$dataset %in% c('Pollen. et al. Primary','Pollen. et al. Organoid'))],assay = 'RNA',features = gene_list,group.by = 'sample')
  
  primary <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Trevino|Polioudakis|Bhaduri',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"avg.exp"]))
  })
  organoid <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Uzquiano|Velasco',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"avg.exp"]))
  })
  
  temp <- data.frame(gene_list = gene_list,primary = unlist(primary),organoid = unlist(organoid),gene_group = g)
  return(temp)
}))

exp_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = exp_matrix$gene_group,fixed = TRUE)
exp_matrix <- gather(data = exp_matrix,key = 'group',value = 'exp','primary','organoid')
exp_matrix$gene_group <- factor(exp_matrix$gene_group,levels = c('HG','PC','SC'))
exp_matrix$group <- factor(exp_matrix$group,levels = c('primary','organoid'))

p <- ggplot(data = exp_matrix,mapping = aes(x = gene_group,y = log(exp),fill = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  coord_cartesian(ylim = c(-2.5,2.5)) + 
  scale_fill_manual(values = c('#FF5733','#00BFFF')) + 
  theme_ArchR() + 
  theme(aspect.ratio = 1.4,
        legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  scale_x_discrete(labels = c('Human-Specific','Primate-Specific','Species-Conserved')) + 
  theme(axis.title.x = element_blank()) + 
  ylab('log(Average expession)') + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test',
                     label.y = c(2.5,2.5,2.5)) + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')

pdf(file = './res/step_113_fig_230326/why_give_up/average_expression_boxplot.pdf',width = 6,height = 6)
p
dev.off()

#expression difference
exp_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(g){
  gene_list <- get(g)
  p <- DotPlot(object = RG_RNA_Seurat[,!(RG_RNA_Seurat$dataset %in% c('Pollen. et al. Primary','Pollen. et al. Organoid'))],assay = 'RNA',features = gene_list,group.by = 'sample')
  
  primary <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Trevino|Polioudakis|Bhaduri',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"avg.exp"]))
  })
  organoid <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Uzquiano|Velasco',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"avg.exp"]))
  })
  
  temp <- data.frame(gene_list = gene_list,primary = unlist(primary),organoid = unlist(organoid),gene_group = g)
  return(temp)
}))

exp_matrix[exp_matrix$gene_group == 'HG_gene_list',"gene_group"] <- 'Human-Specific'
exp_matrix[exp_matrix$gene_group == 'PC_gene_list',"gene_group"] <- 'Primate-Specific'
exp_matrix[exp_matrix$gene_group == 'SC_gene_list',"gene_group"] <- 'Species-Conserved'
exp_matrix$gene_group <- factor(exp_matrix$gene_group,levels = c('Human-Specific','Primate-Specific','Species-Conserved'))
exp_matrix$diff <- exp_matrix$primary - exp_matrix$organoid

p <- ggplot(data = exp_matrix,mapping = aes(x = gene_group,y = diff,fill = gene_group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  coord_cartesian(ylim = c(-2.5,2.5)) + 
  scale_fill_manual(values = as.character(col_param$species[c("human","macaque","mouse")])) + 
  theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red',linewidth = 0.2) + 
  stat_compare_means(comparisons = list(c(2,3),c(1,2),c(1,3)),
                     label = 'p.format',method = 't.test',
                     label.y = c(-1.2,-0.9,-0.6),tip.length = 0.001) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('expression difference') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')

pdf(file = './res/step_113_fig_230326/why_give_up/average_expression_difference_boxplot.pdf',width = 6,height = 6)
p
dev.off()

## express percent ---------------------------------------------------------
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(RG_RNA_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(RG_RNA_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(RG_RNA_Seurat)]

#expression
exp_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(g){
  gene_list <- get(g)
  p <- DotPlot(object = RG_RNA_Seurat[,!(RG_RNA_Seurat$dataset %in% c('Pollen. et al. Primary','Pollen. et al. Organoid'))],assay = 'RNA',features = gene_list,group.by = 'sample')
  
  primary <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Trevino|Polioudakis|Bhaduri',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"pct.exp"]))
  })
  organoid <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Uzquiano|Velasco',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"pct.exp"]))
  })
  
  temp <- data.frame(gene_list = gene_list,primary = unlist(primary),organoid = unlist(organoid),gene_group = g)
  return(temp)
}))

exp_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = exp_matrix$gene_group,fixed = TRUE)
exp_matrix <- gather(data = exp_matrix,key = 'group',value = 'exp','primary','organoid')
exp_matrix$gene_group <- factor(exp_matrix$gene_group,levels = c('HG','PC','SC'))
exp_matrix$group <- factor(exp_matrix$group,levels = c('primary','organoid'))

p <- ggplot(data = exp_matrix,mapping = aes(x = gene_group,y = exp,fill = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  coord_cartesian(ylim = c(0,75)) + 
  scale_fill_manual(values = c('#FF5733','#00BFFF')) + 
  theme_ArchR() + 
  theme(aspect.ratio = 1.4,
        legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  scale_x_discrete(labels = c('Human-Specific','Primate-Specific','Species-Conserved')) + 
  theme(axis.title.x = element_blank()) + 
  ylab('expression percent') + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test',
                     label.y = c(75,75,75)) + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')

pdf(file = './res/step_113_fig_230326/why_give_up/expression_percent_boxplot.pdf',width = 6,height = 6)
p
dev.off()

#expression difference
exp_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(g){
  gene_list <- get(g)
  p <- DotPlot(object = RG_RNA_Seurat[,!(RG_RNA_Seurat$dataset %in% c('Pollen. et al. Primary','Pollen. et al. Organoid'))],assay = 'RNA',features = gene_list,group.by = 'sample')
  
  primary <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Trevino|Polioudakis|Bhaduri',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"pct.exp"]))
  })
  organoid <- base::lapply(X = gene_list,FUN = function(x){
    return(mean(p$data[grepl(pattern = 'Uzquiano|Velasco',x = p$data$id,fixed = FALSE) & p$data$features.plot == x,"pct.exp"]))
  })
  
  temp <- data.frame(gene_list = gene_list,primary = unlist(primary),organoid = unlist(organoid),gene_group = g)
  return(temp)
}))

exp_matrix[exp_matrix$gene_group == 'HG_gene_list',"gene_group"] <- 'Human-Specific'
exp_matrix[exp_matrix$gene_group == 'PC_gene_list',"gene_group"] <- 'Primate-Specific'
exp_matrix[exp_matrix$gene_group == 'SC_gene_list',"gene_group"] <- 'Species-Conserved'
exp_matrix$gene_group <- factor(exp_matrix$gene_group,levels = c('Human-Specific','Primate-Specific','Species-Conserved'))
exp_matrix$diff <- exp_matrix$primary - exp_matrix$organoid

p <- ggplot(data = exp_matrix,mapping = aes(x = gene_group,y = diff,fill = gene_group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  coord_cartesian(ylim = c(-20,40)) + 
  scale_fill_manual(values = as.character(col_param$species[c("human","macaque","mouse")])) + 
  theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red',linewidth = 0.2) + 
  stat_compare_means(comparisons = list(c(2,3),c(1,2),c(1,3)),
                     label = 'p.format',method = 't.test',
                     tip.length = 0.02,label.y = c(25,30,35)) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('expression percent difference') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')

pdf(file = './res/step_113_fig_230326/why_give_up/expression_percent_difference_boxplot.pdf',width = 6,height = 6)
p
dev.off()


## split sample ------------------------------------------------------------
#by dataset
for (i in c('HG_gene_list','PC_gene_list','SC_gene_list')) {
  gene_list <- get(i)
  p <- DotPlot(object = RG_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'dataset')
  p$data$group <- NA
  p$data[grep(pattern = 'Trevino|Polioudakis|Bhaduri|Pollen. et al. Primary',x = p$data$id,fixed = FALSE),"group"] <- 'primary'
  p$data[grep(pattern = 'Uzquiano|Velasco|Pollen. et al. Organoid',x = p$data$id,fixed = FALSE),"group"] <- 'organoid'
  p$data$group <- factor(p$data$group,levels = c('primary','organoid'))
  
  #primary
  p <- ggplot(data = p$data,mapping = aes(x = id,y = avg.exp,fill = id,linetype = group)) + 
    geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
    coord_cartesian(ylim = c(0,5)) + 
    theme_ArchR() + 
    scale_linetype_manual(values = c('solid','dashed')) + 
    theme(axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
    theme(aspect.ratio = 8/length(unique(p$data$id))) + 
    labs(title = i) + CenterTitle() + NoLegend()
  
  if(i != 'HG_gene_list'){
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  
  assign(x = paste('p',i,sep = '_'),value = p)
}

pdf(file = './res/step_113_fig_230326/why_give_up/average_expression_boxplot_split_dataset.pdf',width = 15,height = 6)
p_HG_gene_list + p_PC_gene_list + p_SC_gene_list + plot_layout(ncol = 3)
dev.off()

#by sample
for (i in c('HG_gene_list','PC_gene_list','SC_gene_list')) {
  gene_list <- get(i)
  p <- DotPlot(object = RG_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'sample')
  p$data$group <- NA
  p$data[grep(pattern = 'Trevino|Polioudakis|Bhaduri|Pollen. et al. Primary',x = p$data$id,fixed = FALSE),"group"] <- 'primary'
  p$data[grep(pattern = 'Uzquiano|Velasco|Pollen. et al. Organoid',x = p$data$id,fixed = FALSE),"group"] <- 'organoid'
  p$data$group <- factor(p$data$group,levels = c('primary','organoid'))
  
  #primary
  p <- ggplot(data = p$data,mapping = aes(x = id,y = avg.exp,fill = id,linetype = group)) + 
    geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
    coord_cartesian(ylim = c(0,5)) + 
    theme_ArchR() + 
    scale_linetype_manual(values = c('solid','dashed')) + 
    theme(axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
    theme(aspect.ratio = 8/length(unique(p$data$id))) + 
    labs(title = i) + CenterTitle() + NoLegend()
  
  if(i != 'HG_gene_list'){
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  
  assign(x = paste('p',i,sep = '_'),value = p)
}

pdf(file = './res/step_113_fig_230326/why_give_up/HG_average_expression_boxplot_split_sample.pdf',width = 12,height = 6)
p_HG_gene_list
dev.off()

# ECM signature -----------------------------------------------------------

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

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
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
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat'),"group"] <- 'organoid'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('PFC_RNA_Seurat'),"author"] <- 'Bhaduri. et al'
temp[temp$dataset %in% c('Hp_RNA_Seurat','Ho_RNA_Seurat'),"author"] <- 'Pollen. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'
temp[temp$dataset %in% c('Paola_3_Seurat','Paola_6_Seurat'),"author"] <- 'Velasco. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = data_list)
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Bhaduri. et al','Pollen. et al','Uzquiano. et al','Velasco. et al'))

p1 <- ggplot(data = temp[temp$group == 'primary',],mapping = aes(x = author,y = ECM1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.1,0.15)) + 
  #scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00','#39CCCC')) + 
  theme(aspect.ratio = 7/4.2) + 
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
  scale_fill_manual(values = c('#39CCCC','#0074D9','#2ECC40')) + 
  theme(aspect.ratio = 7/3.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('ECM module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_113_fig_230326/primary_organoid_ECM_signature_boxplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

t.test(x = temp[temp$group == 'primary',"ECM1"],
       y = temp[temp$group == 'organoid',"ECM1"],
       alternative = 'greater')

# HG signature ------------------------------------------------------------
#get gene list
gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat) & gene_list %in% rownames(PD_RNA_Seurat) & gene_list %in% rownames(PFC_RNA_Seurat) & gene_list %in% rownames(Hp_RNA_Seurat) & gene_list %in% rownames(Ho_RNA_Seurat) & gene_list %in% rownames(Ho_1.5_Seurat) & gene_list %in% rownames(Ho_2_Seurat) & gene_list %in% rownames(Ho_3_Seurat) & gene_list %in% rownames(Ho_4_Seurat) & gene_list %in% rownames(Paola_3_Seurat) & gene_list %in% rownames(Paola_6_Seurat)]

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in data_list) {
  assign(x = i,value = AddModuleScore(object = get(i),features = list(gene_list),name = 'HG'))
  print(paste(i,'done!',sep = ' '))
  gc()
}

#plot
temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$re_anno %in% 'RG-neu',]
  meta_data <- meta_data[,c("re_anno","HG1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$group <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat'),"group"] <- 'organoid'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('PFC_RNA_Seurat'),"author"] <- 'Bhaduri. et al'
temp[temp$dataset %in% c('Hp_RNA_Seurat','Ho_RNA_Seurat'),"author"] <- 'Pollen. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'
temp[temp$dataset %in% c('Paola_3_Seurat','Paola_6_Seurat'),"author"] <- 'Velasco. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = data_list)
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Bhaduri. et al','Pollen. et al','Uzquiano. et al','Velasco. et al'))

p1 <- ggplot(data = temp[temp$group == 'primary',],mapping = aes(x = author,y = HG1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.25,1)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00','#39CCCC')) + 
  theme(aspect.ratio = 7/4.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('HG module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"HG1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p2 <- ggplot(data = temp[temp$group == 'organoid',],mapping = aes(x = author,y = HG1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.25,1)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#39CCCC','#0074D9','#2ECC40')) + 
  theme(aspect.ratio = 7/3.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('HG module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"HG1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_113_fig_230326/primary_organoid_HG_signature_boxplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

# PC signature ------------------------------------------------------------
#get gene list
gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat) & gene_list %in% rownames(PD_RNA_Seurat) & gene_list %in% rownames(PFC_RNA_Seurat) & gene_list %in% rownames(Hp_RNA_Seurat) & gene_list %in% rownames(Ho_RNA_Seurat) & gene_list %in% rownames(Ho_1.5_Seurat) & gene_list %in% rownames(Ho_2_Seurat) & gene_list %in% rownames(Ho_3_Seurat) & gene_list %in% rownames(Ho_4_Seurat) & gene_list %in% rownames(Paola_3_Seurat) & gene_list %in% rownames(Paola_6_Seurat)]

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in data_list) {
  assign(x = i,value = AddModuleScore(object = get(i),features = list(gene_list),name = 'PC'))
  print(paste(i,'done!',sep = ' '))
  gc()
}

#plot
temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$re_anno %in% 'RG-neu',]
  meta_data <- meta_data[,c("re_anno","PC1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$group <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat'),"group"] <- 'organoid'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('PFC_RNA_Seurat'),"author"] <- 'Bhaduri. et al'
temp[temp$dataset %in% c('Hp_RNA_Seurat','Ho_RNA_Seurat'),"author"] <- 'Pollen. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'
temp[temp$dataset %in% c('Paola_3_Seurat','Paola_6_Seurat'),"author"] <- 'Velasco. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = data_list)
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Bhaduri. et al','Pollen. et al','Uzquiano. et al','Velasco. et al'))

p1 <- ggplot(data = temp[temp$group == 'primary',],mapping = aes(x = author,y = PC1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.25,1)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00','#39CCCC')) + 
  theme(aspect.ratio = 7/4.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('PC module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"PC1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p2 <- ggplot(data = temp[temp$group == 'organoid',],mapping = aes(x = author,y = PC1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.25,1)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#39CCCC','#0074D9','#2ECC40')) + 
  theme(aspect.ratio = 7/3.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('PC module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"PC1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_113_fig_230326/primary_organoid_PC_signature_boxplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

# SC signature ------------------------------------------------------------
#get gene list
gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat) & gene_list %in% rownames(PD_RNA_Seurat) & gene_list %in% rownames(PFC_RNA_Seurat) & gene_list %in% rownames(Hp_RNA_Seurat) & gene_list %in% rownames(Ho_RNA_Seurat) & gene_list %in% rownames(Ho_1.5_Seurat) & gene_list %in% rownames(Ho_2_Seurat) & gene_list %in% rownames(Ho_3_Seurat) & gene_list %in% rownames(Ho_4_Seurat) & gene_list %in% rownames(Paola_3_Seurat) & gene_list %in% rownames(Paola_6_Seurat)]

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat')
for (i in data_list) {
  assign(x = i,value = AddModuleScore(object = get(i),features = list(gene_list),name = 'SC'))
  print(paste(i,'done!',sep = ' '))
  gc()
}

#plot
temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$re_anno %in% 'RG-neu',]
  meta_data <- meta_data[,c("re_anno","SC1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$group <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','PFC_RNA_Seurat','Hp_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat','Paola_3_Seurat','Paola_6_Seurat'),"group"] <- 'organoid'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('PFC_RNA_Seurat'),"author"] <- 'Bhaduri. et al'
temp[temp$dataset %in% c('Hp_RNA_Seurat','Ho_RNA_Seurat'),"author"] <- 'Pollen. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'
temp[temp$dataset %in% c('Paola_3_Seurat','Paola_6_Seurat'),"author"] <- 'Velasco. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = data_list)
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Bhaduri. et al','Pollen. et al','Uzquiano. et al','Velasco. et al'))

p1 <- ggplot(data = temp[temp$group == 'primary',],mapping = aes(x = author,y = SC1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.25,1)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00','#39CCCC')) + 
  theme(aspect.ratio = 7/4.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('SC module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"SC1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p2 <- ggplot(data = temp[temp$group == 'organoid',],mapping = aes(x = author,y = SC1,fill = author,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + coord_cartesian(ylim = c(-0.25,1)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#39CCCC','#0074D9','#2ECC40')) + 
  theme(aspect.ratio = 7/3.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('SC module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"SC1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_113_fig_230326/primary_organoid_SC_signature_boxplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

# select classic and important RG marker for dot plot ---------------------
#set gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(RG_RNA_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(RG_RNA_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(RG_RNA_Seurat)]
gene_list <- getGOgeneSet(x = 'GO:0031012',OrgDb = org.Hs.eg.db,ont = 'CC',keytype = 'SYMBOL')
classic_marker <- c('SLC1A3','SOX9','SOX6','PAX6','VIM','PTN','HOPX','FAM107A','TNC','FZD2','SFRP1','CD9')
#ECM_gene <- c(HG_gene_list,PC_gene_list)[c(HG_gene_list,PC_gene_list,SC_gene_list) %in% gene_list]
SC_gene_list[SC_gene_list %in% gene_list]
PC_gene_list[PC_gene_list %in% gene_list]
HG_gene_list[HG_gene_list %in% gene_list]

#order gene list
HG_ECM_gene <- c('ANXA2','LGALS3','ITGA2','ANXA6','ANXA5','TIMP1','LRRC17','SPARC')
PC_ECM_gene <- c('ITGB5','BMP7','COL11A1','APOE','LRIG1','LRRTM3','HTRA1','F3','LRRC3B','BCAN','CLU')
SC_ECM_gene <- c('COL4A5','MFGE8','GPC4')
gene_list <- list(core = classic_marker,HG = HG_ECM_gene,PC = PC_ECM_gene,SC = SC_ECM_gene)

#plot
p <- DotPlot(object = RG_RNA_Seurat[,!(RG_RNA_Seurat$sample %in% c('Velasco. et al. 6m'))],assay = 'RNA',features = gene_list,group.by = 'dataset',cols = c('white','#00366CFF'),scale = FALSE)
p$data$avg.exp.scaled[which(p$data$avg.exp.scaled > 1.5)] <- 1.5
p <- p + guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(text = element_text(size = 10,family = 'sans'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  labs(title = '') + 
  scale_y_discrete(labels = c('Trevino. et al.','Polioudakis. et al.','Bhaduri. et al.','Pollen. et al. Primary','Pollen. et al. Organoid','Uzquiano. et al.','Velasco. et al. 3m')) + 
  theme(legend.position = 'bottom')

pdf(file = './res/step_113_fig_230326/Primary_Organoid_ECM_marker_dotplot.pdf',width = 12,height = 5)
p
dev.off()

# QC ----------------------------------------------------------------------

#QC matrix
QC_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_2_Seurat','Paola_3_Seurat'),FUN = function(i){
  temp <- get(i)
  temp <- data.frame(nCount_RNA = temp$nCount_RNA, nFeature_RNA = temp$nFeature_RNA,dataset = i)
  return(temp)
}))

#QC_matrix <- gather(data = QC_matrix,key = 'QC',value = 'value','nCount_RNA','nFeature_RNA')
QC_matrix$dataset <- factor(QC_matrix$dataset,levels = c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_2_Seurat','Paola_3_Seurat'))

#plot
p1 <- ggplot(data = QC_matrix,mapping = aes(x = dataset,y = nCount_RNA,fill = dataset)) + 
  geom_violin() + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank()) + 
  NoLegend()

p2 <- ggplot(data = QC_matrix,mapping = aes(x = dataset,y = nFeature_RNA,fill = dataset)) + 
  geom_violin() + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank()) + 
  NoLegend()

pdf(file = './res/step_113_fig_230326/why_give_up/dataset_QC_vlnplot.pdf',width = 8,height = 4)
p1+p2+plot_layout(ncol = 2)
dev.off()