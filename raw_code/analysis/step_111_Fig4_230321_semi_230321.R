#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Fig4_230321_semi                                                ##
## Data: 2023.03.21                                                                ##
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

# ECM signature in human organoid -----------------------------------------

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')

#re annotate
table(Greenleaf_RNA_Seurat$ReAnno_celltype)
Greenleaf_RNA_Seurat$unify_cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype
Greenleaf_RNA_Seurat$unify_cell_type <- as.character(Greenleaf_RNA_Seurat$unify_cell_type)
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('RG-1'),"unify_cell_type"] <- 'RG-neu'

table(PD_RNA_Seurat$Cluster)
PD_RNA_Seurat$unify_cell_type <- PD_RNA_Seurat$Cluster
PD_RNA_Seurat$unify_cell_type <- as.character(PD_RNA_Seurat$unify_cell_type)
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('oRG','vRG'),"unify_cell_type"] <- 'RG-neu'

table(Hp_RNA_Seurat$cell_type)
Hp_RNA_Seurat$unify_cell_type <- Hp_RNA_Seurat$cell_type
Hp_RNA_Seurat$unify_cell_type <- as.character(Hp_RNA_Seurat$unify_cell_type)
Hp_RNA_Seurat@meta.data[Hp_RNA_Seurat$cell_type %in% c('RG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_RNA_Seurat$cell_type)
Ho_RNA_Seurat$unify_cell_type <- Ho_RNA_Seurat$cell_type
Ho_RNA_Seurat$unify_cell_type <- as.character(Ho_RNA_Seurat$unify_cell_type)
Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type %in% c('RG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_1.5_Seurat$FinalName)
Ho_1.5_Seurat$unify_cell_type <- Ho_1.5_Seurat$FinalName
Ho_1.5_Seurat$unify_cell_type <- as.character(Ho_1.5_Seurat$unify_cell_type)
Ho_1.5_Seurat@meta.data[Ho_1.5_Seurat$FinalName %in% c('aRG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_2_Seurat$FinalName)
Ho_2_Seurat$unify_cell_type <- Ho_2_Seurat$FinalName
Ho_2_Seurat$unify_cell_type <- as.character(Ho_2_Seurat$unify_cell_type)
Ho_2_Seurat@meta.data[Ho_2_Seurat$FinalName %in% c('aRG','oRG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_3_Seurat$FinalName)
Ho_3_Seurat$unify_cell_type <- Ho_3_Seurat$FinalName
Ho_3_Seurat$unify_cell_type <- as.character(Ho_3_Seurat$unify_cell_type)
Ho_3_Seurat@meta.data[Ho_3_Seurat$FinalName %in% c('aRG','oRG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_4_Seurat$FinalName)
Ho_4_Seurat$unify_cell_type <- Ho_4_Seurat$FinalName
Ho_4_Seurat$unify_cell_type <- as.character(Ho_4_Seurat$unify_cell_type)
Ho_4_Seurat@meta.data[Ho_4_Seurat$FinalName %in% c('aRG','oRG'),"unify_cell_type"] <- 'RG-neu'

#normalize data
for (i in c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')) {
  temp <- get(i)
  temp <- NormalizeData(temp)
  assign(x = i,value = temp)
}

#add module gene
gene_list <- getGOgeneSet(x = 'GO:0031012',OrgDb = org.Hs.eg.db,ont = 'CC',keytype = 'SYMBOL')
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat) & gene_list %in% rownames(PD_RNA_Seurat) & gene_list %in% rownames(Hp_RNA_Seurat) & gene_list %in% rownames(Ho_RNA_Seurat) & gene_list %in% rownames(Ho_1.5_Seurat) & gene_list %in% rownames(Ho_2_Seurat) & gene_list %in% rownames(Ho_3_Seurat) & gene_list %in% rownames(Ho_4_Seurat)]

#add module score
data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in data_list) {
  assign(x = i,value = AddModuleScore(object = get(i),features = list(gene_list),name = 'ECM'))
  print(paste(i,'done!',sep = ' '))
  gc()
}

## generate ggplot data type 1 ---------------------------------------------
temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$unify_cell_type %in% 'RG-neu',]
  meta_data <- meta_data[,c("unify_cell_type","ECM1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$group <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"group"] <- 'organoid'

temp[colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$Age == 'pcw16' & Greenleaf_RNA_Seurat$unify_cell_type == 'RG-neu'],"dataset"] <- 'Greenleaf_16'
temp[colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$Age == 'pcw20' & Greenleaf_RNA_Seurat$unify_cell_type == 'RG-neu'],"dataset"] <- 'Greenleaf_20'
temp[colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$Age == 'pcw21' & Greenleaf_RNA_Seurat$unify_cell_type == 'RG-neu'],"dataset"] <- 'Greenleaf_21'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_16','Greenleaf_20','Greenleaf_21'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('Hp_RNA_Seurat','Ho_RNA_Seurat'),"author"] <- 'Pollen. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = c('Greenleaf_16','Greenleaf_20','Greenleaf_21','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'))
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Pollen. et al','Uzquiano. et al'))

p1 <- ggplot(data = temp[temp$author == 'Trevino. et al',],mapping = aes(x = dataset,y = ECM1,fill = dataset,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  facet_grid(~ author,scales = 'free') + 
  theme_ArchR() + ylim(c(-0.1,0.15)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#32CD324C','#32CD3299','#32CD32','#87CEEB','#B57EDC','#FFF44F','#FFC0CB40','#FFC0CB80','#FFC0CBBF','#FFC0CB')) + 
  theme(aspect.ratio = 1.5) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('group') + ylab('ECM module score') + 
  theme(aspect.ratio = 9/3.2,
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'red') + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p2 <- ggplot(data = temp[temp$author == 'P.D. et al',],mapping = aes(x = dataset,y = ECM1,fill = dataset,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  facet_grid(~ author,scales = 'free') + 
  theme_ArchR() + ylim(c(-0.1,0.15)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#87CEEB')) + 
  theme(aspect.ratio = 1.5) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('group') + ylab('ECM module score') + 
  theme(aspect.ratio = 9/1.2,
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'red') + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p3 <- ggplot(data = temp[temp$author == 'Pollen. et al',],mapping = aes(x = dataset,y = ECM1,fill = dataset,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  facet_grid(~ author,scales = 'free') + 
  theme_ArchR() + ylim(c(-0.1,0.15)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#B57EDC','#FFF44F')) + 
  theme(aspect.ratio = 1.5) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('group') + ylab('ECM module score') + 
  theme(aspect.ratio = 9/2.2,
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'red') + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p4 <- ggplot(data = temp[temp$author == 'Uzquiano. et al',],mapping = aes(x = dataset,y = ECM1,fill = dataset,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  facet_grid(~ author,scales = 'free') + 
  theme_ArchR() + ylim(c(-0.1,0.15)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FFC0CB40','#FFC0CB80','#FFC0CBBF','#FFC0CB')) + 
  theme(aspect.ratio = 1.5) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('group') + ylab('ECM module score') + 
  theme(aspect.ratio = 9/4.2,
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'red') + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_111_fig_230321/human_primary_organoid_ECM_expression_boxplot_with_legand.pdf',width = 15,height = 8)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

p1 <- p1+NoLegend()
p2 <- p2+NoLegend()
p3 <- p3+NoLegend()
p4 <- p4+NoLegend()
pdf(file = './res/step_111_fig_230321/human_primary_organoid_ECM_expression_boxplot_without_legand.pdf',width = 10,height = 8)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

## generate ggplot data type 2 ---------------------------------------------
show_col(c('#FF41364C','#FF413699','#FF4136','#FF851B','#FFDC00','#39CCCC','#0074D940','#0074D980','#0074D9BF','#0074D9'))
temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$unify_cell_type %in% 'RG-neu',]
  meta_data <- meta_data[,c("unify_cell_type","ECM1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$group <- NA
temp[temp$dataset %in% c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat'),"group"] <- 'primary'
temp[temp$dataset %in% c('Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"group"] <- 'organoid'

temp[colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$Age == 'pcw16' & Greenleaf_RNA_Seurat$unify_cell_type == 'RG-neu'],"dataset"] <- 'Greenleaf_16'
temp[colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$Age == 'pcw20' & Greenleaf_RNA_Seurat$unify_cell_type == 'RG-neu'],"dataset"] <- 'Greenleaf_20'
temp[colnames(Greenleaf_RNA_Seurat)[Greenleaf_RNA_Seurat$Age == 'pcw21' & Greenleaf_RNA_Seurat$unify_cell_type == 'RG-neu'],"dataset"] <- 'Greenleaf_21'

temp$author <- NA
temp[temp$dataset %in% c('Greenleaf_16','Greenleaf_20','Greenleaf_21'),"author"] <- 'Trevino. et al'
temp[temp$dataset %in% c('PD_RNA_Seurat'),"author"] <- 'P.D. et al'
temp[temp$dataset %in% c('Hp_RNA_Seurat','Ho_RNA_Seurat'),"author"] <- 'Pollen. et al'
temp[temp$dataset %in% c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),"author"] <- 'Uzquiano. et al'

temp$group <- factor(temp$group,levels = c('primary','organoid'))
temp$dataset <- factor(temp$dataset,levels = c('Greenleaf_16','Greenleaf_20','Greenleaf_21','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'))
temp$author <- factor(temp$author,levels = c('Trevino. et al','P.D. et al','Pollen. et al','Uzquiano. et al'))

p1 <- ggplot(data = temp[temp$group == 'primary',],mapping = aes(x = dataset,y = ECM1,fill = dataset,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + ylim(c(-0.1,0.15)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#FF41364C','#FF413699','#FF4136','#FF851B','#FFDC00')) + 
  scale_x_discrete(labels = c('Trevino. et al. GW16','Trevino. et al. GW20','Trevino. et al. GW21','Polioudakis. et al. GW17-18','Pollen. et al. GW11-23')) + 
  theme(aspect.ratio = 7/5.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('ECM module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1)) + 
  geom_hline(yintercept = mean(temp[temp$group == 'primary',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

p2 <- ggplot(data = temp[temp$group == 'organoid',],mapping = aes(x = dataset,y = ECM1,fill = dataset,linetype = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.5) + 
  facet_grid(~ group,scales = 'free') + 
  theme_ArchR() + ylim(c(-0.1,0.15)) + 
  scale_linetype_manual(values = c('primary' = 'solid','organoid' = 'dashed')) + 
  scale_fill_manual(values = c('#39CCCC','#0074D940','#0074D980','#0074D9BF','#0074D9')) + 
  scale_x_discrete(labels = c('Pollen. et al. W5-15','Uzquiano. et al. M1.5','Uzquiano. et al. M2','Uzquiano. et al. M3','Uzquiano. et al. M4')) + 
  theme(aspect.ratio = 7/5.2) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('sample') + ylab('ECM module score') + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 20,hjust = 1,vjust = 1),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  geom_hline(yintercept = mean(temp[temp$group == 'organoid',"ECM1"]),linetype = 'dashed',linewidth = 0.5,color = 'grey')

pdf(file = './res/step_111_fig_230321/human_primary_organoid_ECM_expression_boxplot_with_legand_type_2.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

# featureplot legand ------------------------------------------------------

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

#color param
pal.GEX <- rev(paletteer_c("grDevices::Blues 3", 10))
pal.GEX[1] <- "#e0ecf4"
  
#plot
pdf(file = './res/step_111_fig_230321/featureplot_legand.pdf',width = 6,height = 6)
FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'ITGA2',pt.size = 0.01,cols = pal.GEX,order = TRUE,raster = FALSE) + 
  theme_bw() + NoAxes() + 
  theme(aspect.ratio = 1,
        legend.text = element_blank())
dev.off()

# ECM example gene expression pct in organoid -----------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
PD_RNA_Seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_1.5_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/1.5mo_harmonizedObj_060421.rds')
Ho_2_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/2mo_harmonizedObj_060421.rds')
Ho_3_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/3mo_harm_111120.rds')
Ho_4_Seurat <- readRDS(file = './data/public/Proper_acquisition_of_cell_class_identity_in_organoids_allows_definition_of_fate_specification_programs_of_the_human_cerebral_cortex/4mo_harm_060421.rds')

#re annotate
table(Greenleaf_RNA_Seurat$ReAnno_celltype)
Greenleaf_RNA_Seurat$unify_cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype
Greenleaf_RNA_Seurat$unify_cell_type <- as.character(Greenleaf_RNA_Seurat$unify_cell_type)
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('RG-1'),"unify_cell_type"] <- 'RG-neu'

table(PD_RNA_Seurat$Cluster)
PD_RNA_Seurat$unify_cell_type <- PD_RNA_Seurat$Cluster
PD_RNA_Seurat$unify_cell_type <- as.character(PD_RNA_Seurat$unify_cell_type)
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('oRG','vRG'),"unify_cell_type"] <- 'RG-neu'

table(Hp_RNA_Seurat$cell_type)
Hp_RNA_Seurat$unify_cell_type <- Hp_RNA_Seurat$cell_type
Hp_RNA_Seurat$unify_cell_type <- as.character(Hp_RNA_Seurat$unify_cell_type)
Hp_RNA_Seurat@meta.data[Hp_RNA_Seurat$cell_type %in% c('RG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_RNA_Seurat$cell_type)
Ho_RNA_Seurat$unify_cell_type <- Ho_RNA_Seurat$cell_type
Ho_RNA_Seurat$unify_cell_type <- as.character(Ho_RNA_Seurat$unify_cell_type)
Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type %in% c('RG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_1.5_Seurat$FinalName)
Ho_1.5_Seurat$unify_cell_type <- Ho_1.5_Seurat$FinalName
Ho_1.5_Seurat$unify_cell_type <- as.character(Ho_1.5_Seurat$unify_cell_type)
Ho_1.5_Seurat@meta.data[Ho_1.5_Seurat$FinalName %in% c('aRG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_2_Seurat$FinalName)
Ho_2_Seurat$unify_cell_type <- Ho_2_Seurat$FinalName
Ho_2_Seurat$unify_cell_type <- as.character(Ho_2_Seurat$unify_cell_type)
Ho_2_Seurat@meta.data[Ho_2_Seurat$FinalName %in% c('aRG','oRG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_3_Seurat$FinalName)
Ho_3_Seurat$unify_cell_type <- Ho_3_Seurat$FinalName
Ho_3_Seurat$unify_cell_type <- as.character(Ho_3_Seurat$unify_cell_type)
Ho_3_Seurat@meta.data[Ho_3_Seurat$FinalName %in% c('aRG','oRG'),"unify_cell_type"] <- 'RG-neu'

table(Ho_4_Seurat$FinalName)
Ho_4_Seurat$unify_cell_type <- Ho_4_Seurat$FinalName
Ho_4_Seurat$unify_cell_type <- as.character(Ho_4_Seurat$unify_cell_type)
Ho_4_Seurat@meta.data[Ho_4_Seurat$FinalName %in% c('aRG','oRG'),"unify_cell_type"] <- 'RG-neu'

Ho_1.5_Seurat$Age <- 1.5
Ho_2_Seurat$Age <- 2
Ho_3_Seurat$Age <- 3
Ho_4_Seurat$Age <- 4

#merge polar data
meta_data <- do.call(what = rbind,args = base::lapply(X = c('Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat'),FUN = function(x){
  temp <- get(x)@meta.data[,c('dataset','percent.mito','percent.ribo','S.Score','G2M.Score','Phase','CC.Difference','FinalName','unify_cell_type','Age')]
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

#normalize data
for (i in c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Paola_Seurat')) {
  temp <- get(i)
  temp <- NormalizeData(temp)
  assign(x = i,value = temp)
}

gc()

#group
Greenleaf_RNA_Seurat$group <- paste(Greenleaf_RNA_Seurat$Sample.ID,Greenleaf_RNA_Seurat$Age,Greenleaf_RNA_Seurat$unify_cell_type,sep = '_')
PD_RNA_Seurat$group <- paste(PD_RNA_Seurat$Donor,PD_RNA_Seurat$Gestation_week,PD_RNA_Seurat$unify_cell_type,sep = '_')
Hp_RNA_Seurat$group <- paste(Hp_RNA_Seurat$DonorID,Hp_RNA_Seurat$Age,Hp_RNA_Seurat$unify_cell_type,sep = '_')
Ho_RNA_Seurat$group <- paste(Ho_RNA_Seurat$DonorID,Ho_RNA_Seurat$Age,Ho_RNA_Seurat$unify_cell_type,sep = '_')
Ho_Paola_Seurat$group <- paste(Ho_Paola_Seurat$dataset,Ho_Paola_Seurat$Age,Ho_Paola_Seurat$unify_cell_type,sep = '_')

data_list <- c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Paola_Seurat')
gene_list <- c('ANXA2','ITGA2','LGALS3','ITGB5','F3','BMP7','COL11A1','GPC4','COL4A5','TNC','SFRP1')

temp <- base::do.call(what = rbind,args = base::lapply(X = 1:length(data_list),FUN = function(i){
  p <- DotPlot(object = get(data_list[i]),assay = 'RNA',features = gene_list,group.by = 'group',scale = TRUE,
               cols = c('lightgrey','blue'),scale.min = 0,scale.max = 100)
  p <- p$data
  p$dataset <- data_list[i]
  return(p)
}))
temp <- temp[grep(pattern = '_RG-neu$',x = temp$id,fixed = FALSE),]

for (i in gene_list) {
  pct_matrix <- temp[temp$features.plot == i,]
  pct_matrix %>% dplyr::summarize(mean_pct = mean(pct.exp),sd_pct = sd(pct.exp),.by = dataset) -> df_bar
  df_bar$dataset <- factor(df_bar$dataset,levels = data_list)
  
  p <- ggplot(data = df_bar,mapping = aes(x = dataset,y = mean_pct,fill = dataset)) + 
    geom_bar(stat = 'identity',width = 0.6) + 
    scale_fill_manual(values = c('#FF4136','#FF851B','#FFDC00','#39CCCC','#0074D9')) + 
    geom_errorbar(mapping = aes(ymin = mean_pct - 0.5*sd_pct,ymax = mean_pct + 0.5*sd_pct),width = 0.2,color = 'black',linewidth = 0.2) + 
    theme_ArchR() + theme(aspect.ratio = 1.5) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 8)) + 
    labs(title = i) + CenterTitle() + 
    theme(axis.title = element_blank()) + 
    NoLegend()
  
  assign(x = paste('p',i,sep = '_'),value = p)
}

p_ITGA2 + p_LGALS3 + p_ITGB5 + p_BMP7 + plot_layout(ncol = 4)
