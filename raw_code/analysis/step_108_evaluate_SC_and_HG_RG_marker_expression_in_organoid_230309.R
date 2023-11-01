#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: evaluate SC and HG RG marker expression in organoid             ##
## Data: 2023.03.09                                                                ##
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
library(OpenAI4R)
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize OpenAI
Auth_OpenAI(key = readLines(con = '/content/script/openai_API_key'))
chat <- Init_chat_session()

#initialize ArchR
addArchRThreads(threads = 5)


# load data ---------------------------------------------------------------
#Greenleaf data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

#human primary
Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')

#macaque multiome data
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

#macaque primary data
Mp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Mp_RNA_Seurat_230302.rds')

#human organoid
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Ho_Barbara_Seurat <- readRDS(file = './data/public/Inferring_and_perturbing_cell_fate_regulomes_in_human_brain_organoids/processed_data/Organoid_RNA_Seurat.rds')
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
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('RG-1'),"cell_type"] <- 'RG-neu'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"cell_type"] <- 'Ex'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('InCGE','InMGE'),"cell_type"] <- 'In'
Greenleaf_RNA_Seurat@meta.data[!(Greenleaf_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

#Human organoid
Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type %in% c('RG'),"cell_type"] <- 'RG-neu'
Ho_RNA_Seurat@meta.data[!(Ho_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

#human barbara organoid
Ho_Barbara_Seurat$cell_type <- Ho_Barbara_Seurat$nowakowski_prediction

Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('RG'),"cell_type"] <- 'RG-neu'
Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('EN'),"cell_type"] <- 'Ex'
Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('IN'),"cell_type"] <- 'In'
Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('IPC'),"cell_type"] <- 'IP'
Ho_Barbara_Seurat@meta.data[!(Ho_Barbara_Seurat$cell_type %in% c(c('RG-neu','IP','Ex','In','Cycling'))),"cell_type"] <- 'others'

#human  Paola organoid
Ho_Paola_Seurat$cell_type <- Ho_Paola_Seurat$FinalName
table(Ho_Paola_Seurat$cell_type)

Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('aRG','oRG'),"cell_type"] <- 'RG-neu'
Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('CFuPN','CPN','Newborn CFuPN','Newborn CPN','Newborn DL PN','Newborn PN','PN','Subcortical neurons','Subcortical progenitors'),"cell_type"] <- 'Ex'
Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('Immature IN'),"cell_type"] <- 'In'
Ho_Paola_Seurat@meta.data[Ho_Paola_Seurat$cell_type %in% c('IP'),"cell_type"] <- 'IP'
Ho_Paola_Seurat@meta.data[!(Ho_Paola_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

table(Ho_Paola_Seurat$cell_type)

#human primary
table(Hp_RNA_Seurat$cell_type)

Hp_RNA_Seurat@meta.data[Hp_RNA_Seurat$cell_type %in% c('RG'),"cell_type"] <- 'RG-neu'
Hp_RNA_Seurat@meta.data[!(Hp_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

#PD human
PD_RNA_Seurat$cell_type <- 'others'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('oRG','vRG'),'cell_type'] <- 'RG-neu'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('PgG2M','PgS'),'cell_type'] <- 'Cycling'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('IP'),'cell_type'] <- 'IP'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('ExDp1','ExDp2','ExM','ExM-U','ExN'),'cell_type'] <- 'Ex'
PD_RNA_Seurat@meta.data[PD_RNA_Seurat$Cluster %in% c('InCGE','InMGE'),'cell_type'] <- 'In'

# subset data -------------------------------------------------------------
cell_type_list <- c('RG-neu','IP','Ex','In','others')
for (i in c('Greenleaf_RNA_Seurat','PD_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_Paola_Seurat')) {
  temp <- get(x = i)
  temp <- NormalizeData(temp)
  temp <- temp[,temp$cell_type %in% cell_type_list]
  temp$cell_type <- factor(temp$cell_type,levels = cell_type_list)
  Idents(temp) <- 'cell_type'
  assign(x = i,value = temp)
  gc()
}

# try percentage change ---------------------------------------------------

#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Barbara_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Barbara_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Barbara_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

#calculate mean percentage
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"pct.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"pct.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"pct.exp"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"pct.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix <- gather(data = pct_matrix,key = 'group',value = 'pct','primary','organoid')
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))

#plot
pdf(file = './res/step_108_fig_230309/RG_neu_signature_gene_expression_percentage_between_primary_and_organoid_boxplot.pdf',width = 4,height = 5)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = pct,fill = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  scale_fill_manual(values = c('#FF5733','#00BFFF')) + 
  theme_ArchR() + 
  theme(aspect.ratio = 1.4,
        legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  scale_x_discrete(labels = c('Human-Specific','Primate-Specific','Species-Conserved')) + 
  theme(axis.title.x = element_blank()) + 
  ylab('Percent %') + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')
dev.off()

#try to show the difference between primary and organoid
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"pct.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"pct.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"pct.exp"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"pct.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix[pct_matrix$gene_group == 'HG_gene_list',"gene_group"] <- 'Human-Specific'
pct_matrix[pct_matrix$gene_group == 'PC_gene_list',"gene_group"] <- 'Primate-Specific'
pct_matrix[pct_matrix$gene_group == 'SC_gene_list',"gene_group"] <- 'Species-Conserved'
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('Human-Specific','Primate-Specific','Species-Conserved'))
pct_matrix$diff <- pct_matrix$primary - pct_matrix$organoid

pdf(file = './res/step_108_fig_230309/RG_neu_signature_gene_expression_percentage_difference_between_primary_and_organoid_boxplot.pdf',width = 5,height = 5.5)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = diff,fill = gene_group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  scale_fill_manual(values = as.character(col_param$species[c("human","macaque","mouse")])) + 
  theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red',linewidth = 0.2) + 
  stat_compare_means(comparisons = list(c(2,3),c(1,2),c(1,3)),
                     label = 'p.format',method = 't.test') + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('Percent difference %') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')
dev.off()

#try to show the fold change between primary and organoid
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"pct.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"pct.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"pct.exp"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"pct.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix[pct_matrix$gene_group == 'HG_gene_list',"gene_group"] <- 'Human-Specific'
pct_matrix[pct_matrix$gene_group == 'PC_gene_list',"gene_group"] <- 'Primate-Specific'
pct_matrix[pct_matrix$gene_group == 'SC_gene_list',"gene_group"] <- 'Species-Conserved'
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('Human-Specific','Primate-Specific','Species-Conserved'))
pct_matrix$fc <- pct_matrix$primary / pct_matrix$organoid

ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = fc,fill = gene_group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  ylim(c(0,3)) + 
  scale_fill_manual(values = as.character(col_param$species[c("human","macaque","mouse")])) + 
  theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red',linewidth = 0.2) + 
  stat_compare_means(comparisons = list(c(2,3),c(1,2),c(1,3)),
                     label = 'p.format',method = 't.test',
                     label.y = c(-1.4,-1.2,-1),tip.length = 0.001) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('Percent difference %') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')

## group by sample ---------------------------------------------------------
#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Barbara_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Barbara_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Barbara_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

#calculate mean percentage
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"pct.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"pct.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(c(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"pct.exp"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"pct.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"pct.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = rep(x = get(x),times = 6),gene_group = x,pct = c(primary,organoid),
                     group = c(rep(x = 'primary',times = length(primary)),rep(x = 'organoid',times = length(organoid))),
                     sample = c(rep('Trevino, et al',times = length(get(x))),rep('Polioudakis, et al',times = length(get(x))),rep('Pollen, et al',times = length(get(x))),
                                rep('Pollen, et al',times = length(get(x))),rep('Fleck, et al',times = length(get(x))),rep('Uzquiano, et al',times = length(get(x)))))
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))
pct_matrix$sample <- paste(pct_matrix$group,pct_matrix$sample,sep = '#')
pct_matrix$sample <- factor(pct_matrix$sample,levels = c('primary#Trevino, et al','primary#Polioudakis, et al','primary#Pollen, et al',
                                                         'organoid#Pollen, et al','organoid#Fleck, et al','organoid#Uzquiano, et al'))

#plot
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = pct,fill = sample)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  theme_ArchR()

# try expression change ---------------------------------------------------

## normalized expression ---------------------------------------------------

#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Barbara_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Barbara_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Barbara_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

#calculate mean percentage
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"avg.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix <- gather(data = pct_matrix,key = 'group',value = 'exp','primary','organoid')
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))

#plot
pdf(file = './res/step_108_fig_230309/RG_neu_signature_gene_expression_between_primary_and_organoid_boxplot.pdf',width = 4,height = 5)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = exp,fill = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  scale_fill_manual(values = c('#FF5733','#00BFFF')) + 
  ylim(c(0,3)) + 
  theme_ArchR() + 
  theme(aspect.ratio = 1.4,
        legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  scale_x_discrete(labels = c('Human-Specific','Primate-Specific','Species-Conserved')) + 
  theme(axis.title.x = element_blank()) + 
  ylab('Average expession') + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')
dev.off()

#try to show the difference between primary and organoid
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"avg.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix[pct_matrix$gene_group == 'HG_gene_list',"gene_group"] <- 'Human-Specific'
pct_matrix[pct_matrix$gene_group == 'PC_gene_list',"gene_group"] <- 'Primate-Specific'
pct_matrix[pct_matrix$gene_group == 'SC_gene_list',"gene_group"] <- 'Species-Conserved'
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('Human-Specific','Primate-Specific','Species-Conserved'))
pct_matrix$diff <- pct_matrix$primary - pct_matrix$organoid

pdf(file = './res/step_108_fig_230309/RG_neu_signature_gene_expression_difference_between_primary_and_organoid_boxplot.pdf',width = 5,height = 5.5)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = diff,fill = gene_group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  ylim(c(-1,2)) + 
  scale_fill_manual(values = as.character(col_param$species[c("human","macaque","mouse")])) + 
  theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red',linewidth = 0.2) + 
  stat_compare_means(comparisons = list(c(2,3),c(1,2),c(1,3)),
                     label = 'p.format',method = 't.test',
                     label.y = c(-2.3,-2.1,-1.9),tip.length = 0.001) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('expression difference') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')
dev.off()

#try to show the difference between primary and organoid
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"avg.exp"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"avg.exp"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix[pct_matrix$gene_group == 'HG_gene_list',"gene_group"] <- 'Human-Specific'
pct_matrix[pct_matrix$gene_group == 'PC_gene_list',"gene_group"] <- 'Primate-Specific'
pct_matrix[pct_matrix$gene_group == 'SC_gene_list',"gene_group"] <- 'Species-Conserved'
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('Human-Specific','Primate-Specific','Species-Conserved'))
pct_matrix$fc <- pct_matrix$primary / pct_matrix$organoid

ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = log2(fc),fill = gene_group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  ylim(c(-2,5)) + 
  scale_fill_manual(values = as.character(col_param$species[c("human","macaque","mouse")])) + 
  theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red',linewidth = 0.2) + 
  stat_compare_means(comparisons = list(c(2,3),c(1,2),c(1,3)),
                     label = 'p.format',method = 't.test',
                     label.y = c(3.2,3.6,4),tip.length = 0.001) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('log2(primary/organoid)') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')

## scaled expression ---------------------------------------------------

#gene list
HG_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
HG_gene_list <- HG_gene_list[HG_gene_list %in% rownames(Greenleaf_RNA_Seurat) & HG_gene_list %in% rownames(PD_RNA_Seurat) & HG_gene_list %in% rownames(Hp_RNA_Seurat) & HG_gene_list %in% rownames(Ho_RNA_Seurat) & HG_gene_list %in% rownames(Ho_Barbara_Seurat) & HG_gene_list %in% rownames(Ho_Paola_Seurat)]
SC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
SC_gene_list <- SC_gene_list[SC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & SC_gene_list %in% rownames(PD_RNA_Seurat) & SC_gene_list %in% rownames(Hp_RNA_Seurat) & SC_gene_list %in% rownames(Ho_RNA_Seurat) & SC_gene_list %in% rownames(Ho_Barbara_Seurat) & SC_gene_list %in% rownames(Ho_Paola_Seurat)]
PC_gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
PC_gene_list <- PC_gene_list[PC_gene_list %in% rownames(Greenleaf_RNA_Seurat) & PC_gene_list %in% rownames(PD_RNA_Seurat) & PC_gene_list %in% rownames(Hp_RNA_Seurat) & PC_gene_list %in% rownames(Ho_RNA_Seurat) & PC_gene_list %in% rownames(Ho_Barbara_Seurat) & PC_gene_list %in% rownames(Ho_Paola_Seurat)]

#calculate mean percentage
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp.scaled"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp.scaled"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"avg.exp.scaled"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp.scaled"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"avg.exp.scaled"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"avg.exp.scaled"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix$gene_group <- gsub(pattern = '_gene_list',replacement = '',x = pct_matrix$gene_group,fixed = TRUE)
pct_matrix <- gather(data = pct_matrix,key = 'group',value = 'exp','primary','organoid')
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('HG','PC','SC'))
pct_matrix$group <- factor(pct_matrix$group,levels = c('primary','organoid'))

#plot
pdf(file = './res/step_108_fig_230309/RG_neu_signature_scaled_gene_expression_between_primary_and_organoid_boxplot.pdf',width = 4,height = 5)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = exp,fill = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  scale_fill_manual(values = c('#FF5733','#00BFFF')) + 
  theme_ArchR() + 
  theme(aspect.ratio = 1.4,
        legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  scale_x_discrete(labels = c('Human-Specific','Primate-Specific','Species-Conserved')) + 
  theme(axis.title.x = element_blank()) + 
  ylab('Average expession') + 
  stat_compare_means(aes(group = group),
                     label = 'p.signif',method = 't.test') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')
dev.off()

#try to show the difference between primary and organoid
pct_matrix <- base::do.call(what = rbind,args = base::lapply(X = c('HG_gene_list','PC_gene_list','SC_gene_list'),FUN = function(x){
  p_Greenleaf <- DotPlot(object = Greenleaf_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_PD <- DotPlot(object = PD_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Hp <- DotPlot(object = Hp_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Ho <- DotPlot(object = Ho_RNA_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Barbara <- DotPlot(object = Ho_Barbara_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  p_Paola <- DotPlot(object = Ho_Paola_Seurat,features = get(x),group.by = 'cell_type',idents = cell_type_list)
  
  primary <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Greenleaf$data[p_Greenleaf$data$id == 'RG-neu' & p_Greenleaf$data$features.plot == y,"avg.exp.scaled"],p_PD$data[p_PD$data$id == 'RG-neu' & p_PD$data$features.plot == y,"avg.exp.scaled"],p_Hp$data[p_Hp$data$id == 'RG-neu' & p_Hp$data$features.plot == y,"avg.exp.scaled"]))
  }) %>% unlist
  
  organoid <- base::lapply(X = get(x),FUN = function(y){
    return(mean(p_Ho$data[p_Ho$data$id == 'RG-neu' & p_Ho$data$features.plot == y,"avg.exp.scaled"],p_Barbara$data[p_Barbara$data$id == 'RG-neu' & p_Barbara$data$features.plot == y,"avg.exp.scaled"],p_Paola$data[p_Paola$data$id == 'RG-neu' & p_Paola$data$features.plot == y,"avg.exp.scaled"]))
  }) %>% unlist
  
  temp <- data.frame(gene_name = get(x),gene_group = x,primary = primary,organoid = organoid)
  return(temp)
}))

pct_matrix[pct_matrix$gene_group == 'HG_gene_list',"gene_group"] <- 'Human-Specific'
pct_matrix[pct_matrix$gene_group == 'PC_gene_list',"gene_group"] <- 'Primate-Specific'
pct_matrix[pct_matrix$gene_group == 'SC_gene_list',"gene_group"] <- 'Species-Conserved'
pct_matrix$gene_group <- factor(pct_matrix$gene_group,levels = c('Human-Specific','Primate-Specific','Species-Conserved'))
pct_matrix$diff <- pct_matrix$primary - pct_matrix$organoid

pdf(file = './res/step_108_fig_230309/RG_neu_signature_scaled_gene_expression_difference_between_primary_and_organoid_boxplot.pdf',width = 5,height = 5.5)
ggplot(data = pct_matrix,mapping = aes(x = gene_group,y = diff,fill = gene_group)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.4) + 
  scale_fill_manual(values = as.character(col_param$species[c("human","macaque","mouse")])) + 
  theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red',linewidth = 0.2) + 
  stat_compare_means(comparisons = list(c(2,3),c(1,2),c(1,3)),
                     label = 'p.format',method = 't.test') + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8)) + 
  theme(axis.title.x = element_blank()) + 
  ylab('expression difference') + 
  theme(axis.title = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(size = 14,hjust = 0.5),
        axis.text.x = element_text(angle = 20,hjust = 0.5,vjust = 0.5)) + 
  labs(title = 'RG-neu Signature gene')
dev.off()