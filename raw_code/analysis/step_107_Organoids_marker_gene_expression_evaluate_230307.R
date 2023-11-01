#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Organoids marker gene expression evaluate                       ##
## Data: 2023.03.07                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

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

#parameter
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230211.rds')

# RG-1 and RG-2 identity in organoid --------------------------------------

#re-anno cell type label
Ho_Barbara_Seurat$cell_type <- Ho_Barbara_Seurat$nowakowski_prediction
Ho_1.5_Seurat$cell_type <- Ho_1.5_Seurat$FinalName
Ho_2_Seurat$cell_type <- Ho_2_Seurat$FinalName
Ho_3_Seurat$cell_type <- Ho_3_Seurat$FinalName
Ho_4_Seurat$cell_type <- Ho_4_Seurat$FinalName

#plot
gene_list <- c('SLC1A3','PAX6','PTPRZ1','EGFR','OLIG2')
data_list <- c('Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')

for (j in gene_list) {
  for (i in 1:length(data_list)) {
    assign(x = paste0('p',i),value = DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'cell_type',label = TRUE,repel = TRUE,raster = FALSE) + 
             theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + NoAxes() + NoLegend() + labs(title = data_list[i]))
    assign(x = paste0('p',i+length(data_list)),value = FeaturePlot(object = get(data_list[i]),features = j,pt.size = 0.001,cols = ArchRPalettes$whitePurple[-1],raster = FALSE) + 
             theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + NoAxes() + NoLegend() + labs(title = j))
  }
  char <- paste0('./res/step_107_fig_230307/RG_neu_glia_identity_in_organoid/human_organoid_RG_',j,'_expression_featureplot.pdf')
  pdf(file = char,width = 18,height = 6)
  print(p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+plot_layout(ncol = 6))
  dev.off()
}

gene_list <- c('TRIB2','DUSP6','MAP3K1','PDGFRA','OLIG2','CHST7','EGFR')
data_list <- c('Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')

for (j in gene_list) {
  for (i in 1:length(data_list)) {
    assign(x = paste0('p',i),value = DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'cell_type',label = TRUE,repel = TRUE,raster = FALSE) + 
             theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + NoAxes() + NoLegend() + labs(title = data_list[i]))
    assign(x = paste0('p',i+length(data_list)),value = FeaturePlot(object = get(data_list[i]),features = j,pt.size = 0.001,cols = ArchRPalettes$whitePurple[-1],raster = FALSE) + 
             theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + NoAxes() + NoLegend() + labs(title = j))
  }
  char <- paste0('./res/step_107_fig_230307/RG_glia_identity_in_organoid/human_organoid_RG_',j,'_expression_featureplot.pdf')
  pdf(file = char,width = 18,height = 6)
  print(p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+plot_layout(ncol = 6))
  dev.off()
}

#consider RG-1 only

# re-annotate cell type label ---------------------------------------------
#Greenleaf data
DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'ReAnno_celltype',label = TRUE,repel = TRUE,cols = col_param$celltype) + theme_bw() + theme(aspect.ratio = 1)

Greenleaf_RNA_Seurat$cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('RG-1'),"cell_type"] <- 'RG-neu'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"cell_type"] <- 'Ex'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('InCGE','InMGE'),"cell_type"] <- 'In'
Greenleaf_RNA_Seurat@meta.data[!(Greenleaf_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

#Human organoid
DimPlot(object = Ho_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type %in% c('RG'),"cell_type"] <- 'RG-neu'
Ho_RNA_Seurat@meta.data[!(Ho_RNA_Seurat$cell_type %in% c('RG-neu','IP','Ex','In','Cycling')),"cell_type"] <- 'others'

DimPlot(object = Ho_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

#human barbara organoid
Ho_Barbara_Seurat$cell_type <- Ho_Barbara_Seurat$nowakowski_prediction
DimPlot(object = Ho_Barbara_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('RG'),"cell_type"] <- 'RG-neu'
Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('EN'),"cell_type"] <- 'Ex'
Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('IN'),"cell_type"] <- 'In'
Ho_Barbara_Seurat@meta.data[Ho_Barbara_Seurat$cell_type %in% c('IPC'),"cell_type"] <- 'IP'
Ho_Barbara_Seurat@meta.data[!(Ho_Barbara_Seurat$cell_type %in% c(c('RG-neu','IP','Ex','In','Cycling'))),"cell_type"] <- 'others'

DimPlot(object = Ho_Barbara_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

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

DimPlot(object = Hp_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

# subset data -------------------------------------------------------------
cell_type_list <- c('RG-neu','IP','Ex','In','others')
for (i in c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_Paola_Seurat')) {
  temp <- get(x = i)
  temp <- temp[,temp$cell_type %in% cell_type_list]
  temp$cell_type <- factor(temp$cell_type,levels = cell_type_list)
  Idents(temp) <- 'cell_type'
  assign(x = i,value = temp)
}

# ECM marker gene ---------------------------------------------------------
gene_list <- c('SFRP1','TNC','COL4A5','GPC4','COL11A1','BMP7','F3','ITGB5','LGALS3','ITGA2','ANXA2')
scale.min <- 0
scale.max <- base::lapply(X = list(Greenleaf_RNA_Seurat,Hp_RNA_Seurat,Ho_RNA_Seurat,Ho_Barbara_Seurat,Ho_Paola_Seurat),FUN = function(x){
  temp <- DotPlot(object = x,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE)
  return(max(temp$data$pct.exp))
})
scale.max <- max(unlist(scale.max))

data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_Paola_Seurat')
for (i in 1:length(data_list)) {
  p <- DotPlot(object = get(data_list[i]),assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
               cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max) + 
    coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
    guides(size = guide_legend(title = "Percent Expressed")) + 
    theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) + 
    theme(aspect.ratio = length(gene_list)/5,
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

pdf(file = './res/step_107_fig_230307/human_primary_organoid_ECM_gene_dotplot.pdf',width = 12,height = 5)
p1+p2+p3+p4+p5+plot_layout(ncol = 5)
dev.off()

# feature plot validate cell type identity --------------------------------
Ho_1.5_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_1.5_Seurat',colnames(Ho_1.5_Seurat),sep = '_'),"cell_type"]
Ho_2_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_2_Seurat',colnames(Ho_2_Seurat),sep = '_'),"cell_type"]
Ho_3_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_3_Seurat',colnames(Ho_3_Seurat),sep = '_'),"cell_type"]
Ho_4_Seurat$cell_type <- Ho_Paola_Seurat@meta.data[paste('Ho_4_Seurat',colnames(Ho_4_Seurat),sep = '_'),"cell_type"]

col_list <- col_param$celltype[c("RG","Cycling",'IP',"Ex-1","InMGE","Mic")]
names(col_list) <- c('RG-neu','Cycling','IP','Ex','In','others')

data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'cell_type',label = TRUE,repel = TRUE,raster = FALSE,cols = col_list) + 
    theme_bw() + theme(aspect.ratio = 1) + NoLegend() + CenterTitle() + 
    labs(title = data_list[i]) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  assign(x = paste0('p',i),value = p)
  
  #feature plot
  gene_list <- c('VIM','PAX6','EOMES','NEUROD2','GAD2')
  for (j in 1:length(gene_list)) {
    indi <- FALSE
    if(i %in% c(2,3)){
      indi <- TRUE
    }
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],pt.size = 0.001,order = indi,raster = FALSE,cols = ArchRPalettes$whitePurple[-1]) + 
      theme_bw() + theme(aspect.ratio = 1) + NoLegend() + CenterTitle() + 
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_blank())
    assign(x = paste0('p',i+j*length(data_list)),value = p)
  }
}

char <- paste0('p',1:48) %>% paste(collapse = '+') %>% paste0('+plot_layout(ncol = 8)')

pdf(file = './res/step_107_fig_230307/human_primary_organoid_MKG_featureplot.pdf',width = 24,height = 18)
eval(parse(text = char))
dev.off()

png(filename = './res/step_107_fig_230307/human_primary_organoid_MKG_featureplot.png',width = 2400,height = 1800)
eval(parse(text = char))
dev.off()

# ECM gene express is less in organoid ------------------------------------

#get extra cellular matrix GO gene
gene_list <- getGOgeneSet(x = 'GO:0031012',OrgDb = org.Hs.eg.db,ont = 'CC',keytype = 'SYMBOL')
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat) & gene_list %in% rownames(Hp_RNA_Seurat) & gene_list %in% rownames(Ho_Barbara_Seurat) & gene_list %in% rownames(Ho_Paola_Seurat)]

#add module score
data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in data_list) {
  assign(x = i,value = AddModuleScore(object = get(i),features = list(gene_list),name = 'ECM'))
  print(paste(i,'done!',sep = ' '))
  gc()
}

#add module score to paola data
Ho_Paola_Seurat@meta.data[paste('Ho_1.5_Seurat',colnames(Ho_1.5_Seurat),sep = '_'),'ECM1'] <- Ho_1.5_Seurat$ECM1
Ho_Paola_Seurat@meta.data[paste('Ho_2_Seurat',colnames(Ho_2_Seurat),sep = '_'),'ECM1'] <- Ho_2_Seurat$ECM1
Ho_Paola_Seurat@meta.data[paste('Ho_3_Seurat',colnames(Ho_3_Seurat),sep = '_'),'ECM1'] <- Ho_3_Seurat$ECM1
Ho_Paola_Seurat@meta.data[paste('Ho_4_Seurat',colnames(Ho_4_Seurat),sep = '_'),'ECM1'] <- Ho_4_Seurat$ECM1

#generate ggplot data
data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_Paola_Seurat')
cell_type_list <- c('RG-neu','IP','Ex','In')

temp <- do.call(what = rbind,args = base::lapply(X = data_list,FUN = function(i){
  meta_data <- get(i)@meta.data
  meta_data <- meta_data[meta_data$cell_type %in% cell_type_list,]
  meta_data <- meta_data[,c("cell_type","ECM1")]
  meta_data$dataset <- i
  return(meta_data)
}))

temp$cell_type <- factor(temp$cell_type,levels = cell_type_list)
temp$dataset <- factor(temp$dataset,levels = data_list)

pdf(file = './res/step_107_fig_230307/human_primary_organoid_ECM_module_score_boxplot.pdf',width = 8,height = 5.5)
ggplot(data = temp,mapping = aes(x = cell_type,y = ECM1,fill = dataset)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2) + 
  ylim(c(-0.2,0.2)) + theme_ArchR() + 
  theme(aspect.ratio = 0.6) + 
  scale_fill_manual(values = c('#32CD32','#87CEEB','#B57EDC','#FFF44F','#FFC0CB')) + 
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 8)) + 
  xlab('cell type') + ylab('ECM module score')
dev.off()

#only RG-neu
temp <- temp[temp$cell_type == 'RG-neu',]

pdf(file = './res/step_107_fig_230307/human_primary_organoid_ECM_module_score_boxplot_only_RG_neu.pdf',width = 5.5,height = 5.5)
ggplot(data = temp,mapping = aes(x = dataset,y = ECM1,fill = dataset)) + 
  geom_boxplot(outlier.alpha = 0,size = 0.2,width = 0.6) + 
  ylim(c(-0.25,0.25)) + theme_ArchR() + 
  theme(aspect.ratio = 2) + 
  scale_fill_manual(values = c('#32CD32','#87CEEB','#B57EDC','#FFF44F','#FFC0CB')) + 
  theme(legend.position = 'right',
        legend.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + CenterTitle() + 
  xlab('dataset') + ylab('ECM module score') + labs(title = 'RG-neu') + 
  stat_compare_means(comparisons = list(c(1,2),c(1,3),c(1,4),c(1,5)),
                     label = 'p.signif',method = 't.test',
                     label.y = c(0.07,0.09,0.11,0.13),tip.length = 0.005,vjust = 0.5) 
dev.off()

#expression featureplot
col_list <- col_param$celltype[c("RG","Cycling",'IP',"Ex-1","InMGE","Mic")]
names(col_list) <- c('RG-neu','Cycling','IP','Ex','In','others')

data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_1.5_Seurat','Ho_2_Seurat','Ho_3_Seurat','Ho_4_Seurat')
for (i in 1:length(data_list)) {
  p <- DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'cell_type',label = TRUE,repel = TRUE,raster = FALSE,cols = col_list) + 
    theme_bw() + theme(aspect.ratio = 1) + NoLegend() + CenterTitle() + 
    labs(title = data_list[i]) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  assign(x = paste0('p',i),value = p)
  
  #feature plot
  gene_list <- c('ANXA2','ITGA2','LGALS3','ITGB5','COL11A1')
  for (j in 1:length(gene_list)) {
    if(i %in% c(1,2)){
      indi <- TRUE
    }else{
      indi <- FALSE
    }
    p <- FeaturePlot(object = get(data_list[i]),features = gene_list[j],pt.size = 0.001,order = indi,raster = FALSE,cols = ArchRPalettes$whitePurple[-1]) + 
      theme_bw() + theme(aspect.ratio = 1) + NoLegend() + CenterTitle() + 
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_blank())
    assign(x = paste0('p',i+j*length(data_list)),value = p)
  }
}

char <- paste0('p',1:(8*6)) %>% paste(collapse = '+') %>% paste0('+plot_layout(ncol = 8)')

png(filename = './res/step_107_fig_230307/human_primary_organoid_ECM_gene_featureplot.png',width = 2400,height = 1800)
eval(parse(text = char))
dev.off()

#express percentage barplot
Greenleaf_RNA_Seurat$group <- paste(Greenleaf_RNA_Seurat$Sample.ID,Greenleaf_RNA_Seurat$cell_type,sep = '_')
Hp_RNA_Seurat$group <- paste(Hp_RNA_Seurat$DonorID,Hp_RNA_Seurat$cell_type,sep = '_')
Ho_RNA_Seurat$group <- paste(Ho_RNA_Seurat$DonorID,Ho_RNA_Seurat$cell_type,sep = '_')
Ho_Barbara_Seurat$group <- paste(Ho_Barbara_Seurat$sample,Ho_Barbara_Seurat$cell_type,sep = '_')
Ho_Paola_Seurat$group <- paste(Ho_Paola_Seurat$dataset,Ho_Paola_Seurat$cell_type,sep = '_')

data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_Barbara_Seurat','Ho_Paola_Seurat')
gene_list <- c('ANXA2','ITGA2','LGALS3','ITGB5','F3','BMP7','COL11A1','TNC')
cell_type_list <- c('RG-neu','IP','Ex','In','others')
temp <- base::do.call(what = rbind,args = base::lapply(X = 1:length(data_list),FUN = function(i){
  p <- DotPlot(object = get(data_list[i]),assay = 'RNA',features = gene_list,group.by = 'group',scale = TRUE,
               cols = c('lightgrey',col_param$species['human']),scale.min = 0,scale.max = 100)
  p <- p$data
  p$dataset <- data_list[i]
  return(p)
}))
temp <- temp[grep(pattern = 'RG-neu$',x = temp$id,fixed = FALSE),]

for (i in gene_list) {
  pct_matrix <- temp[temp$features.plot == i,]
  pct_matrix %>% dplyr::summarize(mean_pct = mean(pct.exp),sd_pct = sd(pct.exp),.by = dataset) -> df_bar
  df_bar$dataset <- factor(df_bar$dataset,levels = data_list)
  
  p <- ggplot(data = df_bar,mapping = aes(x = dataset,y = mean_pct,fill = dataset)) + 
    geom_bar(stat = 'identity',width = 0.6) + 
    scale_fill_manual(values = c('#32CD32','#87CEEB','#B57EDC','#FFF44F','#FFC0CB')) + 
    geom_errorbar(mapping = aes(ymin = mean_pct - sd_pct,ymax = mean_pct + sd_pct),width = 0.2,color = 'black',linewidth = 0.2) + 
    theme_ArchR() + theme(aspect.ratio = 2) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 8)) + 
    labs(title = i) + CenterTitle()
  
  if(i != 'TNC'){
    p <- p + NoLegend()
  }
  if(!(i %in% c('ANXA2','F3'))){
    p <- p + theme(axis.title.y = element_blank())
  }
  
  assign(x = paste('p',i,sep = '_'),value = p)
}

paste('p',gene_list,sep = '_') %>% paste(collapse = '+') %>% paste0('+plot_layout(ncol = 4)') -> char
pdf(file = './res/step_107_fig_230307/human_primary_organoid_ECM_gene_express_pct_barplot.pdf',width = 12,height = 9)
eval(parse(text = char))
dev.off()