#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Organoids marker gene expression evaluate                       ##
## Data: 2023.03.03                                                                ##
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
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

Hp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Hp_RNA_Seurat_230302.rds')
Mp_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Mp_RNA_Seurat_230302.rds')
Ho_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Ho_RNA_Seurat_230302.rds')
Co_RNA_Seurat <- readRDS(file = './data/public/Establishing_Cerebral_Organoids_as_Models_of_Human_Specific_Brain_Evolution/processed_data/Co_RNA_Seurat_230302.rds')

Ho_2_RNA_Seurat <- readRDS(file = './data/public/Organoid_single_cell_genomic_atlas_uncovers_human_specific_features_of_brain_development/processed_data/human_filted_Seurat_object.rds')
Co_2_RNA_Seurat <- readRDS(file = './data/public/Organoid_single_cell_genomic_atlas_uncovers_human_specific_features_of_brain_development/processed_data/Chimp_filted_Seurat_object.rds')

#reanno
Ho_2_RNA_Seurat$cell_type <- Ho_2_RNA_Seurat$PredCellType
Ho_2_RNA_Seurat@meta.data[Ho_2_RNA_Seurat$PredCellType == 'EN',"cell_type"] <- 'Ex'
Ho_2_RNA_Seurat@meta.data[Ho_2_RNA_Seurat$PredCellType == 'IN',"cell_type"] <- 'In'
Ho_2_RNA_Seurat@meta.data[Ho_2_RNA_Seurat$PredCellType == 'IPC',"cell_type"] <- 'IP'
Ho_2_RNA_Seurat@meta.data[Ho_2_RNA_Seurat$PredCellType %in% c('OPC','End','Per','VLMC','Mic'),"cell_type"] <- 'others'

Co_2_RNA_Seurat$cell_type <- Co_2_RNA_Seurat$PredCellType
Co_2_RNA_Seurat@meta.data[Co_2_RNA_Seurat$PredCellType == 'EN',"cell_type"] <- 'Ex'
Co_2_RNA_Seurat@meta.data[Co_2_RNA_Seurat$PredCellType == 'IN',"cell_type"] <- 'In'
Co_2_RNA_Seurat@meta.data[Co_2_RNA_Seurat$PredCellType == 'IPC',"cell_type"] <- 'IP'
Co_2_RNA_Seurat@meta.data[Co_2_RNA_Seurat$PredCellType == 'Microglia',"cell_type"] <- 'Mic'
Co_2_RNA_Seurat@meta.data[Co_2_RNA_Seurat$PredCellType %in% c('OPC','End','Per','VLMC','Mic'),"cell_type"] <- 'others'

Greenleaf_RNA_Seurat$cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('RG-1','RG-2'),"cell_type"] <- 'RG'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"cell_type"] <- 'Ex'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('InMGE','InCGE'),"cell_type"] <- 'In'
Greenleaf_RNA_Seurat@meta.data[Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('OPC','End','Per','VLMC','Mic'),"cell_type"] <- 'others'

Hp_RNA_Seurat@meta.data[Hp_RNA_Seurat$cell_type %in% c('OPC','End','Per','VLMC','Mic'),"cell_type"] <- 'others'

Mp_RNA_Seurat@meta.data[Mp_RNA_Seurat$cell_type %in% c('OPC','End','Per','VLMC','Mic'),"cell_type"] <- 'others'

Ho_RNA_Seurat@meta.data[Ho_RNA_Seurat$cell_type %in% c('OPC','End','Per','VLMC','Mic'),"cell_type"] <- 'others'

Co_RNA_Seurat@meta.data[Co_RNA_Seurat$cell_type %in% c('OPC','End','Per','VLMC','Mic'),"cell_type"] <- 'others'

#parameter
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230211.rds')
cell_type_list <- c('RG','IP','Ex','In','others')

#subset data
for (i in c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Mp_RNA_Seurat','Ho_RNA_Seurat','Ho_2_RNA_Seurat','Co_RNA_Seurat','Co_2_RNA_Seurat')) {
  temp <- get(x = i)
  temp <- temp[,temp$cell_type %in% cell_type_list]
  temp$cell_type <- factor(temp$cell_type,levels = cell_type_list)
  Idents(temp) <- 'cell_type'
  assign(x = i,value = temp)
}

# human primary and organoids ECM gene expression dotplot -----------------
gene_list <- c('SFRP1','TNC','COL4A5','GPC4','COL11A1','BMP7','F3','ITGB5','LGALS3','ITGA2','ANXA2')
scale.min <- 0
scale.max <- base::lapply(X = list(Greenleaf_RNA_Seurat,Hp_RNA_Seurat,Ho_RNA_Seurat,Ho_2_RNA_Seurat),FUN = function(x){
  temp <- DotPlot(object = x,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE)
  return(max(temp$data$pct.exp))
})
scale.max <- max(unlist(scale.max))

p1 <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
              cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max) + 
  coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = length(gene_list)/5,
        text = element_text(size = 10,family = 'sans')) + 
  labs(title = 'Greenleaf primary') + 
  NoLegend()

p2 <- DotPlot(object = Hp_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
              cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max) + 
  coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = length(gene_list)/5,
        text = element_text(size = 10,family = 'sans')) + 
  labs(title = 'Human primary') + 
  NoLegend() + theme(axis.text = element_blank())

p3 <- DotPlot(object = Ho_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
              cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max) + 
  coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = length(gene_list)/5,
        text = element_text(size = 10,family = 'sans')) + 
  labs(title = 'Human organoid') + 
  NoLegend() + theme(axis.text = element_blank())

p4 <- DotPlot(object = Ho_2_RNA_Seurat,assay = 'RNA',features = gene_list,group.by = 'cell_type',scale = TRUE,
              cols = c('lightgrey',col_param$species['human']),scale.min = scale.min,scale.max = scale.max) + 
  coord_flip() + RotatedAxis() + guides(color = guide_colorbar(title = "Average Expression")) + 
  guides(size = guide_legend(title = "Percent Expressed")) + 
  theme(panel.border = element_rect(fill = NA,color = 'black',linewidth = 0.8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio = length(gene_list)/5,
        text = element_text(size = 10,family = 'sans')) + 
  labs(title = 'Human organoid 2') + 
  NoLegend() + theme(axis.text = element_blank())

pdf(file = './res/step_106_fig_230303/human_primary_and_human_organoid_ECM_gene_expression_difference.pdf',width = 10,height = 5)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

#feature plot
data_list <- c('Greenleaf_RNA_Seurat','Hp_RNA_Seurat','Ho_RNA_Seurat','Ho_2_RNA_Seurat')
gene_list <- c('SFRP1','TNC','COL4A5','GPC4','COL11A1','BMP7','F3','ITGB5','LGALS3','ITGA2','ANXA2')

for (j in gene_list) {
  for (i in 1:length(data_list)) {
    assign(x = paste0('p',i),value = DimPlot(object = get(data_list[i]),pt.size = 0.001,group.by = 'cell_type',label = TRUE,repel = TRUE,raster = FALSE) + 
             theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + NoAxes() + NoLegend() + labs(title = data_list[i]))
    assign(x = paste0('p',i+length(data_list)),value = FeaturePlot(object = get(data_list[i]),features = j,pt.size = 0.001,cols = ArchRPalettes$whitePurple[-1],raster = FALSE) + 
             theme_bw() + theme(aspect.ratio = 1) + CenterTitle() + NoAxes() + NoLegend() + labs(title = j))
  }
  char <- paste0('./res/step_106_fig_230303/human_primary_and_human_organoid_ECM_',j,'_expression_featureplot.pdf')
  pdf(file = char,width = 12,height = 6)
  print(p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 4))
  dev.off()
}

gene_list <- c('EGFR2')
#ECM gene module score
gene_list <- getGOgeneSet(x = 'GO:0031012',OrgDb = org.Hs.eg.db,ont = 'CC',keytype = 'SYMBOL')
