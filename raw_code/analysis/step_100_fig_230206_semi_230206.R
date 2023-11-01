#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: fig_230206_semi                                                 ##
## Data: 2023.02.06                                                                ##
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

# donor distribution of scATAC-seq data -----------------------------------

#load data
macaque_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_ATAC_ArchR_220912/')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#show the data
p1 <- my_dimplot(embedding = macaque_ATAC_ArchR@embeddings$harmonyUMAP$df,meta_data = macaque_ATAC_ArchR@cellColData,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_bw() + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_ATAC_ArchR@embeddings$harmonyUMAP$df,meta_data = macaque_ATAC_ArchR@cellColData,group.by = 'predictedGroup',label = TRUE,repel = TRUE) + 
  theme_bw() + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#donor distbution
donor_proportion_matrix <- My_Cell_Proportion(meta_data = as.data.frame(macaque_ATAC_ArchR@cellColData),group.by = 'donor',split.by = 'predictedGroup')
donor_list <- names(col_param$sample)
cell_type_list <- names(col_param$celltype)

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$predictedGroup <- factor(donor_proportion_matrix$predictedGroup,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = predictedGroup)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = col_param$sample) + 
  xlab('Cell type') + 
  scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0),name = 'Percentage') + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(linewidth = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(linewidth = 0.8),aspect.ratio = 0.4*15/14,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  guides(fill=guide_legend(title="Donors")) + 
  scale_x_discrete(labels = c('RG-neu','RG-gli','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','VLMC','Mic'))

pdf(file = './res/step_100_fig_230206/macaque_ATAC_ArchR_cell_type_donor_distribution_barplot_no_coord_flip.pdf',width = 6,height = 4)
p
dev.off()

# mouse multiome modality annotation consistency --------------------------
#set working dir
setwd('/home/sunym/temp/')

#load data
mouse_multiome_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_ArchR_221009/')
mouse_multiome_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
mouse_meta_data <- readRDS(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_annotation.rds')

col_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#re-anno
mouse_multiome_ArchR$ReAnno_celltype <- mouse_meta_data[rownames(mouse_multiome_ArchR@cellColData),"ReAnno_celltype"]
mouse_multiome_Seurat$ReAnno_celltype <- mouse_meta_data[rownames(mouse_multiome_Seurat@meta.data),"ReAnno_celltype"]
mouse_multiome_ArchR$predictedGroup <- mouse_meta_data[mouse_multiome_ArchR$predictedCell,"ReAnno_celltype"]

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "ReAnno_celltype", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#dimplot
p <- my_dimplot(embedding = mouse_multiome_ArchR@embeddings$UMAP$df,meta_data = mouse_multiome_ArchR@cellColData,group.by = 'predictedGroup',
                label = FALSE,pt.size = 0.001,raster = FALSE) + 
  scale_color_manual(values = col_param$celltype) + NoAxes() + 
  theme(plot.title = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1,
        legend.position = 'none',
        plot.margin = margin(0,0,0,0, "cm"))

ggsave(filename = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_and_multiome_RNA/mouse_multiome_multiome_RNA_predicted_anno_dimplot.png',plot = p,width = 3.5,height = 3.5)

#confusion heatmap
ori <- mouse_multiome_ArchR$ReAnno_celltype
prd <- mouse_multiome_ArchR$predictedGroup
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic'))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic'))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_gradient(low = '#F1F1F1',high = 'black') + 
  theme(aspect.ratio = 1) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'mouse multiome') + 
  scale_x_discrete(labels = c('RG-neu','Cyc-RG','Cyc-IP','IP','ExN','ExM','ExUp','ExDp','In','End','Per','VLMC','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','Cyc-RG','Cyc-IP','IP','ExN','ExM','ExUp','ExDp','In','End','Per','VLMC','Mic'))

pdf(file = '/content/data/sunym/project/Brain/res/step_100_fig_230206/mouse_multiome_annotation_and_multiome_RNA_predicted_annotation_confusion_heatmap.pdf',width = 5,height = 5)
p
dev.off()