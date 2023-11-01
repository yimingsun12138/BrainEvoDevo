#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: fig_230104_semi                                                 ##
## Data: 2023.01.04                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.1.3/')
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
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# QC nFrags TSS -----------------------------------------------------------

# set working dir
setwd('/home/sunym/temp/')

#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011')
macaque_ATAC_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_ATAC_ArchR_220912')
mouse_multiome_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_ArchR_221009')
col_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#macaque multiome
temp <- macaque_multiome_ArchR@cellColData
temp <- as.data.frame(temp)
temp$donor <- factor(temp$donor,levels = c('A50A','A50B','A82A','A82B'))
col_value <- col_param$sample

p1 <- ggplot(temp,aes(x = donor,y = nFrags,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,80000)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  scale_color_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.line = element_line(color = 'black',linewidth = 1)) + 
  ylab('# nFrags') + labs(title = 'macaque multiome')

p4 <- ggplot(temp,aes(x = donor,y = TSSEnrichment,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,40)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  scale_color_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.ticks = element_blank()) + 
  ylab('# TSS enrichment')

#macaque ATAC
temp <- macaque_ATAC_ArchR@cellColData
temp <- as.data.frame(temp)
temp$donor <- factor(temp$donor,levels = c('A50A','A82A','A82B','A84B','A84C','A68A','A68B'))
col_value <- col_param$sample

p2 <- ggplot(temp,aes(x = donor,y = nFrags,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,80000)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A82A','A82B','A84B','A84C','A68A','A68B')]) + 
  scale_color_manual(values = col_value[c('A50A','A82A','A82B','A84B','A84C','A68A','A68B')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.text.y = element_blank()) + 
  ylab('') + labs(title = 'macaque scATAC')

p5 <- ggplot(temp,aes(x = donor,y = TSSEnrichment,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,40)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A82A','A82B','A84B','A84C','A68A','A68B')]) + 
  scale_color_manual(values = col_value[c('A50A','A82A','A82B','A84B','A84C','A68A','A68B')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()) + 
  ylab('')

#mouse multiome
temp <- mouse_multiome_ArchR@cellColData
temp <- as.data.frame(temp)
temp$donor <- factor(temp$donor,levels = names(col_param$sample.mouse))
col_value <- col_param$sample.mouse

p3 <- ggplot(temp,aes(x = donor,y = nFrags,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,100000)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.line = element_line(color = 'black',linewidth = 1)) + 
  ylab('# nFrags') + labs(title = 'mouse multiome')

p6 <- ggplot(temp,aes(x = donor,y = TSSEnrichment,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,30)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.ticks = element_blank()) + 
  ylab('# TSS enrichment')

#plot with no lengend
pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_ATAC_QC_violin_plot_no_legend.pdf',width = 12,height = 4)
p1+p2+p4+p5+plot_layout(ncol = 2)
dev.off()

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_ATAC_QC_violin_plot_no_legend.pdf',width = 8,height = 4)
p3+p6+plot_layout(ncol = 1)
dev.off()

#plot with legend
temp <- macaque_multiome_ArchR@cellColData
temp <- as.data.frame(temp)
temp$donor <- factor(temp$donor,levels = c('A50A','A50B','A82A','A82B'))
col_value <- col_param$sample

p1 <- ggplot(temp,aes(x = donor,y = nFrags,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,80000)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  scale_color_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.line = element_line(color = 'black',linewidth = 1)) + 
  ylab('# nFrags') + labs(title = 'macaque multiome')

p4 <- ggplot(temp,aes(x = donor,y = TSSEnrichment,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,40)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  scale_color_manual(values = col_value[c('A50A','A50B','A82A','A82B')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.ticks = element_blank()) + 
  ylab('# TSS enrichment')

temp <- macaque_ATAC_ArchR@cellColData
temp <- as.data.frame(temp)
temp$donor <- factor(temp$donor,levels = names(col_param$sample))
col_value <- col_param$sample

p2 <- ggplot(temp,aes(x = donor,y = nFrags,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,80000)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value,drop = FALSE) + 
  scale_color_manual(values = col_value,drop = FALSE) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3.5) + 
  theme(legend.position = 'right',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.text.y = element_blank()) + 
  ylab('') + labs(title = 'macaque scATAC')

p5 <- ggplot(temp,aes(x = donor,y = TSSEnrichment,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,40)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()) + 
  ylab('')

temp <- mouse_multiome_ArchR@cellColData
temp <- as.data.frame(temp)
temp$donor <- factor(temp$donor,levels = names(col_param$sample.mouse))
col_value <- col_param$sample.mouse

p3 <- ggplot(temp,aes(x = donor,y = nFrags,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,100000)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3) + 
  theme(legend.position = 'right',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.line = element_line(color = 'black',linewidth = 1)) + 
  ylab('# nFrags') + labs(title = 'mouse multiome')

p6 <- ggplot(temp,aes(x = donor,y = TSSEnrichment,fill = donor,color = donor)) + 
  geom_violin() + 
  ylim(c(0,30)) + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',linewidth = 1),
        axis.ticks = element_blank()) + 
  ylab('# TSS enrichment')

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_ATAC_QC_violin_plot_with_legend.pdf',width = 12,height = 4)
p1+p2+p4+p5+plot_layout(ncol = 2)
dev.off()

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_ATAC_QC_violin_plot_with_legend.pdf',width = 8,height = 4)
p3+p6+plot_layout(ncol = 1)
dev.off()

# cell type donor distribution --------------------------------------------

#set working dir
setwd(dir = '/home/sunym/temp/')

#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011')
macaque_ATAC_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_ATAC_ArchR_220912')
mouse_multiome_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_ArchR_221009')

macaque_RNA_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_RNA_Seurat_220803.rds')

col_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

mouse_meta_data <- readRDS(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_annotation.rds')

#macaque multiome
donor_proportion_matrix <- My_Cell_Proportion(meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),group.by = 'donor',split.by = 'cell_type')
donor_list <- names(col_param$sample)
cell_type_list <- rev(names(col_param$celltype))

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$cell_type <- factor(donor_proportion_matrix$cell_type,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = cell_type)) +
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
        line = element_line(linewidth = 0.8),aspect.ratio = 2.5) + 
  guides(fill=guide_legend(title="Donors")) + 
  coord_flip()

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_multiome_cell_type_donor_distribution_barplot_coord_flip.pdf',width = 4,height = 5)
p
dev.off()

donor_proportion_matrix <- My_Cell_Proportion(meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),group.by = 'donor',split.by = 'cell_type')
donor_list <- names(col_param$sample)
cell_type_list <- names(col_param$celltype)

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$cell_type <- factor(donor_proportion_matrix$cell_type,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = cell_type)) +
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
        line = element_line(linewidth = 0.8),aspect.ratio = 0.4,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  guides(fill=guide_legend(title="Donors"))

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_multiome_cell_type_donor_distribution_barplot_no_coord_flip.pdf',width = 6,height = 4)
p
dev.off()

#macaque RNA
donor_proportion_matrix <- My_Cell_Proportion(meta_data = macaque_RNA_Seurat@meta.data,group.by = 'donor',split.by = 'cell_type')
donor_list <- names(col_param$sample)
cell_type_list <- rev(names(col_param$celltype))

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$cell_type <- factor(donor_proportion_matrix$cell_type,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = cell_type)) +
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
        line = element_line(linewidth = 0.8),aspect.ratio = 2.5) + 
  guides(fill=guide_legend(title="Donors")) + 
  coord_flip()

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_RNA_cell_type_donor_distribution_barplot_coord_flip.pdf',width = 4,height = 5)
p
dev.off()

donor_proportion_matrix <- My_Cell_Proportion(meta_data = macaque_RNA_Seurat@meta.data,group.by = 'donor',split.by = 'cell_type')
donor_list <- names(col_param$sample)
cell_type_list <- names(col_param$celltype)

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$cell_type <- factor(donor_proportion_matrix$cell_type,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = cell_type)) +
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
        line = element_line(linewidth = 0.8),aspect.ratio = 0.4,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  guides(fill=guide_legend(title="Donors"))

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_RNA_cell_type_donor_distribution_barplot_no_coord_flip.pdf',width = 6,height = 4)
p
dev.off()

#macaque ATAC
donor_proportion_matrix <- My_Cell_Proportion(meta_data = as.data.frame(macaque_ATAC_ArchR@cellColData),group.by = 'donor',split.by = 'cell_type')
donor_list <- names(col_param$sample)
cell_type_list <- rev(c('RG','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End/Mic'))

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$cell_type <- factor(donor_proportion_matrix$cell_type,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = cell_type)) +
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
        line = element_line(linewidth = 0.8),aspect.ratio = 5/3) + 
  guides(fill=guide_legend(title="Donors")) + 
  coord_flip()

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_ATAC_cell_type_donor_distribution_barplot_coord_flip.pdf',width = 4,height = 3.33)
p
dev.off()

donor_proportion_matrix <- My_Cell_Proportion(meta_data = as.data.frame(macaque_ATAC_ArchR@cellColData),group.by = 'donor',split.by = 'cell_type')
donor_list <- names(col_param$sample)
cell_type_list <- c('RG','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End/Mic')

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$cell_type <- factor(donor_proportion_matrix$cell_type,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = cell_type)) +
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
        line = element_line(linewidth = 0.8),aspect.ratio = 0.6,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  guides(fill=guide_legend(title="Donors"))

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_ATAC_cell_type_donor_distribution_barplot_no_coord_flip.pdf',width = 4,height = 4)
p
dev.off()

#mouse multiome
donor_proportion_matrix <- My_Cell_Proportion(meta_data = mouse_meta_data,group.by = 'donor',split.by = 'ReAnno_celltype')
donor_list <- names(col_param$sample.mouse)
cell_type_list <- rev(c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic'))

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$ReAnno_celltype <- factor(donor_proportion_matrix$ReAnno_celltype,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = ReAnno_celltype)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = col_param$sample.mouse) + 
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
        line = element_line(linewidth = 0.8),aspect.ratio = 2.5/15*13) + 
  guides(fill=guide_legend(title="Donors")) + 
  coord_flip()

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_cell_type_donor_distribution_barplot_coord_flip.pdf',width = 4,height = 4.33)
p
dev.off()

donor_proportion_matrix <- My_Cell_Proportion(meta_data = mouse_meta_data,group.by = 'donor',split.by = 'ReAnno_celltype')
donor_list <- names(col_param$sample.mouse)
cell_type_list <- c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic')

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = donor_list)
donor_proportion_matrix$ReAnno_celltype <- factor(donor_proportion_matrix$ReAnno_celltype,levels = cell_type_list)

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = ReAnno_celltype)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = col_param$sample.mouse) + 
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
        line = element_line(linewidth = 0.8),aspect.ratio = 15/13/2.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  guides(fill=guide_legend(title="Donors"))

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_cell_type_donor_distribution_barplot_no_coord_flip.pdf',width = 5.2,height = 4)
p
dev.off()

# mouse marker dotplot ----------------------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

mouse_meta_data <- readRDS(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_annotation.rds')

col_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#re-anno
table(rownames(mouse_meta_data) %in% rownames(mouse_multiome_Seurat@meta.data))
mouse_meta_data <- mouse_meta_data[rownames(mouse_multiome_Seurat@meta.data),]
mouse_multiome_Seurat$ReAnno_celltype <- mouse_meta_data$ReAnno_celltype

donor_seq <- c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic')

dotplot_matrix <- my_dotplot(mouse_multiome_Seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(RG=c('Sox2','Pax6','Hes5'),
                                             Cyc = c('Top2a','Mki67'),
                                             IP = c('Eomes','Btg2','Neurog2'),
                                             Ex = c('Neurod2','Neurod6','Nrp1','Satb2','Cux1','Ldb2','Fezf2','Bcl11b','Tle4','Foxp2'),
                                             In = c('Gad2','Dlx5','Lhx6','Adarb2'),
                                             End = c('Cldn5','Pecam1'),
                                             Per = c('Pdgfrb'),
                                             VLMC = c('Col1a1','Lum'),
                                             Mic = c('Cx3cr1','P2ry12','C1qb')),
                             group.by = 'ReAnno_celltype', cols = c('#3B4992FF','white','#EE0000FF'),
                             return_data_plot = TRUE)
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = donor_seq)

#plot
pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_Seurat_marker_dotplot.pdf',width = 13,height = 6)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom') + 
  xlab('') + ylab('')
dev.off()

# ATAC predicted annotation -----------------------------------------------

## macaque multiome and snRNA --------------------------------------------------------
#set working dir
setwd('/home/sunym/temp/')

#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_multiome_ArchR_221011')
macaque_RNA_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_RNA_Seurat_220803.rds')

col_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#constrained integration
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cols = col_param$celltype)

group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('End','Mic','OPC','Per','VLMC')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C1','C2','C10')]
  ),
  npc = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('Cycling','IP','RG-1','RG-2')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C11','C12','C14')]
  ),
  Ex_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C3','C4','C5','C6','C7','C8','C9','C13')]
  ),
  In_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C15','C16','C17','C18','C19','C20','C21')]
  )
)

macaque_multiome_ArchR <- addImputeWeights(macaque_multiome_ArchR)
macaque_multiome_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = macaque_RNA_Seurat,
  addToArrow = TRUE, 
  groupList = group_list,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#dimplot
p <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = macaque_multiome_ArchR@cellColData,group.by = 'predictedGroup',
                label = FALSE,pt.size = 0.001,raster = FALSE) + 
  scale_color_manual(values = col_param$celltype) + NoAxes() + 
  theme(plot.title = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1,
        legend.position = 'none',
        plot.margin = margin(0,0,0,0, "cm"))

ggsave(filename = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_multiome_and_snRNA/macaque_multiome_snRNA_predicted_anno_dimplot.png',plot = p,width = 3.5,height = 3.5)

#confusion heatmap
ori <- macaque_multiome_ArchR$cell_type
prd <- macaque_multiome_ArchR$predictedGroup
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = names(col_param$celltype))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(col_param$celltype))

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
  xlab('multiome annotation') + ylab('snRNA predicted annotation') + labs(title = 'macaque multiome')

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_multiome_and_snRNA/macaque_multiome_snRNA_predicted_anno_confusion_heatmap.pdf',width = 5,height = 5)
p
dev.off()

## macaque multiome and multiome RNA --------------------------------------------------------
#set working dir
setwd('/home/sunym/temp/')

#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_multiome_ArchR_221011')
macaque_multiome_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

col_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#constrained integration
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cols = col_param$celltype)

group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('End','Mic','OPC','Per','VLMC')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C1','C2','C10')]
  ),
  npc = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('Cycling','IP','RG-1','RG-2')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C11','C12','C14')]
  ),
  Ex_neuron = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C3','C4','C5','C6','C7','C8','C9','C13')]
  ),
  In_neuron = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C15','C16','C17','C18','C19','C20','C21')]
  )
)

macaque_multiome_ArchR <- addImputeWeights(macaque_multiome_ArchR)
macaque_multiome_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = macaque_multiome_Seurat,
  addToArrow = TRUE, 
  groupList = group_list,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)
my_send_sms('predict done!')

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#dimplot
p <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = macaque_multiome_ArchR@cellColData,group.by = 'predictedGroup',
                label = FALSE,pt.size = 0.001,raster = FALSE) + 
  scale_color_manual(values = col_param$celltype) + NoAxes() + 
  theme(plot.title = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1,
        legend.position = 'none',
        plot.margin = margin(0,0,0,0, "cm"))

ggsave(filename = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_multiome_and_multiome_RNA/macaque_multiome_multiome_RNA_predicted_anno_dimplot.png',plot = p,width = 3.5,height = 3.5)

#confusion heatmap
ori <- macaque_multiome_ArchR$cell_type
prd <- macaque_multiome_ArchR$predictedGroup
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = names(col_param$celltype))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(col_param$celltype))

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
  xlab('multiome annotation') + ylab('multiome RNA predicted annotation') + labs(title = 'macaque multiome')

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_multiome_and_multiome_RNA/macaque_multiome_multiome_RNA_predicted_anno_confusion_heatmap.pdf',width = 5,height = 5)
p
dev.off()

## macaque ATAC and snRNA --------------------------------------------------------
#set working dir
setwd('/home/sunym/temp/')

#load data
macaque_ATAC_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_ATAC_ArchR_220912')
macaque_RNA_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_RNA_Seurat_220803.rds')

col_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#constrained integration
p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "Sample", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p1+p2+p3+plot_layout(ncol = 3)

DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cols = col_param$celltype)

group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('End','Mic','OPC','Per','VLMC')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C18','C19')]
  ),
  npc = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('Cycling','IP','RG-1','RG-2')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C20','C21')]
  ),
  Ex_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C1','C2','C3','C4','C5','C6','C7','C8')]
  ),
  In_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C9','C10','C11','C12','C13','C14','C15','C16','C17')]
  )
)

macaque_ATAC_ArchR <- addImputeWeights(macaque_ATAC_ArchR)
macaque_ATAC_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_ATAC_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = macaque_RNA_Seurat,
  addToArrow = TRUE, 
  groupList = group_list,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)
my_send_sms('predict done!')

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "harmonyUMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")

p1+p2+p3+plot_layout(ncol = 3)

#dimplot
p <- my_dimplot(embedding = macaque_ATAC_ArchR@embeddings$UMAP$df,meta_data = macaque_ATAC_ArchR@cellColData,group.by = 'predictedGroup',
                label = FALSE,pt.size = 0.001,raster = FALSE) + 
  scale_color_manual(values = col_param$celltype) + NoAxes() + 
  theme(plot.title = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1,
        legend.position = 'none',
        plot.margin = margin(0,0,0,0, "cm"))

ggsave(filename = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_ATAC_and_snRNA/macaque_ATAC_snRNA_predicted_anno_dimplot.png',plot = p,width = 3.5,height = 3.5)

#confusion heatmap
ori <- macaque_ATAC_ArchR$cell_type
prd <- macaque_ATAC_ArchR$predictedGroup
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = c('RG','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End/Mic'))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(col_param$celltype))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_gradient(low = '#F1F1F1',high = 'black') + 
  theme(aspect.ratio = 1.5) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('scATAC annotation') + ylab('snRNA predicted annotation') + labs(title = 'macaque scATAC')

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/macaque_ATAC_and_snRNA/macaque_ATAC_snRNA_predicted_anno_confusion_heatmap.pdf',width = 5,height = 5)
p
dev.off()

## mouse multiome and multiome RNA --------------------------------------------------------
#set working dir
setwd('/home/sunym/temp/')

#load data
mouse_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/mouse_multiome_ArchR_221009')
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
  xlab('multiome annotation') + ylab('multiome RNA predicted annotation') + labs(title = 'mouse multiome')

pdf(file = '/content/data/sunym/project/Brain/res/step_96_fig_230104/mouse_multiome_and_multiome_RNA/mouse_multiome_multiome_RNA_predicted_anno_confusion_heatmap.pdf',width = 5,height = 5)
p
dev.off()