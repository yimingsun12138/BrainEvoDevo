#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: fig_221228_semi                                                 ##
## Data: 2022.12.28                                                                ##
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
addArchRThreads(threads = 3)

# set working dir ---------------------------------------------------------
setwd('/home/sunym/temp/')

# macaque multiome correspondance -----------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
macaque_multiome_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')
color_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#dimplot
p1 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = macaque_multiome_ArchR@cellColData,
                 group.by = 'cell_type',label = TRUE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = 'cell type annotation')

p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = macaque_multiome_ArchR@cellColData,
                 group.by = 'predictedGroup',label = TRUE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = 'cell type prediction')

pdf(file = '/content/data/sunym/project/Brain/res/step_95_fig_221228/macaque_multiome_cell_type_dimplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = macaque_multiome_ArchR@cellColData,
                 group.by = 'cell_type',label = FALSE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = '')

p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = macaque_multiome_ArchR@cellColData,
                 group.by = 'predictedGroup',label = FALSE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = '')

png(filename = '/content/data/sunym/project/Brain/res/step_95_fig_221228/macaque_multiome_cell_type_dimplot.png',width = 1200,height = 600)
p1+p2+plot_layout(ncol = 2)
dev.off()

#correspondance
ori <- macaque_multiome_ArchR$cell_type
prd <- macaque_multiome_ArchR$predictedGroup
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = names(color_param$celltype))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(color_param$celltype))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_viridis() + theme(aspect.ratio = 1) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('cell type annotation') + ylab('predicted cell type') + labs(title = 'multiome correspondance')

pdf(file = '/content/data/sunym/project/Brain/res/step_95_fig_221228/macaque_multiome_cell_type_confusion_heatmap.pdf',width = 5,height = 5)
p
dev.off()

#gene score gene expression correlation
gene_score_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneScoreMatrix',verbose = TRUE)
gene_express_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneExpressionMatrix',verbose = TRUE)

temp <- gene_score_matrix@assays@data$GeneScoreMatrix
idx <- which(!duplicated(gene_score_matrix@elementMetadata$name))
temp <- temp[idx,]
rownames(temp) <- gene_score_matrix@elementMetadata$name[idx]
temp <- CreateSeuratObject(counts = temp,project = 'temp',assay = 'RNA',meta.data = as.data.frame(gene_score_matrix@colData),min.cells = 0,min.features = 0)
gene_score_matrix <- AverageExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
gene_score_matrix <- log1p(gene_score_matrix$RNA*100)

temp <- gene_express_matrix@assays@data$GeneExpressionMatrix
idx <- which(!duplicated(gene_express_matrix@elementMetadata$name))
temp <- temp[idx,]
rownames(temp) <- gene_express_matrix@elementMetadata$name[idx]
temp <- CreateSeuratObject(counts = temp,project = 'temp',assay = 'RNA',meta.data = as.data.frame(gene_express_matrix@colData),min.cells = 0,min.features = 0)
gene_express_matrix <- AverageExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
gene_express_matrix <- log1p(gene_express_matrix$RNA*100)

gene_list <- dplyr::intersect(x = rownames(gene_score_matrix),y = rownames(gene_express_matrix))
gene_list <- dplyr::intersect(x = gene_list,y = VariableFeatures(object = macaque_multiome_Seurat))
gene_score_matrix <- gene_score_matrix[gene_list,]
gene_express_matrix <- gene_express_matrix[gene_list,]

cor_matrix <- cor(gene_score_matrix,gene_express_matrix)

#heatmap plot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% colnames(cor_matrix)]
cor_matrix <- cor_matrix[cell_type_list,cell_type_list]

col_func <- colorRamp2(breaks = c(0.25,0.75),colors = ArchRPalettes$whiteRed)
p <- Heatmap(matrix = cor_matrix,cluster_columns = FALSE,cluster_rows = FALSE,
             column_split = factor(c('NPC','NPC','NPC','NPC','Ex','Ex','Ex','Ex','In','In','Others','Others','Others','Others','Others'),
                                   levels = c('NPC','Ex','In','Others')),
             row_split = factor(c('NPC','NPC','NPC','NPC','Ex','Ex','Ex','Ex','In','In','Others','Others','Others','Others','Others'),
                                levels = c('NPC','Ex','In','Others')),
             border = TRUE,col = col_func,height = unit(3.5,'inches'),width = unit(3.5,'inches'),
             column_title = 'gene express',row_title = 'gene score',name = 'correlation')

pdf(file = '/content/data/sunym/project/Brain/res/step_95_fig_221228/macaque_multiome_cell_type_gene_express_gene_score_cor_heatmap.pdf',width = 5.5,height = 5.5)
p
dev.off()

# mouse multiome correspondance -------------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
mouse_multiome_ArchR <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_ArchR_221009/')
color_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_221212.rds')

#dimplot
p1 <- my_dimplot(embedding = mouse_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = mouse_multiome_ArchR@cellColData,
                 group.by = 'Gex_macaque_cell_type',label = TRUE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = 'cell type annotation')

p2 <- my_dimplot(embedding = mouse_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = mouse_multiome_ArchR@cellColData,
                 group.by = 'predictedGroup',label = TRUE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = 'cell type prediction')

pdf(file = '/content/data/sunym/project/Brain/res/step_95_fig_221228/mouse_multiome_cell_type_dimplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- my_dimplot(embedding = mouse_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = mouse_multiome_ArchR@cellColData,
                 group.by = 'Gex_macaque_cell_type',label = FALSE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = '')

p2 <- my_dimplot(embedding = mouse_multiome_ArchR@embeddings$UMAP$df,
                 meta_data = mouse_multiome_ArchR@cellColData,
                 group.by = 'predictedGroup',label = FALSE,repel = TRUE,
                 cols = color_param$celltype,pt.size = 0.1) + 
  theme_cowplot() + theme(aspect.ratio = 1) + NoLegend() + NoAxes() + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(title = '')

png(filename = '/content/data/sunym/project/Brain/res/step_95_fig_221228/mouse_multiome_cell_type_dimplot.png',width = 1200,height = 600)
p1+p2+plot_layout(ncol = 2)
dev.off()

#correspondance
ori <- mouse_multiome_ArchR$Gex_macaque_cell_type
prd <- mouse_multiome_ArchR$predictedGroup
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = names(color_param$celltype))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(color_param$celltype))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_viridis() + theme(aspect.ratio = 1) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('cell type annotation') + ylab('predicted cell type') + labs(title = 'multiome correspondance')

pdf(file = '/content/data/sunym/project/Brain/res/step_95_fig_221228/mouse_multiome_cell_type_confusion_heatmap.pdf',width = 5,height = 5)
p
dev.off()

#gene score gene expression correlation
gene_score_matrix <- getMatrixFromProject(ArchRProj = mouse_multiome_ArchR,useMatrix = 'GeneScoreMatrix',verbose = TRUE)
gene_express_matrix <- getMatrixFromProject(ArchRProj = mouse_multiome_ArchR,useMatrix = 'GeneExpressionMatrix',verbose = TRUE)

temp <- gene_score_matrix@assays@data$GeneScoreMatrix
idx <- which(!duplicated(gene_score_matrix@elementMetadata$name))
temp <- temp[idx,]
rownames(temp) <- gene_score_matrix@elementMetadata$name[idx]
temp <- CreateSeuratObject(counts = temp,project = 'temp',assay = 'RNA',meta.data = as.data.frame(gene_score_matrix@colData),min.cells = 0,min.features = 0)
gene_score_matrix <- AverageExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'Gex_macaque_cell_type',slot = 'counts',verbose = TRUE)
gene_score_matrix <- log1p(gene_score_matrix$RNA*100)

temp <- gene_express_matrix@assays@data$GeneExpressionMatrix
idx <- which(!duplicated(gene_express_matrix@elementMetadata$name))
temp <- temp[idx,]
rownames(temp) <- gene_express_matrix@elementMetadata$name[idx]
temp <- CreateSeuratObject(counts = temp,project = 'temp',assay = 'RNA',meta.data = as.data.frame(gene_express_matrix@colData),min.cells = 0,min.features = 0)
gene_express_matrix <- AverageExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'Gex_macaque_cell_type',slot = 'counts',verbose = TRUE)
gene_express_matrix <- log1p(gene_express_matrix$RNA*100)

gene_list <- dplyr::intersect(x = rownames(gene_score_matrix),y = rownames(gene_express_matrix))
gene_list <- dplyr::intersect(x = gene_list,y = VariableFeatures(object = mouse_multiome_Seurat))
gene_score_matrix <- gene_score_matrix[gene_list,]
gene_express_matrix <- gene_express_matrix[gene_list,]

cor_matrix <- cor(gene_score_matrix,gene_express_matrix)

#heatmap plot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% colnames(cor_matrix)]
cor_matrix <- cor_matrix[cell_type_list,cell_type_list]

col_func <- colorRamp2(breaks = c(0.2,0.8),colors = ArchRPalettes$whiteRed)
p <- Heatmap(matrix = cor_matrix,cluster_columns = FALSE,cluster_rows = FALSE,
             column_split = factor(c('NPC','NPC','NPC','Ex','Ex','Ex','Ex','In','In','Others','Others','Others','Others'),
                                   levels = c('NPC','Ex','In','Others')),
             row_split = factor(c('NPC','NPC','NPC','Ex','Ex','Ex','Ex','In','In','Others','Others','Others','Others'),
                                levels = c('NPC','Ex','In','Others')),
             border = TRUE,col = col_func,height = unit(3.5,'inches'),width = unit(3.5,'inches'),
             column_title = 'gene express',row_title = 'gene score',name = 'correlation')

pdf(file = '/content/data/sunym/project/Brain/res/step_95_fig_221228/mouse_multiome_cell_type_gene_express_gene_score_cor_heatmap.pdf',width = 5.5,height = 5.5)
p
dev.off()