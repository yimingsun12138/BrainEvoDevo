#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque snRNA_seq quality control fig                           ##
## Data: 2022.09.15                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

# notice:
#donor sequence: A50A/A50B (E80), A82A/A82B (E90), A84B/A84C (E84), A68A/A68B (E92)

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#general setting
setwd('/data/User/sunym/project/Brain/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# project macaque E100 public data ----------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')
macaque_DFC_Seurat <- readRDS(file = './data/public/Spatiotemporal_transcriptomic_divergence_across_human_and_macaque_brain_development/processed_data/macaque_DFC_Seurat_220915.rds')

#shared variable feature
temp <- macaque_RNA_Seurat
temp <- FindVariableFeatures(object = temp,assay = 'RNA',selection.method = 'vst',nfeatures = 10000)
gene_list <- intersect(x = VariableFeatures(temp),y = rownames(macaque_DFC_Seurat@assays$RNA@counts))

macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 29,resolution = 1,group.by = 'cell_type',label = TRUE)
macaque_DFC_Seurat <- my_process_seurat(object = macaque_DFC_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#projection
temp <- projectMatrix_SeuratUMAP(X_scaled = macaque_DFC_Seurat@assays$RNA@scale.data,object = macaque_RNA_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
macaque_DFC_Seurat@reductions$projected_PCA <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
macaque_DFC_Seurat@reductions$projected_UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')

DimPlot(object = macaque_DFC_Seurat,reduction = 'projected_UMAP',label = TRUE,repel = TRUE,group.by = 'cell_type')

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_RNA_Seurat,query = macaque_DFC_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'RNA',l2.norm = FALSE,dims = 1:29,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_RNA_Seurat$cell_type,l2.norm = FALSE,dims = 1:29,verbose = TRUE)
macaque_DFC_Seurat <- AddMetaData(object = macaque_DFC_Seurat,metadata = predictions)
macaque_DFC_Seurat$new_cell_type <- macaque_DFC_Seurat$predicted.id
macaque_DFC_Seurat$new_cell_type_score <- macaque_DFC_Seurat$prediction.score.max

DimPlot(object = macaque_DFC_Seurat,reduction = 'projected_UMAP',label = TRUE,repel = TRUE,group.by = 'new_cell_type')

#macaque integrated data
macaque_integrated_Seurat <- intersect(x = rownames(macaque_RNA_Seurat@assays$RNA@counts),y = rownames(macaque_DFC_Seurat@assays$RNA@counts))
macaque_integrated_Seurat <- cbind(macaque_RNA_Seurat@assays$RNA@counts[macaque_integrated_Seurat,],macaque_DFC_Seurat@assays$RNA@counts[macaque_integrated_Seurat,])
macaque_integrated_Seurat <- CreateSeuratObject(counts = macaque_integrated_Seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)

#add meta data
macaque_integrated_Seurat$species <- 'macaque'
macaque_integrated_Seurat$dataset <- NA
macaque_integrated_Seurat$Age <- NA
macaque_integrated_Seurat$cell_type <- NA
macaque_integrated_Seurat$ori_cell_type <- NA

macaque_integrated_Seurat@meta.data[colnames(macaque_DFC_Seurat),"dataset"] <- 'public'
macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_Seurat),"dataset"] <- 'snRNA'
macaque_integrated_Seurat@meta.data[colnames(macaque_DFC_Seurat),"cell_type"] <- macaque_DFC_Seurat$new_cell_type
macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_Seurat),"cell_type"] <- macaque_RNA_Seurat$cell_type
macaque_integrated_Seurat@meta.data[colnames(macaque_DFC_Seurat),"ori_cell_type"] <- as.character(macaque_DFC_Seurat$cell_type)
macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_Seurat),"ori_cell_type"] <- macaque_RNA_Seurat$cell_type

macaque_RNA_Seurat$Age <- NA
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A50A','A50B'),"Age"] <- 'E80'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A82A','A82B'),"Age"] <- 'E90'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A84B','A84C'),"Age"] <- 'E84'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A68A','A68B'),"Age"] <- 'E92'

macaque_integrated_Seurat@meta.data[colnames(macaque_DFC_Seurat),"Age"] <- 'E110'
macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_Seurat),"Age"] <- macaque_RNA_Seurat$Age
macaque_integrated_Seurat@meta.data[1:3,]

#add pca
temp_1 <- macaque_RNA_Seurat@reductions$PCA@cell.embeddings
colnames(temp_1) <- paste('PC',1:50,sep = '_')
temp_2 <- macaque_DFC_Seurat@reductions$projected_PCA@cell.embeddings
colnames(temp_2) <- paste('PC',1:50,sep = '_')
temp <- rbind(temp_1,temp_2)
macaque_integrated_Seurat@reductions$pca <- CreateDimReducObject(embeddings = temp[colnames(macaque_integrated_Seurat),],assay = 'RNA',key = 'PC_')

temp_1 <- macaque_RNA_Seurat@reductions$umap@cell.embeddings
colnames(temp_1) <- paste('UMAP',1:2,sep = '_')
temp_2 <- macaque_DFC_Seurat@reductions$projected_UMAP@cell.embeddings
colnames(temp_2) <- paste('UMAP',1:2,sep = '_')
temp <- rbind(temp_1,temp_2)
macaque_integrated_Seurat@reductions$umap <- CreateDimReducObject(embeddings = temp[colnames(macaque_integrated_Seurat),],assay = 'RNA',key = 'UMAP_')

#display cell type
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype
cell_type_seq <- c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End','Per','VLMC','Mic')

x_lower <- min(macaque_integrated_Seurat@reductions$umap@cell.embeddings[,1])
x_upper <- max(macaque_integrated_Seurat@reductions$umap@cell.embeddings[,1])
y_lower <- min(macaque_integrated_Seurat@reductions$umap@cell.embeddings[,2])
y_upper <- max(macaque_integrated_Seurat@reductions$umap@cell.embeddings[,2])

p1 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[cell_type_seq])

p2 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'cell_type',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[cell_type_seq])

p3 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'cell_type',label = FALSE) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[cell_type_seq])

pdf(file = './res/step_58_fig_220915/macaque_integrated_Seurat_by_PCA_projection_cell_type_dimplot.pdf',width = 8,height = 8)
p1
p2
p3
dev.off()

#donor display
macaque_integrated_Seurat$donor_cell_type <- macaque_integrated_Seurat$cell_type
macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_Seurat),"donor_cell_type"] <- 'own'
temp_col <- append(temp_col,c('own' = '#e9e9e9'))

p1 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'donor_cell_type',label = TRUE,repel = TRUE,
                 order = cell_type_seq) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[c(cell_type_seq,'own')])

p2 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'donor_cell_type',label = FALSE,
                 order = cell_type_seq) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[c(cell_type_seq,'own')])

p3 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'donor_cell_type',label = FALSE,
                 order = cell_type_seq) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[c(cell_type_seq,'own')])

pdf(file = './res/step_58_fig_220915/macaque_integrated_Seurat_by_PCA_projection_cell_type_dimplot_on_own_data.pdf',width = 8,height = 8)
p1
p2
p3
dev.off()

#ori cell type display
macaque_integrated_Seurat$donor_ori_cell_type <- macaque_integrated_Seurat$ori_cell_type
macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_Seurat),"donor_ori_cell_type"] <- 'own'
temp_col <- c('Astro' = '#faa818','Blood' = '#41a30d','Endo' = '#fbdf72','ExN' = '#367d7d','IntN' = '#d33502',
              'Microglia' = '#6ebcbc','Oligo' = '#37526d','OPC' = '#916848','Peri' = '#f5b390','own' = '#e9e9e9')

p1 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'donor_ori_cell_type',label = TRUE,repel = TRUE,
                 order = c('IntN','Astro','Blood','Endo','ExN','Microglia','Oligo','OPC','Peri')) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col)

p2 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'donor_ori_cell_type',label = FALSE,
                 order = c('IntN','Astro','Blood','Endo','ExN','Microglia','Oligo','OPC','Peri')) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col)

p3 <- my_dimplot(embedding = macaque_integrated_Seurat@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(macaque_integrated_Seurat@meta.data),
                 group.by = 'donor_ori_cell_type',label = FALSE,
                 order = c('IntN','Astro','Blood','Endo','ExN','Microglia','Oligo','OPC','Peri')) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col)

pdf(file = './res/step_58_fig_220915/macaque_integrated_Seurat_by_PCA_projection_cell_type_dimplot_on_own_data.pdf',width = 8,height = 8)
p1
p2
p3
dev.off()

#macaque plot
temp <- macaque_integrated_Seurat[,macaque_integrated_Seurat$dataset == 'snRNA']

col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype
cell_type_seq <- c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End','Per','VLMC','Mic')

x_lower <- min(temp@reductions$umap@cell.embeddings[,1])
x_upper <- max(temp@reductions$umap@cell.embeddings[,1])
y_lower <- min(temp@reductions$umap@cell.embeddings[,2])
y_upper <- max(temp@reductions$umap@cell.embeddings[,2])

p1 <- my_dimplot(embedding = temp@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(temp@meta.data),
                 group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[cell_type_seq])

p2 <- my_dimplot(embedding = temp@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(temp@meta.data),
                 group.by = 'cell_type',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[cell_type_seq])

p3 <- my_dimplot(embedding = temp@reductions$umap@cell.embeddings,
                 meta_data = as.data.frame(temp@meta.data),
                 group.by = 'cell_type',label = FALSE) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[cell_type_seq])

pdf(file = './res/step_58_fig_220915/macaque_RNA_Seurat_by_PCA_projection_cell_type_dimplot.pdf',width = 8,height = 8)
p1
p2
p3
dev.off()

#confusion heatmap
scibet::Confusion_heatmap(ori = macaque_DFC_Seurat$cell_type,prd = macaque_DFC_Seurat$new_cell_type)

ori <- macaque_DFC_Seurat$cell_type
prd <- macaque_DFC_Seurat$new_cell_type
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End','Per','VLMC','Mic'))
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = c('Astro','ExN','IntN','OPC','Oligo','Endo','Peri','Blood','Microglia'))
p <- cross.validation.filt %>% ggplot(aes(ori, prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + scale_fill_viridis() + 
  theme(aspect.ratio = 15/9)

pdf(file = './res/step_58_fig_220915/macaque_integrated_Seurat_by_PCA_projection_cell_type_confusionheatmap.pdf',width = 5,height = 4)
p
dev.off()