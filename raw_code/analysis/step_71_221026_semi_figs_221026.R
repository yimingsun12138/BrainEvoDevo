#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: 221026_semi_figs                                                ##
## Data: 2022.10.26                                                                ##
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
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(SuperExactTest)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# marker gene enrichment exact test ---------------------------------------
#laod data
RG.1_marker <- readRDS(file = './res/step_71_fig_221026/ExciatoryRGIP-FindAllMarker-FilterBasedPerct_RG-1_expanded-mkgeneset_strict-humanspecific.rds')
IP_marker <- readRDS(file = './res/step_71_fig_221026/ExciatoryRGIP-FindAllMarker-FilterBasedPerct_IP_expanded-mkgeneset_strict-humanspecific.rds')
Ex.4_marker <- readRDS(file = './res/step_71_fig_221026/ExciatoryRGIP-FindAllMarker-FilterBasedPerct_Ex-4_expanded-mkgeneset_strict-humanspecific.rds')

#all genes
human_data <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_data <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_human_symbol_220802.rds')
mouse_data <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_human_symbol_221009.rds')

gene_list <- SuperExactTest::intersect(rownames(human_data),rownames(macaque_data),rownames(mouse_data))


## super test --------------------------------------------------------------

#RG.1
marker_list <- RG.1_marker
res <- supertest(x = marker_list,n = length(gene_list))

pdf(file = './res/step_71_fig_221026/RG.1_marker_intersect_supertest.pdf',width = 3.5,height = 4)
plot(res, Layout="landscape", degree=2:3,sort.by="size",margin=c(0.5,5,1,2))
dev.off()

temp <- data.frame(set = c('HR','HM','RM','HRM'),
                   size = res$overlap.sizes[c('110','101','011','111')],
                   expect = res$overlap.expected[c('110','101','011','111')],
                   p = -log10(res$P.value[c('110','101','011','111')]))
temp$FE <- temp$size/temp$expect
temp$set <- factor(temp$set,levels = temp$set)

pdf(file = './res/step_71_fig_221026/RG.1_marker_intersect_supertest_dotplot.pdf',width = 3.6,height = 4)
ggplot(data = temp,aes(x = set,y = size,size = FE,color = p)) + 
  geom_point() + 
  scale_color_gradient(low = 'blue',high = 'red') + 
  theme_bw() + 
  guides(color = guide_colorbar(title = '-log10 P-value')) + 
  guides(size = guide_legend(title = 'Fold enrichment')) + 
  theme(aspect.ratio = 2,
        axis.text = element_text(size = 10),
        axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
  labs(title = 'RG-1 marker')
dev.off()

#IP
marker_list <- IP_marker
res <- supertest(x = marker_list,n = length(gene_list))

pdf(file = './res/step_71_fig_221026/IP_marker_intersect_supertest.pdf',width = 3.5,height = 4)
plot(res, Layout="landscape", degree=2:3,sort.by="size",margin=c(0.5,5,1,2))
dev.off()

temp <- data.frame(set = c('HR','HM','RM','HRM'),
                   size = res$overlap.sizes[c('110','101','011','111')],
                   expect = res$overlap.expected[c('110','101','011','111')],
                   p = -log10(res$P.value[c('110','101','011','111')]))
temp$FE <- temp$size/temp$expect
temp$set <- factor(temp$set,levels = temp$set)

pdf(file = './res/step_71_fig_221026/IP_marker_intersect_supertest_dotplot.pdf',width = 3.6,height = 4)
ggplot(data = temp,aes(x = set,y = size,size = FE,color = p)) + 
  geom_point() + 
  scale_color_gradient(low = 'blue',high = 'red') + 
  theme_bw() + 
  guides(color = guide_colorbar(title = '-log10 P-value')) + 
  guides(size = guide_legend(title = 'Fold enrichment')) + 
  theme(aspect.ratio = 2,
        axis.text = element_text(size = 10),
        axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
  labs(title = 'IP marker')
dev.off()

#Ex-4
marker_list <- Ex.4_marker
res <- supertest(x = marker_list,n = length(gene_list))

pdf(file = './res/step_71_fig_221026/Ex.4_marker_intersect_supertest.pdf',width = 3.5,height = 4)
plot(res, Layout="landscape", degree=2:3,sort.by="size",margin=c(0.5,5,1,2))
dev.off()

temp <- data.frame(set = c('HR','HM','RM','HRM'),
                   size = res$overlap.sizes[c('110','101','011','111')],
                   expect = res$overlap.expected[c('110','101','011','111')],
                   p = -log10(res$P.value[c('110','101','011','111')]))
temp$FE <- temp$size/temp$expect
temp$set <- factor(temp$set,levels = temp$set)

pdf(file = './res/step_71_fig_221026/Ex.4_marker_intersect_supertest_dotplot.pdf',width = 3.6,height = 4)
ggplot(data = temp,aes(x = set,y = size,size = FE,color = p)) + 
  geom_point() + 
  scale_color_gradient(low = 'blue',high = 'red') + 
  theme_bw() + 
  guides(color = guide_colorbar(title = '-log10 P-value')) + 
  guides(size = guide_legend(title = 'Fold enrichment')) + 
  theme(aspect.ratio = 2,
        axis.text = element_text(size = 10),
        axis.title = element_text(face = 'bold',size = 14),
        plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
  labs(title = 'Ex-4 marker')
dev.off()

# global marker featureplot -----------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#migrating neuron
p1 <- FeaturePlot(features = 'NRP1',order = TRUE,object = Greenleaf_RNA_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p2 <- FeaturePlot(features = 'NRP1',order = TRUE,object = macaque_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p3 <- FeaturePlot(features = 'Nrp1',order = TRUE,object = mouse_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

png(filename = './res/step_71_fig_221026/NRP1_three_species_featureplot.png',width = 1500,height = 500)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#upperlayer neuron
p1 <- FeaturePlot(features = 'SATB2',order = TRUE,object = Greenleaf_RNA_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p2 <- FeaturePlot(features = 'SATB2',order = TRUE,object = macaque_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p3 <- FeaturePlot(features = 'Satb2',order = TRUE,object = mouse_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

png(filename = './res/step_71_fig_221026/SATB24_three_species_featureplot.png',width = 1500,height = 500)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#deeplayer neuron
p1 <- FeaturePlot(features = 'FOXP2',order = TRUE,object = Greenleaf_RNA_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p2 <- FeaturePlot(features = 'FOXP2',order = TRUE,object = macaque_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p3 <- FeaturePlot(features = 'Foxp2',order = TRUE,object = mouse_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

png(filename = './res/step_71_fig_221026/FOXP2_three_species_featureplot.png',width = 1500,height = 500)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#newborn neuron
p1 <- FeaturePlot(features = 'NEUROD2',order = TRUE,object = Greenleaf_RNA_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p2 <- FeaturePlot(features = 'NEUROD2',order = TRUE,object = macaque_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p3 <- FeaturePlot(features = 'Neurod2',order = TRUE,object = mouse_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

png(filename = './res/step_71_fig_221026/Neurod2_three_species_featureplot.png',width = 1500,height = 500)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#Inter neuron
p1 <- FeaturePlot(features = 'GAD2',order = TRUE,object = Greenleaf_RNA_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p2 <- FeaturePlot(features = 'GAD2',order = TRUE,object = macaque_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p3 <- FeaturePlot(features = 'Gad2',order = TRUE,object = mouse_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

png(filename = './res/step_71_fig_221026/GAD2_three_species_featureplot.png',width = 1500,height = 500)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

# OPC in mouse ------------------------------------------------------------
#load data
mouse_nature_E16_before <- readRDS(file = './processed_data/221008_summary/mouse_Nature_E14_E15_E16_RNA_Seurat_221017.rds')

#umap
x_lower <- min(mouse_nature_E16_before@reductions$umap@cell.embeddings[,1])
x_upper <- max(mouse_nature_E16_before@reductions$umap@cell.embeddings[,1])
y_lower <- min(mouse_nature_E16_before@reductions$umap@cell.embeddings[,2])
y_upper <- max(mouse_nature_E16_before@reductions$umap@cell.embeddings[,2])

p1 <- DimPlot(object = mouse_nature_E16_before,group.by = 'New_cellType',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15))) + 
  labs(title = 'cell type')

p2 <- FeaturePlot(features = 'Olig2',order = FALSE,object = mouse_nature_E16_before,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none') + 
  NoAxes() + 
  theme(plot.background = element_rect(size = 1,fill = NA,color = 'black'))

png(filename = './res/step_71_fig_221026/mouse_nature_before_E16_OPC_marker.png',width = 1000,height = 500)
p1+p2+plot_layout(ncol = 2)
dev.off()

#data after E17
mouse_nature_data <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/mouse_RNA_seurat.rds')
mouse_nature_data <- mouse_nature_data[,mouse_nature_data$donor_id %in% c('mouse_E17','mouse_E18_S1','mouse_E18_S3')]

#plot
x_lower <- min(mouse_nature_data@reductions$umap@cell.embeddings[,1])
x_upper <- max(mouse_nature_data@reductions$umap@cell.embeddings[,1])
y_lower <- min(mouse_nature_data@reductions$umap@cell.embeddings[,2])
y_upper <- max(mouse_nature_data@reductions$umap@cell.embeddings[,2])

p1 <- DimPlot(object = mouse_nature_data,group.by = 'New_cellType',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15))) + 
  labs(title = 'cell type')

p2 <- FeaturePlot(features = 'Olig2',order = FALSE,object = mouse_nature_data,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none') + 
  NoAxes() + 
  theme(plot.background = element_rect(size = 1,fill = NA,color = 'black'))

p3 <- FeaturePlot(features = 'Egfr',order = FALSE,object = mouse_nature_data,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none') + 
  NoAxes() + 
  theme(plot.background = element_rect(size = 1,fill = NA,color = 'black'))

png(filename = './res/step_71_fig_221026/mouse_nature_after_E17_OPC_marker.png',width = 1500,height = 500)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

# progenitor cell UMAP ----------------------------------------------------
macaque_progenitor <- macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% c('RG-1','RG-2','IP','OPC','Cycling')]
macaque_progenitor <- my_process_seurat(object = macaque_progenitor,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_progenitor <- my_process_seurat(object = macaque_progenitor,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 10,resolution = 1,group.by = 'cell_type',label = TRUE)

mouse_progenitor <- mouse_nature_data[,mouse_nature_data$New_cellType %in% c('Cycling glial cells','Intermediate progenitors','Oligodendrocytes')]

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
#convert to macaque genes
meta_data <- mouse_progenitor@meta.data
mouse_progenitor <- mouse_progenitor@assays$RNA@counts
mouse_to_macaque_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_Mmul10.csv')
mouse_to_macaque_anno <- mouse_to_macaque_anno[,c(2,1,4,3)]

mouse_progenitor <- My_Convert_Homology_Gene_ID(express_matrix = mouse_progenitor,anno = mouse_to_macaque_anno,filter_anno = TRUE,future.globals.maxSize = 200*(1024^3),workers = 5)

#recreate mouse seurat
mouse_progenitor <- CreateSeuratObject(counts = mouse_progenitor,project = 'mouse',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#reprocess macaque
macaque_progenitor <- FindVariableFeatures(object = macaque_progenitor,nfeatures = 2000)
gene_list <- VariableFeatures(macaque_progenitor)[VariableFeatures(macaque_progenitor) %in% rownames(mouse_progenitor)]

macaque_progenitor <- my_process_seurat(object = macaque_progenitor,assay = 'RNA',reduction.name = 'pca',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_progenitor <- my_process_seurat(object = macaque_progenitor,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 10,resolution = 1,group.by = 'cell_type',label = TRUE)

mouse_progenitor <- my_process_seurat(object = mouse_progenitor,assay = 'RNA',reduction.name = 'pca',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#projection
temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_progenitor@assays$RNA@scale.data,object = macaque_progenitor,assayUsed = 'RNA',missing_gene = FALSE)
mouse_progenitor@reductions$projected_pca <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
mouse_progenitor@reductions$projected_umap <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')
DimPlot(mouse_progenitor,group.by = 'New_cellType',label = TRUE,repel = TRUE,reduction = 'projected_umap')

FeaturePlot(object = mouse_progenitor,features = 'OLIG2',reduction = 'projected_umap')

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_progenitor,query = mouse_progenitor,
                                  ref_reduction = 'pca',query_reduction = 'projected_pca',
                                  ref_assay = 'RNA',query_assay = 'RNA',l2.norm = FALSE,
                                  dims = 1:10,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_progenitor$cell_type,l2.norm = FALSE,dims = 1:10,verbose = TRUE)
mouse_progenitor <- AddMetaData(object = mouse_progenitor,metadata = predictions)

DimPlot(mouse_progenitor,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'projected_umap')

#plot integration
x_lower <- min(macaque_progenitor@reductions$umap@cell.embeddings[,1])
x_upper <- max(macaque_progenitor@reductions$umap@cell.embeddings[,1])
y_lower <- min(macaque_progenitor@reductions$umap@cell.embeddings[,2])
y_upper <- max(macaque_progenitor@reductions$umap@cell.embeddings[,2])

p1 <- DimPlot(object = macaque_progenitor,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = color_param$celltype) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15))) + 
  NoLegend()

x_lower <- min(mouse_progenitor@reductions$projected_umap@cell.embeddings[,1])
x_upper <- max(mouse_progenitor@reductions$projected_umap@cell.embeddings[,1])
y_lower <- min(mouse_progenitor@reductions$projected_umap@cell.embeddings[,2])
y_upper <- max(mouse_progenitor@reductions$projected_umap@cell.embeddings[,2])

p2 <- DimPlot(object = mouse_progenitor,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'projected_umap') + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = color_param$celltype[c('RG-1','RG-2','Cycling','IP','OPC')]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

pdf(file = './res/step_71_fig_221026/progenitor_cell_umap.pdf',width = 5,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

mouse_progenitor$predicted.id <- factor(mouse_progenitor$predicted.id,levels = c('RG-1','RG-2','Cycling','IP','OPC'))
pdf(file = './res/step_71_fig_221026/progenitor_cell_vlnplot.pdf',width = 6,height = 5)
VlnPlot(object = mouse_progenitor,features = c('PAX6','EOMES','OLIG2','EGFR','TOP2A'),group.by = 'predicted.id',cols = color_param$celltype[c('RG-1','RG-2','Cycling','IP','OPC')])
dev.off()

#proportion
temp_human <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
temp_macaque <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
temp_mouse <- mouse_nature_data
temp_mouse@meta.data[colnames(mouse_progenitor),"New_cellType"] <- as.character(mouse_progenitor$predicted.id)

temp <- data.frame(species = c('human','human','human','macaque','macaque','macaque','macaque','mouse','mouse'),
                   age = c('pcw16','pcw20','pcw21','E80','E80','E90','E90','E17.5','E18.5'),
                   percent = c(sum(temp_human$Age == 'pcw16' & temp_human$ReAnno_celltype %in% c('OPC','RG-2'))/sum(temp_human$Age == 'pcw16'),
                               sum(temp_human$Age == 'pcw20' & temp_human$ReAnno_celltype %in% c('OPC','RG-2'))/sum(temp_human$Age == 'pcw20'),
                               sum(temp_human$Age == 'pcw21' & temp_human$ReAnno_celltype %in% c('OPC','RG-2'))/sum(temp_human$Age == 'pcw21'),
                               sum(temp_macaque$donor == 'A50A' & temp_macaque$cell_type %in% c('OPC','RG-2'))/sum(temp_macaque$donor == 'A50A'),
                               sum(temp_macaque$donor == 'A50B' & temp_macaque$cell_type %in% c('OPC','RG-2'))/sum(temp_macaque$donor == 'A50B'),
                               sum(temp_macaque$donor == 'A82A' & temp_macaque$cell_type %in% c('OPC','RG-2'))/sum(temp_macaque$donor == 'A82A'),
                               sum(temp_macaque$donor == 'A82B' & temp_macaque$cell_type %in% c('OPC','RG-2'))/sum(temp_macaque$donor == 'A82B'),
                               sum(temp_mouse$orig.ident == 'E17' & temp_mouse$New_cellType %in% c('OPC','RG-2'))/sum(temp_mouse$orig.ident == 'E17'),
                               sum(temp_mouse$orig.ident == 'E18' & temp_mouse$New_cellType %in% c('OPC','RG-2'))/sum(temp_mouse$orig.ident == 'E18')))
temp$percent <- temp$percent*100
temp$species <- factor(temp$species,levels = unique(temp$species))
temp$age <- factor(temp$age,levels = unique(temp$age))

pdf(file = './res/step_71_fig_221026/OPC_pre_OPC_proportion_dotplot.pdf',width = 4.5,height = 4)
ggplot(data = temp,aes(x = species,y = percent,color = age,shape = species)) + 
  geom_point(size = 4) + 
  theme_bw() + 
  scale_color_manual(values = c('pcw16' = '#a6cde2','pcw20' = '#62a2cb','pcw21' = '#1e78b4',
                                'E80' = '#74c476','E90' = '#34a047',
                                'E17.5' = '#cab2d6','E18.5' = '#6a3e98')) + 
  theme(aspect.ratio = 1)
dev.off()


# UMAP ----------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = '/home/sunym/temp/merged_Greenleaf_ATAC_ArchR_221018/')
mouse_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/mouse_multiome_ArchR_221009/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#umap
p1 <- my_dimplot(embedding = mouse_multiome_ArchR@embeddings$UMAP$df,meta_data = mouse_multiome_ArchR@cellColData,
                 group.by = 'Gex_macaque_cell_type',label = FALSE) + 
  scale_color_manual(values = color_param$celltype) + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1) + 
  labs(title = 'mouse')

p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = macaque_multiome_ArchR@cellColData,
                 group.by = 'cell_type',label = FALSE) + 
  scale_color_manual(values = color_param$celltype) + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque')

p3 <- my_dimplot(embedding = Greenleaf_ATAC_ArchR@embeddings$UMAP$df,meta_data = Greenleaf_ATAC_ArchR@cellColData,
                 group.by = 'cell_type',label = FALSE) + 
  scale_color_manual(values = color_param$celltype) + 
  NoAxes() + 
  theme(aspect.ratio = 1) + 
  labs(title = 'human')

png(filename = './res/step_71_fig_221026/ATAC_data_umapplot.png',width = 1000,height = 500)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#add peakset
peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = peakset$human,force = TRUE)
macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = peakset$macaque,force = TRUE)
mouse_multiome_ArchR <- addPeakSet(ArchRProj = mouse_multiome_ArchR,peakSet = peakset$mouse,force = TRUE)

Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR,force = TRUE)
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,force = TRUE)
mouse_multiome_ArchR <- addPeakMatrix(ArchRProj = mouse_multiome_ArchR,force = TRUE)

#get marker
human_markersPeaks <- getMarkerFeatures(
  ArchRProj = Greenleaf_ATAC_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
macaque_markersPeaks <- getMarkerFeatures(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
mouse_markersPeaks <- getMarkerFeatures(
  ArchRProj = mouse_multiome_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "Gex_macaque_cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

human_markerList <- getMarkers(human_markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
macaque_markerList <- getMarkers(macaque_markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
mouse_markerList <- getMarkers(mouse_markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = macaque_markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")


# PPP1R17 -----------------------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#PPP1R17
p1 <- FeaturePlot(features = 'PPP1R17',order = TRUE,object = Greenleaf_RNA_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p2 <- FeaturePlot(features = 'PPP1R17',order = TRUE,object = macaque_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

p3 <- FeaturePlot(features = 'Ppp1r17',order = TRUE,object = mouse_multiome_Seurat,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

png(filename = './res/step_71_fig_221026/PPP1R17_three_species_featureplot.png',width = 500,height = 1500)
p1+p2+p3+plot_layout(ncol = 1)
dev.off()

# mouse nature data tRG marker --------------------------------------------
#load data
mouse_nature_data <- readRDS(file = './processed_data/221008_summary/mouse_Nature_E14_E15_E16_RNA_Seurat_221017.rds')

#feature plot
p1 <- FeaturePlot(features = 'Cryab',order = FALSE,object = mouse_nature_data,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')
p2 <- FeaturePlot(features = 'Fndc1',order = FALSE,object = mouse_nature_data,
                  cols = as.character(ArchRPalettes$whitePurple[-1]),
                  pt.size = 0.001,slot = 'data') + 
  NoAxes() + 
  theme(plot.title = element_blank(),
        aspect.ratio = 1,legend.position = 'none')

png(filename = './res/step_71_fig_221026/mouse_nature_tRG_marker.png',width = 1000,height = 500)
p1+p2+plot_layout(ncol = 2)
dev.off()

#load data
Brain_RNA_Seurat <- readRDS(file = '/content/load_data/Greenleaf_macaque-multiome_mouse-multiome_Seurat_human_symbol_progenitor_221017.rds')

png(filename = './res/step_71_fig_221026/Brain_progenitor_CRYAB_featureplot.png',width = 1500,height = 500)
FeaturePlot(features = 'CRYAB',order = TRUE,object = Brain_RNA_Seurat,
            cols = as.character(ArchRPalettes$whitePurple[-1]),
            pt.size = 0.001,slot = 'data',split.by = 'species')
dev.off()

png(filename = './res/step_71_fig_221026/Brain_progenitor_FNDC1_featureplot.png',width = 1500,height = 500)
FeaturePlot(features = 'FNDC1',order = FALSE,object = Brain_RNA_Seurat,
            cols = as.character(ArchRPalettes$whitePurple[-1]),
            pt.size = 0.001,slot = 'data',split.by = 'species')
dev.off()

#test jupyter