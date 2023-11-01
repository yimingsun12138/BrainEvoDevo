#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: re_integrate three species RNA_seq to prove param insensitive   ##
## Data: 2023.02.07                                                                ##
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

#library
library(Seurat)
library(ggplot2)
library(dplyr)
library(parallel)
library(rliger)
library(harmony)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

# load data ---------------------------------------------------------------
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_human_symbol_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_human_symbol_221009.rds')

#add meta data
Greenleaf_RNA_Seurat$species <- 'human'
macaque_multiome_Seurat$species <- 'macaque'
mouse_multiome_Seurat$species <- 'mouse'

#create seurat list
seurat.list <- list(Greenleaf_RNA_Seurat,macaque_multiome_Seurat,mouse_multiome_Seurat)

# CCA k=5 ndim=40 ---------------------------------------------------------

#set param
k_anchors <- 5
ndim <- 40

#preprocess
seurat.list <- base::lapply(X = seurat.list,FUN = function(x){
  x <- NormalizeData(object = x,normalization.method = "LogNormalize",scale.factor = 10000)
  x <- FindVariableFeatures(object = x,selection.method = "vst",nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = seurat.list)

#do integration
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list,assay = 'RNA',anchor.features = features,
                                         normalization.method = 'LogNormalize',reduction = 'cca',
                                         dims = 1:ndim,k.anchor = k_anchors,verbose = TRUE)
seurat.combined <- IntegrateData(anchorset = seurat.anchors,new.assay.name = 'integrated',
                                 normalization.method = 'LogNormalize',dims = 1:ndim,verbose = TRUE)

#preprocess integrated data
DefaultAssay(seurat.combined) <- 'integrated'
seurat.combined <- ScaleData(object = seurat.combined)
seurat.combined <- RunPCA(object = seurat.combined,npcs = 50,reduction.name = "pca",verbose = TRUE)
seurat.combined <- FindNeighbors(object = seurat.combined,dims = 1:ndim,reduction = "pca")
seurat.combined <- FindClusters(object = seurat.combined,resolution = 1)
seurat.combined <- RunUMAP(object = seurat.combined,dims = 1:ndim,reduction = "pca",reduction.name = "umap")

#save data
char <- paste0('/home/sunym/temp/','CCA_k_',k_anchors,'_ndim_',ndim,'_seurat_obj.rds')
saveRDS(object = seurat.combined,file = char)
print('CCA integration done!')

# CCA k=5 ndim=50 ---------------------------------------------------------

#set param
k_anchors <- 5
ndim <- 50

#preprocess
seurat.list <- base::lapply(X = seurat.list,FUN = function(x){
  x <- NormalizeData(object = x,normalization.method = "LogNormalize",scale.factor = 10000)
  x <- FindVariableFeatures(object = x,selection.method = "vst",nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = seurat.list)

#do integration
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list,assay = 'RNA',anchor.features = features,
                                         normalization.method = 'LogNormalize',reduction = 'cca',
                                         dims = 1:ndim,k.anchor = k_anchors,verbose = TRUE)
seurat.combined <- IntegrateData(anchorset = seurat.anchors,new.assay.name = 'integrated',
                                 normalization.method = 'LogNormalize',dims = 1:ndim,verbose = TRUE)

#preprocess integrated data
DefaultAssay(seurat.combined) <- 'integrated'
seurat.combined <- ScaleData(object = seurat.combined)
seurat.combined <- RunPCA(object = seurat.combined,npcs = 50,reduction.name = "pca",verbose = TRUE)
seurat.combined <- FindNeighbors(object = seurat.combined,dims = 1:ndim,reduction = "pca")
seurat.combined <- FindClusters(object = seurat.combined,resolution = 1)
seurat.combined <- RunUMAP(object = seurat.combined,dims = 1:ndim,reduction = "pca",reduction.name = "umap")

#save data
char <- paste0('/home/sunym/temp/','CCA_k_',k_anchors,'_ndim_',ndim,'_seurat_obj.rds')
saveRDS(object = seurat.combined,file = char)
print('CCA integration done!')

# liger integration k 30 lambda 5 -------------------------------------------------------

#set param
k_dim <- 30
lambda_liger <- 5

#create liger object
Brain_liger <- createLiger(list(human = Greenleaf_RNA_Seurat@assays$RNA@counts,
                                macaque = macaque_multiome_Seurat@assays$RNA@counts,
                                mouse = mouse_multiome_Seurat@assays$RNA@counts))

#preprocess liger object
Brain_liger <- normalize(object = Brain_liger)
Brain_liger <- selectGenes(object = Brain_liger)
Brain_liger <- scaleNotCenter(object = Brain_liger)
gc()

#run NMF
Brain_liger <- optimizeALS(object = Brain_liger,k = k_dim,lambda = lambda_liger,verbose = TRUE)

#quantile normalize and joint clustering
Brain_liger <- quantile_norm(object = Brain_liger)
Brain_liger <- louvainCluster(object = Brain_liger,resolution = 0.5)

# #visualize by UMAP
# Brain_liger <- runUMAP(object = Brain_liger,distance = 'cosine',n_neighbors = 30,min_dist = 0.3)
# all.plots <- plotByDatasetAndCluster(Brain_liger,axis.labels = c('UMAP 1','UMAP 2'),return.plots = TRUE)
# all.plots[[1]] + all.plots[[2]]

#convert to seurat object
Brain_Seurat <- ligerToSeurat(object = Brain_liger)
Brain_Seurat <- RunUMAP(object = Brain_Seurat,dims = 1:k_dim,reduction = 'inmf')

#add meta data
Brain_Seurat$cell_type <- NA
Brain_Seurat@meta.data[paste('human',colnames(Greenleaf_RNA_Seurat),sep = '_'),"cell_type"] <- Greenleaf_RNA_Seurat$ReAnno_celltype
Brain_Seurat@meta.data[paste('macaque',colnames(macaque_multiome_Seurat),sep = '_'),"cell_type"] <- macaque_multiome_Seurat$cell_type
Brain_Seurat@meta.data[paste('mouse',colnames(mouse_multiome_Seurat),sep = '_'),"cell_type"] <- mouse_multiome_Seurat$macaque_cell_type

# DimPlot(Brain_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

#save data
char <- paste0('/home/sunym/temp/liger_k_',k_dim,'_lambda_',lambda_liger,'_seurat_obj.rds')
saveRDS(object = Brain_Seurat,file = char)
print('all done!')

# liger integration k 25 lambda 20 -------------------------------------------------------

#set param
k_dim <- 25
lambda_liger <- 20

#create liger object
Brain_liger <- createLiger(list(human = Greenleaf_RNA_Seurat@assays$RNA@counts,
                                macaque = macaque_multiome_Seurat@assays$RNA@counts,
                                mouse = mouse_multiome_Seurat@assays$RNA@counts))

#preprocess liger object
Brain_liger <- normalize(object = Brain_liger)
Brain_liger <- selectGenes(object = Brain_liger)
Brain_liger <- scaleNotCenter(object = Brain_liger)
gc()

#run NMF
Brain_liger <- optimizeALS(object = Brain_liger,k = k_dim,lambda = lambda_liger,verbose = TRUE)

#quantile normalize and joint clustering
Brain_liger <- quantile_norm(object = Brain_liger)
Brain_liger <- louvainCluster(object = Brain_liger,resolution = 0.5)

#visualize by UMAP
Brain_liger <- runUMAP(object = Brain_liger,distance = 'cosine',n_neighbors = 30,min_dist = 0.3)
all.plots <- plotByDatasetAndCluster(Brain_liger,axis.labels = c('UMAP 1','UMAP 2'),return.plots = TRUE)
all.plots[[1]] + all.plots[[2]]

#convert to seurat object
Brain_Seurat <- ligerToSeurat(object = Brain_liger)
Brain_Seurat <- RunUMAP(object = Brain_Seurat,dims = 1:k_dim,reduction = 'inmf')

#add meta data
Brain_Seurat$cell_type <- NA
Brain_Seurat@meta.data[paste('human',colnames(Greenleaf_RNA_Seurat),sep = '_'),"cell_type"] <- Greenleaf_RNA_Seurat$ReAnno_celltype
Brain_Seurat@meta.data[paste('macaque',colnames(macaque_multiome_Seurat),sep = '_'),"cell_type"] <- macaque_multiome_Seurat$cell_type
Brain_Seurat@meta.data[paste('mouse',colnames(mouse_multiome_Seurat),sep = '_'),"cell_type"] <- mouse_multiome_Seurat$macaque_cell_type

DimPlot(Brain_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

#save data
char <- paste0('/home/sunym/temp/liger_k_',k_dim,'_lambda_',lambda_liger,'_seurat_obj.rds')
saveRDS(object = Brain_Seurat,file = char)
print('all done!')

# harmony integration -----------------------------------------------------

#find variable features
macaque_multiome_Seurat <- FindVariableFeatures(object = macaque_multiome_Seurat,selection.method = 'vst',assay = 'RNA',nfeatures = 2100,verbose = TRUE)
gene_list <- VariableFeatures(object = macaque_multiome_Seurat)
gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)]
gene_list <- gene_list[gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts)]

#process macaque multiome Seurat
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 32,resolution = 1,group.by = 'cell_type',label = TRUE)

#harmony integration
#sigma = 0.2,theta = 0.18
Brain_harmony <- my_harmony_integration(named_seurat_list = list(human = Greenleaf_RNA_Seurat,macaque = macaque_multiome_Seurat,mouse = mouse_multiome_Seurat),
                                        assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(human = NULL,macaque = c('nCount_RNA','donor'),mouse = NULL),
                                        npcs = 50,reference_loading = 'macaque',integration_var = 'species',harmony_input_dim = 30,max.iter.harmony = 100,
                                        UMAP_dim = 30,resolution = 1,sigma = 0.2,theta = 0.18)

#add meta data
Brain_harmony@meta.data[1:3,]
Brain_harmony$cell_type <- NA
Brain_harmony@meta.data[colnames(Greenleaf_RNA_Seurat),"cell_type"] <- Greenleaf_RNA_Seurat$ReAnno_celltype
Brain_harmony@meta.data[colnames(macaque_multiome_Seurat),"cell_type"] <- macaque_multiome_Seurat$cell_type
Brain_harmony@meta.data[colnames(mouse_multiome_Seurat),"cell_type"] <- mouse_multiome_Seurat$macaque_cell_type

DimPlot(object = Brain_harmony,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)
DimPlot(object = Brain_harmony,group.by = 'dataset',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

# #save data
# saveRDS(object = Brain_harmony,file = '/content/data/sunym/temp/harmony_integration_230209.rds')

# final param -------------------------------------------------------------

## CCA k 5 ndim 40 ---------------------------------------------------------

#load data
Brain_Seurat <- readRDS(file = '/content/data/sunym/temp/CCA_k_5_ndim_40_seurat_obj.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#merged species plot
p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_40/merged_species_dimplot.png',plot = p,width = 3.5,height = 3.5)

p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'CCA k=5 ndim=40')
pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_40/merged_species_dimplot.pdf',width = 4,height = 4)
p
dev.off()

#predicted cell type dimplot, split by species
temp_macaque <- Brain_Seurat[,Brain_Seurat$species == 'macaque']
temp_human <- Brain_Seurat[,Brain_Seurat$species == 'human']
temp_mouse <- Brain_Seurat[,Brain_Seurat$species == 'mouse']
temp_macaque$cell_type <- macaque_multiome_Seurat@meta.data[rownames(temp_macaque@meta.data),"cell_type"]

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integrated',query_assay = 'integrated',l2.norm = TRUE,dims = 1:40,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = TRUE,dims = 1:40,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_mouse,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integrated',query_assay = 'integrated',l2.norm = TRUE,dims = 1:40,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = TRUE,dims = 1:40,verbose = TRUE)
temp_mouse <- AddMetaData(object = temp_mouse,metadata = predictions)

temp_macaque$cell_type <- factor(temp_macaque$cell_type,levels = names(col_param$celltype))
temp_human$predicted.id <- factor(temp_human$predicted.id,levels = names(col_param$celltype))
temp_mouse$predicted.id <- factor(temp_mouse$predicted.id,levels = names(col_param$celltype))

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_40/human_predicted_cell_type_dimplot.png',plot = p1,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_40/macaque_predicted_cell_type_dimplot.png',plot = p2,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_40/mouse_predicted_cell_type_dimplot.png',plot = p3,width = 3.5,height = 3.5)

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'human predicted cell type')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque cell type')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$celltype,drop = FALSE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'mouse predicted cell type')

pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_40/predicted_cell_type_dimplot_splited_by_species.pdf',width = 12,height = 4)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#predicted cell type confusion heatmap
temp_human$ori_cell_type <- Greenleaf_RNA_Seurat@meta.data[rownames(temp_human@meta.data),"ReAnno_celltype"]
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
temp_mouse$ori_cell_type <- meta_data[rownames(temp_mouse@meta.data),"ReAnno_celltype"]

ori <- temp_human$ori_cell_type
prd <- temp_human$predicted.id
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
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'Greenleaf scRNAseq') + 
  scale_x_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','Mic'))

pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_40/human_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()

ori <- temp_mouse$ori_cell_type
prd <- temp_mouse$predicted.id
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic'))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(col_param$celltype))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_gradient(low = '#F1F1F1',high = 'black') + 
  theme(aspect.ratio = 15/13) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'mouse multiome') + 
  scale_x_discrete(labels = c('RG-neu','Cyc-RG','Cyc-IP','IP','ExN','ExM','ExUp','ExDp','In','End','Per','VLMC','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','VLMC','Mic'))

pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_40/mouse_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()

## CCA k 5 ndim 50 ---------------------------------------------------------

#load data
Brain_Seurat <- readRDS(file = '/content/data/sunym/temp/CCA_k_5_ndim_50_seurat_obj.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#merged species plot
p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_50/merged_species_dimplot.png',plot = p,width = 3.5,height = 3.5)

p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'CCA k=5 ndim=50')
pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_50/merged_species_dimplot.pdf',width = 4,height = 4)
p
dev.off()

#predicted cell type dimplot, split by species
temp_macaque <- Brain_Seurat[,Brain_Seurat$species == 'macaque']
temp_human <- Brain_Seurat[,Brain_Seurat$species == 'human']
temp_mouse <- Brain_Seurat[,Brain_Seurat$species == 'mouse']
temp_macaque$cell_type <- macaque_multiome_Seurat@meta.data[rownames(temp_macaque@meta.data),"cell_type"]

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integrated',query_assay = 'integrated',l2.norm = TRUE,dims = 1:50,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_mouse,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integrated',query_assay = 'integrated',l2.norm = TRUE,dims = 1:50,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = TRUE,dims = 1:50,verbose = TRUE)
temp_mouse <- AddMetaData(object = temp_mouse,metadata = predictions)

temp_macaque$cell_type <- factor(temp_macaque$cell_type,levels = names(col_param$celltype))
temp_human$predicted.id <- factor(temp_human$predicted.id,levels = names(col_param$celltype))
temp_mouse$predicted.id <- factor(temp_mouse$predicted.id,levels = names(col_param$celltype))

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_50/human_predicted_cell_type_dimplot.png',plot = p1,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_50/macaque_predicted_cell_type_dimplot.png',plot = p2,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/CCA_k_5_ndim_50/mouse_predicted_cell_type_dimplot.png',plot = p3,width = 3.5,height = 3.5)

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'human predicted cell type')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque cell type')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$celltype,drop = FALSE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'mouse predicted cell type')

pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_50/predicted_cell_type_dimplot_splited_by_species.pdf',width = 12,height = 4)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#predicted cell type confusion heatmap
temp_human$ori_cell_type <- Greenleaf_RNA_Seurat@meta.data[rownames(temp_human@meta.data),"ReAnno_celltype"]
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
temp_mouse$ori_cell_type <- meta_data[rownames(temp_mouse@meta.data),"ReAnno_celltype"]

ori <- temp_human$ori_cell_type
prd <- temp_human$predicted.id
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
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'Greenleaf scRNAseq') + 
  scale_x_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','Mic'))

pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_50/human_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()

ori <- temp_mouse$ori_cell_type
prd <- temp_mouse$predicted.id
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic'))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(col_param$celltype))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_gradient(low = '#F1F1F1',high = 'black') + 
  theme(aspect.ratio = 15/13) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'mouse multiome') + 
  scale_x_discrete(labels = c('RG-neu','Cyc-RG','Cyc-IP','IP','ExN','ExM','ExUp','ExDp','In','End','Per','VLMC','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','VLMC','Mic'))

pdf(file = './res/step_101_fig_230207/CCA_k_5_ndim_50/mouse_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()

## liger k 35 lambda 20 ----------------------------------------------------

#load data
Brain_Seurat <- readRDS(file = '/content/data/sunym/temp/liger_k_35_lambda_20_seurat_obj.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

Brain_Seurat$species <- Brain_Seurat$orig.ident
Brain_Seurat <- RenameCells(object = Brain_Seurat,new.names = sub(pattern = '^human_',replacement = '',x = colnames(Brain_Seurat@assays$RNA@counts),fixed = FALSE))
Brain_Seurat <- RenameCells(object = Brain_Seurat,new.names = sub(pattern = '^macaque_',replacement = '',x = colnames(Brain_Seurat@assays$RNA@counts),fixed = FALSE))
Brain_Seurat <- RenameCells(object = Brain_Seurat,new.names = sub(pattern = '^mouse_',replacement = '',x = colnames(Brain_Seurat@assays$RNA@counts),fixed = FALSE))

#merged species plot
p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/liger/merged_species_dimplot.png',plot = p,width = 3.5,height = 3.5)

p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'liger k 35 lambda 20')
pdf(file = './res/step_101_fig_230207/liger/merged_species_dimplot.pdf',width = 4,height = 4)
p
dev.off()

#predicted cell type dimplot, split by species
temp_macaque <- Brain_Seurat[,Brain_Seurat$species == 'macaque']
temp_human <- Brain_Seurat[,Brain_Seurat$species == 'human']
temp_mouse <- Brain_Seurat[,Brain_Seurat$species == 'mouse']
temp_macaque$cell_type <- macaque_multiome_Seurat@meta.data[rownames(temp_macaque@meta.data),"cell_type"]

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'inmf',query_reduction = 'inmf',
                                  ref_assay = 'RNA',query_assay = 'RNA',l2.norm = TRUE,dims = 1:35,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = TRUE,dims = 1:35,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_mouse,ref_reduction = 'inmf',query_reduction = 'inmf',
                                  ref_assay = 'RNA',query_assay = 'RNA',l2.norm = TRUE,dims = 1:35,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = TRUE,dims = 1:35,verbose = TRUE)
temp_mouse <- AddMetaData(object = temp_mouse,metadata = predictions)

temp_macaque$cell_type <- factor(temp_macaque$cell_type,levels = names(col_param$celltype))
temp_human$predicted.id <- factor(temp_human$predicted.id,levels = names(col_param$celltype))
temp_mouse$predicted.id <- factor(temp_mouse$predicted.id,levels = names(col_param$celltype))

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/liger/human_predicted_cell_type_dimplot.png',plot = p1,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/liger/macaque_predicted_cell_type_dimplot.png',plot = p2,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/liger/mouse_predicted_cell_type_dimplot.png',plot = p3,width = 3.5,height = 3.5)

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'human predicted cell type')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque cell type')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$celltype,drop = FALSE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'mouse predicted cell type')

pdf(file = './res/step_101_fig_230207/liger/predicted_cell_type_dimplot_splited_by_species.pdf',width = 12,height = 4)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#predicted cell type confusion heatmap
temp_human$ori_cell_type <- Greenleaf_RNA_Seurat@meta.data[rownames(temp_human@meta.data),"ReAnno_celltype"]
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
temp_mouse$ori_cell_type <- meta_data[rownames(temp_mouse@meta.data),"ReAnno_celltype"]

ori <- temp_human$ori_cell_type
prd <- temp_human$predicted.id
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
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'Greenleaf scRNAseq') + 
  scale_x_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','Mic'))

pdf(file = './res/step_101_fig_230207/liger/human_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()

ori <- temp_mouse$ori_cell_type
prd <- temp_mouse$predicted.id
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic'))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(col_param$celltype))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_gradient(low = '#F1F1F1',high = 'black') + 
  theme(aspect.ratio = 15/13) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'mouse multiome') + 
  scale_x_discrete(labels = c('RG-neu','Cyc-RG','Cyc-IP','IP','ExN','ExM','ExUp','ExDp','In','End','Per','VLMC','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','VLMC','Mic'))

pdf(file = './res/step_101_fig_230207/liger/mouse_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()

## harmony sigma = 0.2 theta = 0.18 ----------------------------------------

#load data
Brain_Seurat <- readRDS(file = '/content/data/sunym/temp/harmony_integration_230209.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

Brain_Seurat$species <- Brain_Seurat$dataset

#merged species plot
p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/harmony/merged_species_dimplot.png',plot = p,width = 3.5,height = 3.5)

p <- DimPlot(object = Brain_Seurat,group.by = 'species',label = FALSE,pt.size = 0.001,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$species) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'harmony sigma = 0.2 theta = 0.18')
pdf(file = './res/step_101_fig_230207/harmony/merged_species_dimplot.pdf',width = 4,height = 4)
p
dev.off()

#predicted cell type dimplot, split by species
temp_macaque <- Brain_Seurat[,Brain_Seurat$species == 'macaque']
temp_human <- Brain_Seurat[,Brain_Seurat$species == 'human']
temp_mouse <- Brain_Seurat[,Brain_Seurat$species == 'mouse']
temp_macaque$cell_type <- macaque_multiome_Seurat@meta.data[rownames(temp_macaque@meta.data),"cell_type"]

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:30,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:30,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)

anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_mouse,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:30,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,l2.norm = FALSE,dims = 1:30,verbose = TRUE)
temp_mouse <- AddMetaData(object = temp_mouse,metadata = predictions)

temp_human@meta.data[temp_human$integration_snn_res.1 %in% c('12','13') & temp_human$predicted.id == 'RG-2',"predicted.id"] <- 'RG-1'
temp_mouse@meta.data[temp_mouse$integration_snn_res.1 %in% c('12','13') & temp_mouse$predicted.id == 'RG-2',"predicted.id"] <- 'RG-1'

temp_macaque$cell_type <- factor(temp_macaque$cell_type,levels = names(col_param$celltype))
temp_human$predicted.id <- factor(temp_human$predicted.id,levels = names(col_param$celltype))
temp_mouse$predicted.id <- factor(temp_mouse$predicted.id,levels = names(col_param$celltype))

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = '')
ggsave(filename = './res/step_101_fig_230207/harmony/human_predicted_cell_type_dimplot.png',plot = p1,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/harmony/macaque_predicted_cell_type_dimplot.png',plot = p2,width = 3.5,height = 3.5)
ggsave(filename = './res/step_101_fig_230207/harmony/mouse_predicted_cell_type_dimplot.png',plot = p3,width = 3.5,height = 3.5)

p1 <- DimPlot(object = temp_human,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'human predicted cell type')
p2 <- DimPlot(object = temp_macaque,pt.size = 0.001,reduction = 'umap',group.by = 'cell_type',label = FALSE,raster = FALSE) + 
  NoAxes() + NoLegend() + 
  scale_color_manual(values = col_param$celltype) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'macaque cell type')
p3 <- DimPlot(object = temp_mouse,pt.size = 0.001,reduction = 'umap',group.by = 'predicted.id',label = FALSE,raster = FALSE) + 
  NoAxes() + 
  scale_color_manual(values = col_param$celltype,drop = FALSE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'mouse predicted cell type')

pdf(file = './res/step_101_fig_230207/harmony/predicted_cell_type_dimplot_splited_by_species.pdf',width = 12,height = 4)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#predicted cell type confusion heatmap
temp_human$ori_cell_type <- Greenleaf_RNA_Seurat@meta.data[rownames(temp_human@meta.data),"ReAnno_celltype"]
meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
temp_mouse$ori_cell_type <- meta_data[rownames(temp_mouse@meta.data),"ReAnno_celltype"]

ori <- temp_human$ori_cell_type
prd <- temp_human$predicted.id
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
  theme(aspect.ratio = 15/14) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'Greenleaf scRNAseq') + 
  scale_x_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','VLMC','Mic'))

pdf(file = './res/step_101_fig_230207/harmony/human_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()

ori <- temp_mouse$ori_cell_type
prd <- temp_mouse$predicted.id
cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
cross.validation.filt[is.na(cross.validation.filt)] = 0
cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
cross.validation.filt$ori <- factor(cross.validation.filt$ori,levels = c('RG-1','Cyc-RG','Cyc-IP','IP','Ex-1','Ex-2','Ex-3','Ex-4','In','End','Per','VLMC','Mic'))
cross.validation.filt$prd <- factor(cross.validation.filt$prd,levels = names(col_param$celltype))

p <- ggplot(data = cross.validation.filt,aes(x = ori, y = prd, fill = Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + theme(panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                       axis.ticks = element_blank(), axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_gradient(low = '#F1F1F1',high = 'black') + 
  theme(aspect.ratio = 15/13) + 
  theme(axis.title = element_text(size = 12,face = 'bold'),
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  xlab('Original Annotation') + ylab('Predicted Annotation') + labs(title = 'mouse multiome') + 
  scale_x_discrete(labels = c('RG-neu','Cyc-RG','Cyc-IP','IP','ExN','ExM','ExUp','ExDp','In','End','Per','VLMC','Mic')) + 
  scale_y_discrete(labels = c('RG-neu','RG-gli','Cycling','IP','ExN','ExM','ExUp','ExDp','InMGE','InCGE','OPC','End','Per','VLMC','Mic'))

pdf(file = './res/step_101_fig_230207/harmony/mouse_predicted_cell_type_annotation_confusion_heatmap.pdf',width = 8,height = 5)
p
dev.off()
