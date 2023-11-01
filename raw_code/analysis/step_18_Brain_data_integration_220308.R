#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Brain data integration                                          ##
## Data: 2022.03.08                                                                ##
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
setwd('/data/User/sunym/project/Brain/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')

#library
library(Rmisc)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(sctransform)
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
library(patchwork)
library(networkD3)
library(htmlwidgets)
library(circlize)
library(harmony)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# integrate with Greenleaf RNA-seq ----------------------------------------
macaque_integration_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
Brain_RNA_Seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_greenleaf_RNA_harmony_MNN_seurat_211102.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
table(colnames(Greenleaf_RNA_Seurat) %in% colnames(Brain_RNA_Seurat))

Brain_RNA_Seurat <- Brain_RNA_Seurat[,!(Brain_RNA_Seurat$predict_label %in% c('Unknown'))]
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,colnames(Greenleaf_RNA_Seurat) %in% colnames(Brain_RNA_Seurat)]
Greenleaf_RNA_Seurat$predict_label <- Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"predict_label"]
#convert gene id
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')

temp <- macaque_integration_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 40*(1024^3),workers = 5)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_integration_Seurat[['converted']] <- temp

#find variable feature
macaque_integration_Seurat <- my_process_seurat(object = macaque_integration_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(VariableFeatures(macaque_integration_Seurat),rownames(Greenleaf_RNA_Seurat@assays$converted@counts))

macaque_integration_Seurat$species <- 'macaque'
Greenleaf_RNA_Seurat$species <- 'human'

#harmony
Brain_RNA_Seurat <- my_harmony_integration(named_seurat_list = list(macaque=macaque_integration_Seurat,human=Greenleaf_RNA_Seurat),assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque=c('nCount_converted','donor','batch','tech'),human=c('Age','Batch','nCount_converted')),npcs = 50,reference_loading = 'macaque',
                                           integration_var = 'species',harmony_input_dim = 30,max.iter.harmony = 50,UMAP_dim = 30,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                           yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')
my_send_sms('harmony done!')

Brain_RNA_Seurat$cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(macaque_integration_Seurat),"cell_type"] <- macaque_integration_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"cell_type"] <- Greenleaf_RNA_Seurat$predict_label
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')

#label transfer
temp_human <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'human']
temp_macaque <- Brain_RNA_Seurat[,Brain_RNA_Seurat$dataset == 'macaque']
anchors <- my_FindTransferAnchors(reference = temp_macaque,query = temp_human,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:30,features = gene_list,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = temp_macaque$cell_type,weight.reduction = 'pcaproject',l2.norm = FALSE,dims = 1:30,verbose = TRUE)
temp_human <- AddMetaData(object = temp_human,metadata = predictions)
DimPlot(temp_human,group.by = 'predicted.id',label = TRUE,repel = TRUE) + 
  DimPlot(temp_macaque,group.by = 'cell_type',label = TRUE,repel = TRUE)

scibet::Confusion_heatmap(ori = temp_human$cell_type,prd = temp_human$predicted.id)
Brain_RNA_Seurat@meta.data[colnames(temp_human),"cell_type"] <- temp_human$predicted.id

p1 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = Brain_RNA_Seurat,group.by = 'dataset',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
pdf(file = './res/step_18_fig_220308/macaque_integration_greenleaf_RNA_harmony_label_transfer_dimplot.pdf',width = 13,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#save data
saveRDS(Brain_RNA_Seurat,file = './processed_data/220305_summary/macaque_integration_greenleaf_RNA_harmony_label_transfer_seurat_220309.rds')
