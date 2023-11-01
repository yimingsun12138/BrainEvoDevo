#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque NPC lineage analysis                                    ##
## Data: 2023.08.04                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.1/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(dplyr)
library(Seurat)
library(ggplot2)
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
library(CellChat)
library(OpenAI4R)
library(paletteer)
library(ggpattern)
library(ggrepel)
library(openxlsx)
library(readxl)
library(monocle)
library(destiny)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0613')

# load and filter data ----------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'ReAnno_celltype',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

macaque_NPC_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% c('RG-1','IP','Ex-1','RG-2','OPC')]
Greenleaf_NPC_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('RG-1','IP','Ex-1','RG-2','OPC')]

# redo UMAP ---------------------------------------------------------------
Greenleaf_NPC_Seurat <- my_process_seurat(object = Greenleaf_NPC_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 1500,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
Greenleaf_NPC_Seurat <- my_process_seurat(object = Greenleaf_NPC_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 9,resolution = 0.3,group.by = 'ReAnno_celltype',label = TRUE)

var_gene <- intersect(x = VariableFeatures(Greenleaf_NPC_Seurat),y = rownames(macaque_NPC_Seurat@assays$RNA@counts))
Greenleaf_NPC_Seurat$species <- 'human'
macaque_NPC_Seurat$species <- 'macaque'

Brain_NPC_Seurat <- my_harmony_integration(named_seurat_list = list(human = Greenleaf_NPC_Seurat,macaque = macaque_NPC_Seurat),
                                           assay = 'RNA',variable_feature = var_gene,var_to_regress_list = NULL,npcs = 50,
                                           reference_loading = 'human',integration_var = 'species',harmony_input_dim = 10,max.iter.harmony = 10,
                                           UMAP_dim = 10,resolution = 0.5)
Brain_NPC_Seurat$cell_type <- NA
Brain_NPC_Seurat@meta.data[colnames(Greenleaf_NPC_Seurat),"cell_type"] <- Greenleaf_NPC_Seurat$ReAnno_celltype
Brain_NPC_Seurat@meta.data[colnames(macaque_NPC_Seurat),"cell_type"] <- macaque_NPC_Seurat$cell_type

DimPlot(object = Brain_NPC_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)
DimPlot(object = Brain_NPC_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')
DimPlot(object = Brain_NPC_Seurat,group.by = 'dataset',label = TRUE,repel = TRUE) + theme_bw() + theme(aspect.ratio = 1)

#I do not think this will give a good result...

# run monocle2 ------------------------------------------------------------
#monocle dataset
macaque_expression <- macaque_NPC_Seurat@assays$RNA@counts
macaque_metadata <- macaque_NPC_Seurat@meta.data
macaque_gene <- data.frame(gene = rownames(macaque_expression),gene_short_name = rownames(macaque_expression))
rownames(macaque_gene) <- macaque_gene$gene

macaque_metadata <- new('AnnotatedDataFrame',macaque_metadata)
macaque_gene <- new('AnnotatedDataFrame',macaque_gene)
macaque_NPC_monocle <- newCellDataSet(cellData = macaque_expression,phenoData = macaque_metadata,featureData = macaque_gene,expressionFamily = negbinomial.size())

#Estimate size factors and dispersions
macaque_NPC_monocle <- estimateSizeFactors(object = macaque_NPC_monocle)
macaque_NPC_monocle <- estimateDispersions(object = macaque_NPC_monocle)

#find variable gene
diff_test_res <- differentialGeneTest(macaque_NPC_monocle,fullModelFormulaStr = "~ cell_type")
diff_test_res <- diff_test_res[order(diff_test_res$qval,decreasing = FALSE),]
ordering_genes <- row.names(diff_test_res)[1:800]

macaque_NPC_monocle <- setOrderingFilter(macaque_NPC_monocle, ordering_genes)
plot_ordering_genes(macaque_NPC_monocle)

#DDRtree
macaque_NPC_monocle <- reduceDimension(macaque_NPC_monocle,max_components = 2,method = 'DDRTree')

#order cells
macaque_NPC_monocle <- orderCells(macaque_NPC_monocle)
plot_cell_trajectory(macaque_NPC_monocle, color_by = "cell_type")

# RG-2 OPC lineage --------------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
macaque_OPC_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% c('RG-2','OPC')]

#create monocle object
macaque_expression <- macaque_OPC_Seurat@assays$RNA@counts
macaque_metadata <- macaque_OPC_Seurat@meta.data
macaque_gene <- data.frame(gene = rownames(macaque_expression),gene_short_name = rownames(macaque_expression))
rownames(macaque_gene) <- macaque_gene$gene

macaque_metadata <- new('AnnotatedDataFrame',macaque_metadata)
macaque_gene <- new('AnnotatedDataFrame',macaque_gene)
macaque_OPC_monocle <- newCellDataSet(cellData = macaque_expression,phenoData = macaque_metadata,featureData = macaque_gene,expressionFamily = negbinomial.size())

#Estimate size factors and dispersions
macaque_OPC_monocle <- estimateSizeFactors(object = macaque_OPC_monocle)
macaque_OPC_monocle <- estimateDispersions(object = macaque_OPC_monocle)

#find variable gene
diff_test_res <- differentialGeneTest(macaque_OPC_monocle,fullModelFormulaStr = "~ cell_type")
diff_test_res <- diff_test_res[order(diff_test_res$qval,decreasing = FALSE),]
ordering_genes <- row.names(diff_test_res)[1:50]

macaque_OPC_monocle <- setOrderingFilter(macaque_OPC_monocle, ordering_genes)
plot_ordering_genes(macaque_OPC_monocle)

#DDRtree
macaque_OPC_monocle <- reduceDimension(macaque_OPC_monocle,max_components = 2,method = 'DDRTree')

#order cells
macaque_OPC_monocle <- orderCells(macaque_OPC_monocle,reverse = TRUE)
plot_cell_trajectory(macaque_OPC_monocle, color_by = "cell_type")
plot_cell_trajectory(macaque_OPC_monocle, color_by = "Pseudotime")

#Finding Genes that Change as a Function of Pseudotime
diff_test_res <- differentialGeneTest(macaque_OPC_monocle,fullModelFormulaStr = "~ sm.ns(Pseudotime)")
