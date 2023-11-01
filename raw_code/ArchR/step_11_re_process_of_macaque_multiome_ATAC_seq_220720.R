#####################################################################################
## Project: macaque fetal Brain multiome ATAC_seq                                  ##
## Script Purpose: re_process of macaque multiome ATAC_seq                         ##
## Data: 2022.07.20                                                                ##
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
setwd('/data/User/sunym/project/Brain/ArchR/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(parallel)
library(ArchR)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(harmony)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(rtracklayer)
library(org.Mmu.eg.db)
library(clusterProfiler)
library(scibet)
library(viridis)
library(networkD3)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------

#load macaque multiome data
temp <- loadArchRProject(path = './processed_data/macaque_multiome_ArchR_220411/')

#load Arrow file
ArrowFiles <- list.files(path = './arrow_file/macaque_multiome_data/')
ArrowFiles <- paste('./arrow_file/macaque_multiome_data',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

#create ArchR project
macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/project/Brain/processed_data/220718_summary/macaque_multiome_ArchR_220720',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

#add metadata
macaque_multiome_ArchR@cellColData$species <- 'macaque'
macaque_multiome_ArchR@cellColData$donor <- as.character(macaque_multiome_ArchR$Sample)
table(macaque_multiome_ArchR$donor)

#filter cells
temp <- readRDS(file = '../processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')
table(colnames(temp) %in% rownames(macaque_multiome_ArchR@cellColData))
macaque_multiome_ArchR <- macaque_multiome_ArchR[colnames(temp)]

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# QC plot -----------------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')

#QC plot
df <- macaque_multiome_ArchR@cellColData[,c('nFrags','TSSEnrichment','Sample')]
df$nFrags <- log10(df$nFrags)
colnames(df) <- replace(x = colnames(df),list = c(colnames(df) == 'nFrags'),values = 'log10(nFrags)')
df

for (i in unique(df$Sample)) {
  p <- df[df$Sample == i,]
  assign(x = paste('QC',as.character(i),sep = '_'),
         value = ggPoint(
           x = p[,1], 
           y = p[,2], 
           colorDensity = TRUE,
           continuousSet = "sambaNight",
           xlabel = "Log10 Unique Fragments",
           ylabel = "TSS Enrichment",
           title = as.character(i),
           xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
           ylim = c(0, quantile(df[,2], probs = 0.99))) + 
           theme(plot.title = element_text(hjust = 0.5)) + 
           geom_hline(yintercept = 4, lty = "dashed") + 
           geom_vline(xintercept = 3, lty = "dashed"))
}

pdf(file = './res/step_11_fig_220720/TSS_enrichment_vs_nFrag_dotplot_after_filter.pdf',width = 16,height = 5)
QC_A50A+QC_A50B+QC_A82A+QC_A82B+plot_layout(ncol = 4)
dev.off()

#tss enrichment
p <- plotGroups(ArchRProj = macaque_multiome_ArchR, 
                groupBy = "Sample", 
                colorBy = "cellColData", 
                name = "TSSEnrichment",
                plotAs = "violin", 
                alpha = 0.4, 
                addBoxPlot = TRUE) + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12,face = 'bold')) + xlab('') + 
  scale_y_continuous(breaks = seq(0,30,5))

pdf(file = './res/step_11_fig_220720/TSS_enrichment_vlnplot_after_filter.pdf',width = 5,height = 5)
p
dev.off()

#nFrag distribution
p <- plotGroups(ArchRProj = macaque_multiome_ArchR, 
                groupBy = "Sample", 
                colorBy = "cellColData", 
                name = "log10(nFrags)",
                plotAs = "violin", 
                alpha = 0.4, 
                addBoxPlot = TRUE) + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12,face = 'bold')) + xlab('')

pdf(file = './res/step_11_fig_220720/nFrag_vlnplot_after_filter.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_11_fig_220720/frag_length_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_11_fig_220720/TSS_accessibility_by_sample_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

# dim reduction -----------------------------------------------------------
#LSI
macaque_multiome_ArchR <- addIterativeLSI(ArchRProj = macaque_multiome_ArchR,
                                          useMatrix = "TileMatrix", 
                                          name = "IterativeLSI", 
                                          iterations = 3, 
                                          clusterParams = list(
                                            resolution = c(0.6), 
                                            sampleCells = 10000, 
                                            n.start = 10,
                                            maxClusters = 30
                                          ), 
                                          varFeatures = 20000, 
                                          totalFeatures = 500000, 
                                          dimsToUse = 1:30, 
                                          corCutOff = 0.5, 
                                          force = TRUE)

#cluster
macaque_multiome_ArchR <- addClusters(
  input = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:15, 
  force = TRUE
)

cM <- confusionMatrix(i = macaque_multiome_ArchR$Clusters,j = macaque_multiome_ArchR$Sample)
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

#umap
macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:15, 
  force = TRUE
)

#add cell type
temp <- readRDS(file = '../processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')
macaque_multiome_ArchR$cell_type <- temp@meta.data[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]
macaque_multiome_ArchR$sub_cell_type <- temp@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub_cell_type"]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
pdf(file = './res/step_11_fig_220720/macaque_multiome_ArchR_cell_type_dimplot.pdf',width = 18,height = 8)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "sub_cell_type", embedding = "UMAP")
pdf(file = './res/step_11_fig_220720/macaque_multiome_ArchR_sub_cell_type_dimplot.pdf',width = 18,height = 8)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

scibet::Confusion_heatmap(ori = macaque_multiome_ArchR$Clusters,prd = macaque_multiome_ArchR$cell_type)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# label transfer ----------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')
macaque_multiome_Seurat <- readRDS(file = '../processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')

DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#constrained integration
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('End','Mic','OPC','Per','VLMC')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C1','C2','C3')]
  ),
  npc = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('Cycling','IP','RG-1','RG-2')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C4','C5','C6','C7')]
  ),
  Ex_neuron = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C13','C14','C15','C16','C17','C18')]
  ),
  In_neuron = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C8','C9','C10','C11','C12')]
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

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

scibet::Confusion_heatmap(ori = macaque_multiome_ArchR$cell_type,prd = macaque_multiome_ArchR$predictedGroup)

#sankey plot
confusion_matrix <- my_confusion_matrix(ori = paste(macaque_multiome_ArchR$cell_type,'RNA',sep = '.'),prd = paste(macaque_multiome_ArchR$predictedGroup,'ATAC',sep = '.'))
confusion_matrix <- confusion_matrix[confusion_matrix$value > 20,]
nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

p <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                   Source = "IDori", Target = "IDprd",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, nodeWidth= 20, fontSize=13, nodePadding=20)
p

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# add gene expression matrix ----------------------------------------------
#load data

#ATAC
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')
#RNA
A50A_seRNA <- import10xFeatureMatrix(input = '/data/User/sunym/project/Brain/data/Multiome/A50A/outs/filtered_feature_bc_matrix.h5',
                                     names = 'A50A')
A50B_seRNA <- import10xFeatureMatrix(input = '/data/User/sunym/project/Brain/data/Multiome/A50B/outs/filtered_feature_bc_matrix.h5',
                                     names = 'A50B')
A82A_seRNA <- import10xFeatureMatrix(input = '/data/User/sunym/project/Brain/data/Multiome/A82A/outs/filtered_feature_bc_matrix.h5',
                                     names = 'A82A')
A82B_seRNA <- import10xFeatureMatrix(input = '/data/User/sunym/project/Brain/data/Multiome/A82B/outs/filtered_feature_bc_matrix.h5',
                                     names = 'A82B')

#modify rowRanges
temp <- rowRanges(A50A_seRNA)
name_list <- seqlevels(temp)
name_list <- base::lapply(X = name_list,FUN = function(x){
  if(x %in% c(as.character(1:20),'X','Y')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
})
name_list <- base::unlist(name_list)
names(name_list) <- seqlevels(temp)
temp <- renameSeqlevels(x = temp,value = name_list)

#modify expression matrix
gene_list <- rownames(assay(A50A_seRNA))
seRNA <- cbind(assay(A50A_seRNA)[gene_list,],assay(A50B_seRNA)[gene_list,],assay(A82A_seRNA)[gene_list,],assay(A82B_seRNA)[gene_list,])
seRNA <- seRNA[,rownames(macaque_multiome_ArchR@cellColData)]
seRNA <- SummarizedExperiment(assays = list(counts=seRNA),rowRanges = temp)

macaque_multiome_ArchR <- addGeneExpressionMatrix(input = macaque_multiome_ArchR,seRNA = seRNA,verbose = TRUE,force = TRUE)
getAvailableMatrices(macaque_multiome_ArchR)

temp <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneExpressionMatrix',verbose = TRUE)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# compare multiome and label transfer -------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')
macaque_multiome_Seurat <- readRDS(file = '../processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')

DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#add gex cell type
macaque_multiome_ArchR$Gex_cell_type <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]
table(macaque_multiome_ArchR$Gex_cell_type == macaque_multiome_ArchR$cell_type)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
pdf(file = './res/step_11_fig_220720/macaque_multiome_RNA_cell_type_label_transfer_dimplot.pdf',width = 16,height = 7)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

temp <- function(ori,prd){
  cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
  cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
  return(cross.validation.filt)
}

temp <- temp(ori = macaque_multiome_ArchR$Gex_cell_type,prd = macaque_multiome_ArchR$predictedGroup)
temp$ori <- factor(temp$ori,levels = rev(c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-2','RG-1','OPC','Mic','End','Per','VLMC')))
temp$prd <- factor(temp$prd,levels = rev(c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-2','RG-1','OPC','Mic','End','Per','VLMC')))

pdf(file = './res/step_11_fig_220720/macaque_multiome_RNA_cell_type_label_transfer_confusion_heatmap.pdf',width = 5.5,height = 5.5)
ggplot(data = temp,aes(x=ori,y=prd,fill=Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), 
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_viridis() + 
  theme(aspect.ratio = 1,axis.title = element_text(size = 12,face = 'bold')) + 
  xlab('Gex_cell_type') + ylab('predict label')
dev.off()

#add gex sub cell type
DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
macaque_multiome_ArchR$Gex_sub_cell_type <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub_cell_type"]
table(macaque_multiome_ArchR$sub_cell_type == macaque_multiome_ArchR$Gex_sub_cell_type)

scibet::Confusion_heatmap(ori = macaque_multiome_ArchR$Gex_cell_type,prd = macaque_multiome_ArchR$Gex_sub_cell_type)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_sub_cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# cell type correlation ---------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')
macaque_multiome_Seurat <- readRDS(file = '../processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')

#gene express matrix
gene_express_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneExpressionMatrix',verbose = TRUE)
temp <- gene_express_matrix@colData
gene_list <- gene_express_matrix@elementMetadata$name
gene_express_matrix <- gene_express_matrix@assays@data$GeneExpressionMatrix
rownames(gene_express_matrix) <- gene_list
gene_express_matrix <- CreateSeuratObject(counts = gene_express_matrix,project = 'temp',assay = 'RNA',meta.data = as.data.frame(temp),min.cells = 0,min.features = 0)

#gene integratin matrix
gene_integration_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneIntegrationMatrix',verbose = TRUE)
temp <- as.data.frame(gene_integration_matrix@colData)
gene_list <- gene_integration_matrix@elementMetadata$name
gene_list <- gene_list[duplicated(gene_list)]
gene_list <- unique(gene_list)
gene_list <- !(gene_integration_matrix@elementMetadata$name %in% gene_list)
gene_integration_matrix <- gene_integration_matrix[gene_list,]
table(duplicated(gene_integration_matrix@elementMetadata$name))
gene_list <- gene_integration_matrix@elementMetadata$name
gene_integration_matrix <- gene_integration_matrix@assays@data$GeneIntegrationMatrix
rownames(gene_integration_matrix) <- gene_list
gene_integration_matrix <- CreateSeuratObject(counts = gene_integration_matrix,project = 'temp',assay = 'RNA',meta.data = temp,min.cells = 0,min.features = 0)

#create pseudo-bulk matrix
gene_express_matrix <- AverageExpression(object = gene_express_matrix,assays = 'RNA',return.seurat = FALSE,group.by = 'Gex_cell_type',slot = 'counts',verbose = TRUE)
gene_express_matrix <- log1p(gene_express_matrix$RNA)
gene_integration_matrix <- AverageExpression(object = gene_integration_matrix,assays = 'RNA',return.seurat = FALSE,group.by = 'predictedGroup',slot = 'counts',verbose = TRUE)
gene_integration_matrix <- log1p(gene_integration_matrix$RNA)

#define gene list
gene_list <- VariableFeatures(macaque_multiome_Seurat)
gene_list <- dplyr::intersect(gene_list,rownames(gene_express_matrix))
gene_list <- dplyr::intersect(gene_list,rownames(gene_integration_matrix))

temp <- matrix(nrow = dim(gene_express_matrix)[2],ncol = dim(gene_integration_matrix)[2])
colnames(temp) <- colnames(gene_integration_matrix)
rownames(temp) <- colnames(gene_express_matrix)
temp <- as.data.frame(temp)
for (i in rownames(temp)) {
  for (j in colnames(temp)) {
    temp[i,j] <- cor(x = gene_express_matrix[gene_list,i],y = gene_integration_matrix[gene_list,j],method = 'pearson')
  }
}

cell_list_col <- c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-2','RG-1','OPC','Mic','End','Per','VLMC')
cell_list_row <- c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-2','RG-1','OPC','Mic','End','Per','VLMC')

pdf(file = './res/step_11_fig_220720/macaque_multiome_RNA_prediction_expression_cor_on_variable_gene.pdf',width = 6,height = 6)
Heatmap(matrix = temp[cell_list_row,cell_list_col],cluster_rows = FALSE,cluster_columns = FALSE,row_title = 'multiome express',column_title = 'prediction express',
        height = unit(4,'inches'),width = unit(4,'inches'),
        row_split = factor(c('In','In','Ex','Ex','Ex','Ex','NPC','NPC','NPC','NPC','NN','NN','NN','NN','NN'),levels = c('In','Ex','NPC','NN')),
        column_split = factor(c('In','In','Ex','Ex','Ex','Ex','NPC','NPC','NPC','NPC','NN','NN','NN','NN','NN'),levels = c('In','Ex','NPC','NN')),
        border = TRUE,name = 'Cor')
dev.off()

#re-fine cell type
macaque_multiome_ArchR$cell_type <- macaque_multiome_ArchR$Gex_cell_type
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Gex_cell_type %in% c('RG-1','RG-2'),"cell_type"] <- 'RG'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Gex_cell_type %in% c('End','Per','VLMC'),"cell_type"] <- 'End/Per/VLMC'

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_sub_cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

p1+p2+p3+p4+plot_layout(ncol = 4)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# marker gene validation --------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')

marker_gene <- c('PECAM1', #End
                 'PDGFRB', #Per
                 'CX3CR1', #Mic
                 'DNAH11','GJA1','SPATA17','LRRC71','SPAG17','FOXJ1', #Ependymal
                 'SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3', #RG
                 'SOX10','OLIG2','EGFR', #OPC
                 'TOP2A','MKI67','CLSPN','AURKA', #Cyc
                 'EOMES','PPP1R17', #IP
                 'NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2', #Ex
                 'DLX5','GAD2','GAD1','DLX2', #In
                 'LHX6','SST', #InMGE
                 'SP8','NR2F2', #InCGE
                 'MEIS2','ETV1' #InPSB
)

p <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_multiome_ArchR)
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

#End/Per checked
p1+p2+p$PECAM1+plot_layout(ncol = 3)
p1+p2+p$PDGFRB+plot_layout(ncol = 3)

#Mic checked
p1+p2+p$CX3CR1+plot_layout(ncol = 3)

#RG checked
p1+p2+p$SOX9+plot_layout(ncol = 3)
p1+p2+p$PAX6+plot_layout(ncol = 3)
p1+p2+p$VIM+plot_layout(ncol = 3)

#OPC checked
p1+p2+p$SOX10+plot_layout(ncol = 3)
p1+p2+p$OLIG2+plot_layout(ncol = 3)
p1+p2+p$EGFR+plot_layout(ncol = 3)

#Cycling checked
p1+p2+p$TOP2A+plot_layout(ncol = 3)
p1+p2+p$MKI67+plot_layout(ncol = 3)

#IP checked
p1+p2+p$EOMES+plot_layout(ncol = 3)
p1+p2+p$PPP1R17+plot_layout(ncol = 3)

#Ex checked
p1+p2+p$NEUROD2+plot_layout(ncol = 3)
p1+p2+p$NEUROD6+plot_layout(ncol = 3)
p1+p2+p$TBR1+plot_layout(ncol = 3)
p1+p2+p$SATB2+plot_layout(ncol = 3)
p1+p2+p$FEZF2+plot_layout(ncol = 3)

#In checked
p1+p2+p$GAD2+plot_layout(ncol = 3)
p1+p2+p$GAD1+plot_layout(ncol = 3)
p1+p2+p$DLX5+plot_layout(ncol = 3)
p1+p2+p$DLX2+plot_layout(ncol = 3)

#InMGE checked
p1+p2+p$LHX6+plot_layout(ncol = 3)
p1+p2+p$SST+plot_layout(ncol = 3)

#InCGE checked
p1+p2+p$SP8+plot_layout(ncol = 3)
p1+p2+p$NR2F2+plot_layout(ncol = 3)

#all checked

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# call peak ---------------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+p4+plot_layout(ncol = 4)

#notice point:
#maxpeaks = 150k
#macs2 FDR = 0.05
#modify mincells and max cells while doing pseudobulk
#mincells: 120

#add pseudo bulk
temp <- paste(macaque_multiome_ArchR$donor,macaque_multiome_ArchR$cell_type,sep = '_')
table(temp) %>% max()
table(temp) %>% min()
macaque_multiome_ArchR <- addGroupCoverages(ArchRProj = macaque_multiome_ArchR, 
                                            groupBy = "cell_type", 
                                            minCells = 120, 
                                            maxCells = max(table(temp)), 
                                            maxFragments = 10^10, 
                                            minReplicates = 4, 
                                            maxReplicates = 4, 
                                            sampleRatio = 0.8)

temp <- as.data.frame(macaque_multiome_ArchR@projectMetadata$GroupCoverages$cell_type$coverageMetadata)

#call peak
pathToMacs2 <- '/data/User/sunym/env/MACS2/bin/macs2'

macaque_multiome_ArchR <- addReproduciblePeakSet(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  pathToMacs2 = pathToMacs2,
  maxPeaks = 150000,
  genomeSize = 2.7e+09,
  cutOff = 0.05,
  force = TRUE
)

my_send_sms('call peak done!')

getPeakSet(macaque_multiome_ArchR)

#add peak matrix
macaque_multiome_ArchR <- addPeakMatrix(macaque_multiome_ArchR)
getAvailableMatrices(macaque_multiome_ArchR)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# peak co-accessibility ---------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')

#add co-accessibility
macaque_multiome_ArchR <- addCoAccessibility(
  ArchRProj = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI"
)

#see the results
cA <- getCoAccessibility(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA

markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(macaque_multiome_ArchR,resolution = 2000),
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','RG','OPC')
)

#RG
grid::grid.newpage()
grid::grid.draw(p$SOX9)

grid::grid.newpage()
grid::grid.draw(p$HOPX)

grid::grid.newpage()
grid::grid.draw(p$PAX6)

grid::grid.newpage()
grid::grid.draw(p$VIM)

#OPC
grid::grid.newpage()
grid::grid.draw(p$OLIG2)

grid::grid.newpage()
grid::grid.draw(p$EGFR)

#IP
grid::grid.newpage()
grid::grid.draw(p$EOMES)

grid::grid.newpage()
grid::grid.draw(p$PPP1R17)

#Ex
grid::grid.newpage()
grid::grid.draw(p$NEUROD2)

grid::grid.newpage()
grid::grid.draw(p$NEUROD6)

grid::grid.newpage()
grid::grid.draw(p$SATB2)

grid::grid.newpage()
grid::grid.draw(p$FEZF2)

#In
grid::grid.newpage()
grid::grid.draw(p$DLX5)

grid::grid.newpage()
grid::grid.draw(p$GAD2)

grid::grid.newpage()
grid::grid.draw(p$SST)

grid::grid.newpage()
grid::grid.draw(p$SP8)

grid::grid.newpage()
grid::grid.draw(p$LHX6)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)

# peak gene expression correlation ----------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')

#add P2G
macaque_multiome_ArchR <- addPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI",
  useMatrix = 'GeneExpressionMatrix'
)

p2g <- getPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g

markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','RG','OPC')
)

grid::grid.newpage()
grid::grid.draw(p$SOX9)

grid::grid.newpage()
grid::grid.draw(p$PAX6)

grid::grid.newpage()
grid::grid.draw(p$HOPX)

grid::grid.newpage()
grid::grid.draw(p$OLIG2)

grid::grid.newpage()
grid::grid.draw(p$EGFR)

grid::grid.newpage()
grid::grid.draw(p$EOMES)

grid::grid.newpage()
grid::grid.draw(p$NEUROD2)

grid::grid.newpage()
grid::grid.draw(p$NEUROD6)

grid::grid.newpage()
grid::grid.draw(p$DLX5)

p <- plotPeak2GeneHeatmap(ArchRProj = macaque_multiome_ArchR, groupBy = "cell_type", k = 10)
ComplexHeatmap::draw(p, heatmap_legend_side = "bot", annotation_legend_side = "bot", width = unit(2,'inches'), heigh = unit(8,'inches'))

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/', load = FALSE, overwrite = TRUE)
