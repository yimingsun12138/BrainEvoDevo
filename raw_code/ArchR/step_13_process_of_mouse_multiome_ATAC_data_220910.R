#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: process of mouse multiome ATAC data                             ##
## Data: 2022.09.10                                                                ##
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
library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(org.Mm.eg.db)
library(clusterProfiler)
library(viridis)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmusculus.UCSC.mm10)

geneAnnotation <- rtracklayer::import('/data/User/sunym/project/Brain/data/reference/ensembl_gtf_for_mapping/Mus_musculus.GRCm38.102.gtf')
geneAnnotation <- rtracklayer::as.data.frame(geneAnnotation)
geneAnnotation$seqnames <- as.character(geneAnnotation$seqnames)
geneAnnotation <- geneAnnotation[geneAnnotation$seqnames %in% c(as.character(1:19),'X','Y','MT'),]
geneAnnotation$seqnames <- paste0('chr',geneAnnotation$seqnames)
geneAnnotation <- GenomicRanges::makeGRangesFromDataFrame(geneAnnotation,keep.extra.columns = TRUE)

TxDb_object <- GenomicFeatures::makeTxDbFromGRanges(geneAnnotation)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb_object, OrgDb = org.Mm.eg.db)

# create arrow file -------------------------------------------------------
inputFiles <- list.files(path = './fragment_file/mouse_multiome_data/220907/')
temp <- sub(pattern = '_fragments.tsv.gz',replacement = '',x = inputFiles,fixed = TRUE)
inputFiles <- paste('./fragment_file/mouse_multiome_data/220907',inputFiles,sep = '/')
names(inputFiles) <- temp
file.exists(inputFiles)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  TileMatParams = list(tileSize = 5000),
  subThreading = FALSE
)

my_send_sms('create arrow file done!')

# doublet inference -------------------------------------------------------
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# create ArchR project from arrow file ------------------------------------
ArrowFiles <- list.files(path = './arrow_file/mouse_multiome_data/220910/')
ArrowFiles <- paste('./arrow_file/mouse_multiome_data/220910',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

mouse_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = './processed_data/mouse_multiome_ArchR_220910',
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

#add meta data
mouse_multiome_ArchR$species <- 'mouse'
mouse_multiome_ArchR$donor <- mouse_multiome_ArchR$Sample
mouse_multiome_ArchR$Age <- base::unlist(base::lapply(X = mouse_multiome_ArchR$donor,FUN = function(x){
  temp <- strsplit(x = x,split = '_')
  return(temp[[1]][1])
}))
mouse_multiome_ArchR$batch <- NA
mouse_multiome_ArchR@cellColData[mouse_multiome_ArchR$donor %in% c('E145_1','E155_1','E155_2'),"batch"] <- '220707'
mouse_multiome_ArchR@cellColData[mouse_multiome_ArchR$donor %in% c('E145_2','E155_3','E165_2'),"batch"] <- '220907'

# QC plot -----------------------------------------------------------------
#Tss enrichment vs. nFrags
df <- mouse_multiome_ArchR@cellColData[,c('nFrags','TSSEnrichment','Sample')]
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

pdf(file = './res/step_13_fig_220910/TSS_enrichment_vs_nFrag_dotplot.pdf',width = 12,height = 10)
QC_E145_1+QC_E145_2+QC_E155_1+QC_E155_2+QC_E155_3+QC_E165_2+plot_layout(ncol = 3)
dev.off()

#tss enrichment
p <- plotGroups(ArchRProj = mouse_multiome_ArchR, 
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

pdf(file = './res/step_13_fig_220910/TSS_enrichment_vlnplot.pdf',width = 4,height = 4)
p
dev.off()

#nFrag distribution
p <- plotGroups(ArchRProj = mouse_multiome_ArchR, 
                groupBy = "Sample", 
                colorBy = "cellColData", 
                name = "log10(nFrags)",
                plotAs = "violin", 
                alpha = 0.4, 
                addBoxPlot = TRUE) + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12,face = 'bold')) + xlab('')

pdf(file = './res/step_13_fig_220910/nFrag_vlnplot_.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_13_fig_220910/frag_length_.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_13_fig_220910/TSS_accessibility_by_sample.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220910/", load = FALSE, overwrite = TRUE)

# dim reduction -----------------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220910/')
#LSI
mouse_multiome_ArchR <- addIterativeLSI(ArchRProj = mouse_multiome_ArchR,
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
mouse_multiome_ArchR <- addClusters(
  input = mouse_multiome_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

cM <- confusionMatrix(i = mouse_multiome_ArchR$Clusters,j = mouse_multiome_ArchR$Sample)
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

#umap
mouse_multiome_ArchR <- addUMAP(
  ArchRProj = mouse_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:30, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

# InCGE? ------------------------------------------------------------------
mouse_multiome_ArchR <- addImputeWeights(ArchRProj = mouse_multiome_ArchR)

marker_gene <- c('Dlx5','Dlx2','Gad1','Gad2','Lhx6','Meis2','Sp8','Cxcl14','Htr3a','Adarb2','Neurod6')

p <- plotEmbedding(
  ArchRProj = mouse_multiome_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(mouse_multiome_ArchR)
)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "donor", embedding = "UMAP")

p1+p2+p$Dlx5+plot_layout(ncol = 3)
p1+p2+p$Dlx2+plot_layout(ncol = 3)
p1+p2+p$Gad1+plot_layout(ncol = 3)
p1+p2+p$Gad2+plot_layout(ncol = 3)
p1+p2+p$Meis2+plot_layout(ncol = 3)
p1+p2+p$Sp8+plot_layout(ncol = 3)
p1+p2+p$Adarb2+plot_layout(ncol = 3)
p1+p2+p$Neurod6+plot_layout(ncol = 3)

#hard to tell whether InCGE exists

#load meta data
meta_data <- readRDS(file = '/data/User/sunym/project/Brain/res/step_55_fig_220907/mouse_multiome_meta_data_220910_0038_temp_data.rds')
table(rownames(meta_data) %in% rownames(mouse_multiome_ArchR@cellColData))
mouse_multiome_ArchR$filted <- 'filted'
mouse_multiome_ArchR@cellColData[rownames(mouse_multiome_ArchR@cellColData) %in% rownames(meta_data),"filted"] <- 'preserved'

plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "filted", embedding = "UMAP")

# save meta data ----------------------------------------------------------
meta_data <- mouse_multiome_ArchR@cellColData
meta_data <- as.data.frame(meta_data)
saveRDS(object = meta_data,file = './res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')

# after filter ------------------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220910/')
meta_data <- readRDS(file = '/data/User/sunym/project/Brain/res/step_55_fig_220907/mouse_multiome_meta_data_220911_1211_temp_data.rds')
mouse_multiome_ArchR <- mouse_multiome_ArchR[rownames(meta_data)]

# dim reduction -----------------------------------------------------------
#LSI
mouse_multiome_ArchR <- addIterativeLSI(ArchRProj = mouse_multiome_ArchR,
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
mouse_multiome_ArchR <- addClusters(
  input = mouse_multiome_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:20, 
  force = TRUE
)

cM <- confusionMatrix(i = mouse_multiome_ArchR$Clusters,j = mouse_multiome_ArchR$Sample)
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

#umap
mouse_multiome_ArchR <- addUMAP(
  ArchRProj = mouse_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:20, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

#add meta data
mouse_multiome_ArchR$macaque_cell_type <- meta_data[rownames(mouse_multiome_ArchR@cellColData),"macaque_cell_type"]
p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "macaque_cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

#cluster 2 cell type
table(mouse_multiome_ArchR@cellColData[mouse_multiome_ArchR$Clusters == 'C2',"macaque_cell_type"])
#mainly cycling cells

# load mouse multiome seurat ----------------------------------------------
#load data
E145_1 <- Seurat::Read10X(data.dir = '../data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix/')
E145_2 <- Seurat::Read10X(data.dir = '../data/mouse_multiome/220907/E145_2/outs/filtered_feature_bc_matrix/')
E155_1 <- Seurat::Read10X(data.dir = '../data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix/')
E155_2 <- Seurat::Read10X(data.dir = '../data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix/')
E155_3 <- Seurat::Read10X(data.dir = '../data/mouse_multiome/220907/E155_3/outs/filtered_feature_bc_matrix/')
E165_2 <- Seurat::Read10X(data.dir = '../data/mouse_multiome/220907/E165_2/outs/filtered_feature_bc_matrix/')

E145_1 <- E145_1$`Gene Expression`
E145_2 <- E145_2$`Gene Expression`
E155_1 <- E155_1$`Gene Expression`
E155_2 <- E155_2$`Gene Expression`
E155_3 <- E155_3$`Gene Expression`
E165_2 <- E165_2$`Gene Expression`

colnames(E145_1) <- paste('E145_1#',colnames(E145_1),sep = '')
colnames(E145_2) <- paste('E145_2#',colnames(E145_2),sep = '')
colnames(E155_1) <- paste('E155_1#',colnames(E155_1),sep = '')
colnames(E155_2) <- paste('E155_2#',colnames(E155_2),sep = '')
colnames(E155_3) <- paste('E155_3#',colnames(E155_3),sep = '')
colnames(E165_2) <- paste('E165_2#',colnames(E165_2),sep = '')

E145_1 <- CreateSeuratObject(counts = E145_1,project = 'E145_1',assay = 'RNA',min.cells = 0,min.features = 0)
E145_2 <- CreateSeuratObject(counts = E145_2,project = 'E145_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_1 <- CreateSeuratObject(counts = E155_1,project = 'E155_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_2 <- CreateSeuratObject(counts = E155_2,project = 'E155_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_3 <- CreateSeuratObject(counts = E155_3,project = 'E155_3',assay = 'RNA',min.cells = 0,min.features = 0)
E165_2 <- CreateSeuratObject(counts = E165_2,project = 'E165_2',assay = 'RNA',min.cells = 0,min.features = 0)

MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
E145_1[['percent.mt']] <- PercentageFeatureSet(object = E145_1,features = MT_gene_list)
E145_2[['percent.mt']] <- PercentageFeatureSet(object = E145_2,features = MT_gene_list)
E155_1[['percent.mt']] <- PercentageFeatureSet(object = E155_1,features = MT_gene_list)
E155_2[['percent.mt']] <- PercentageFeatureSet(object = E155_2,features = MT_gene_list)
E155_3[['percent.mt']] <- PercentageFeatureSet(object = E155_3,features = MT_gene_list)
E165_2[['percent.mt']] <- PercentageFeatureSet(object = E165_2,features = MT_gene_list)

mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

temp <- unlist(base::lapply(X = colnames(mouse_multiome_Seurat),FUN = function(x){
  temp <- strsplit(x = x,split = '#',fixed = TRUE)
  temp <- temp[[1]][1]
  return(temp)
}))
mouse_multiome_Seurat$donor <- temp

temp <- unlist(base::lapply(X = mouse_multiome_Seurat$donor,FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  temp <- temp[[1]][1]
  return(temp)
}))
mouse_multiome_Seurat$Age <- temp

mouse_multiome_Seurat$batch <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_1','E155_1','E155_2'),"batch"] <- '220707'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_2','E155_3','E165_2'),"batch"] <- '220907'

#load meta data
meta_data <- readRDS(file = '../res/step_55_fig_220907/mouse_multiome_meta_data_220911_1211_temp_data.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]
mouse_multiome_Seurat$human_area <- meta_data$human_area
mouse_multiome_Seurat$macaque_cell_type <- meta_data$macaque_cell_type
mouse_multiome_Seurat$ArchR_Cluster <- NA
mouse_multiome_Seurat@meta.data[rownames(mouse_multiome_ArchR@cellColData),"ArchR_Cluster"] <- mouse_multiome_ArchR$Clusters

#save memory
remove(list = c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2'))
gc()

#pre-process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 20
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1),group.by = 'RNA_snn_res.0.5',label = TRUE)

#display
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)
p4 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'ArchR_Cluster',label = TRUE,repel = TRUE)
p1+p2+p3+p4+plot_layout(ncol = 2)

#show c2
DimPlot(object = mouse_multiome_Seurat,group.by = 'ArchR_Cluster',label = TRUE,repel = TRUE,cells.highlight = colnames(mouse_multiome_Seurat)[mouse_multiome_Seurat$ArchR_Cluster == 'C2'])
#c2 can be filted

# filter C2 ---------------------------------------------------------------
mouse_multiome_Seurat <- mouse_multiome_Seurat[,colnames(mouse_multiome_Seurat)[!(mouse_multiome_Seurat$ArchR_Cluster %in% c('C2'))]]
mouse_multiome_ArchR <- mouse_multiome_ArchR[colnames(mouse_multiome_Seurat)]

#redo process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = c('nCount_RNA','batch','Age','donor'),npcs = 50,preprocess = TRUE)
ndims <- 20
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1),group.by = 'RNA_snn_res.0.5',label = TRUE)

mouse_multiome_ArchR <- addIterativeLSI(ArchRProj = mouse_multiome_ArchR,
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

mouse_multiome_ArchR <- addClusters(
  input = mouse_multiome_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:20, 
  force = TRUE
)

cM <- confusionMatrix(i = mouse_multiome_ArchR$Clusters,j = mouse_multiome_ArchR$Sample)
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

mouse_multiome_ArchR <- addUMAP(
  ArchRProj = mouse_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:20, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "macaque_cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

#save meta data
meta_data <- as.data.frame(mouse_multiome_ArchR@cellColData)
saveRDS(object = meta_data,file = './res/step_13_fig_220910/mouse_multiome_ArchR_filted_meta_data_220912_1420.rds')
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220910/", load = FALSE, overwrite = TRUE)

# finally check mouse multiome Seurat -------------------------------------
mouse_multiome_Seurat@meta.data[rownames(mouse_multiome_ArchR@cellColData),"ArchR_Cluster"] <- mouse_multiome_ArchR$Clusters
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)
p4 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'ArchR_Cluster',label = TRUE,repel = TRUE)
p1+p2+p3+p4+plot_layout(ncol = 2)

#convert gene id
mouse_to_human <- read.csv(file = '../data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human <- mouse_to_human[,c(2,1,4,3)]
temp <- mouse_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)
mouse_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

# load data
macaque_multiome_Seurat <- readRDS(file = '../processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')
gene_list <- dplyr::intersect(x = VariableFeatures(macaque_multiome_Seurat),y = rownames(mouse_multiome_Seurat@assays$converted@counts))
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'PCA',preprocess = FALSE,dim_to_use = 32,resolution = 1,group.by = 'cell_type',label = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = gene_list,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

#project
temp <- projectMatrix_SeuratUMAP(X_scaled = mouse_multiome_Seurat@assays$converted@scale.data,object = macaque_multiome_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
mouse_multiome_Seurat[['projected_PCA']] <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'converted')
mouse_multiome_Seurat[['projected_UMAP']] <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'converted')

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = mouse_multiome_Seurat,ref_reduction = 'PCA',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'converted',l2.norm = TRUE,dims = 1:32,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$cell_type,l2.norm = TRUE,dims = 1:32,verbose = TRUE)
mouse_multiome_Seurat <- AddMetaData(object = mouse_multiome_Seurat,metadata = predictions)
mouse_multiome_Seurat$macaque_cell_type <- mouse_multiome_Seurat$predicted.id
mouse_multiome_Seurat$macaque_cell_type_score <- mouse_multiome_Seurat$prediction.score.max

#display
DimPlot(object = mouse_multiome_Seurat,reduction = 'umap',group.by = 'macaque_cell_type',label = TRUE,repel = TRUE)
FeaturePlot(object = mouse_multiome_Seurat,features = 'EOMES',reduction = 'umap')

#save data
saveRDS(object = mouse_multiome_Seurat,file = '/data/User/sunym/temp/mouse_multiome_Seurat_220912.rds')

# update ArchR annotation -------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220910/')
mouse_multiome_ArchR$macaque_cell_type <- mouse_multiome_Seurat@meta.data[rownames(mouse_multiome_ArchR@cellColData),"macaque_cell_type"]
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'donor',label = FALSE,reduction = 'umap')
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE,reduction = 'umap')
p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE,reduction = 'umap')
p4 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'ArchR_Cluster',label = TRUE,repel = TRUE,reduction = 'umap')
p1+p2+p3+p4+plot_layout(ncol = 2)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "macaque_cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220910/", load = FALSE, overwrite = TRUE)
