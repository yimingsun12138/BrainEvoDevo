#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: process of mouse multiome ATAC data                             ##
## Data: 2022.10.03                                                                ##
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

# create ArchRProject from arrow file -------------------------------------
ArrowFiles <- list.files(path = './arrow_file/mouse_multiome_data/220910/')
ArrowFiles <- paste('./arrow_file/mouse_multiome_data/220910',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

mouse_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003',
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

# add RNA data ------------------------------------------------------------
#filter data
mouse_multiome_Seurat <- readRDS(file = '../processed_data/220802_summary/mouse_multiome_Seurat_220922.rds')
table(colnames(mouse_multiome_Seurat) %in% rownames(mouse_multiome_ArchR@cellColData))
mouse_multiome_ArchR <- mouse_multiome_ArchR[colnames(mouse_multiome_Seurat)]

#load RNA
E145_1_seRNA <- import10xFeatureMatrix(input = '../data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E145_1')
E145_2_seRNA <- import10xFeatureMatrix(input = '../data/mouse_multiome/220907/E145_2/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E145_2')
E155_1_seRNA <- import10xFeatureMatrix(input = '../data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E155_1')
E155_2_seRNA <- import10xFeatureMatrix(input = '../data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E155_2')
E155_3_seRNA <- import10xFeatureMatrix(input = '../data/mouse_multiome/220907/E155_3/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E155_3')
E165_2_seRNA <- import10xFeatureMatrix(input = '../data/mouse_multiome/220907/E165_2/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E165_2')

#modify rowRanges
temp <- rowRanges(E145_1_seRNA)
name_list <- seqlevels(temp)
name_list <- base::lapply(X = name_list,FUN = function(x){
  if(x %in% c(as.character(1:19),'X','Y')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
})
name_list <- base::unlist(name_list)
names(name_list) <- seqlevels(temp)
temp <- renameSeqlevels(x = temp,value = name_list)

#modify expression matrix
gene_list <- rownames(assay(E145_1_seRNA))
seRNA <- cbind(assay(E145_1_seRNA)[gene_list,],assay(E145_2_seRNA)[gene_list,],assay(E155_1_seRNA)[gene_list,],assay(E155_2_seRNA)[gene_list,],assay(E155_3_seRNA)[gene_list,],assay(E165_2_seRNA)[gene_list,])
seRNA <- seRNA[,rownames(mouse_multiome_ArchR@cellColData)]
seRNA <- SummarizedExperiment(assays = list(counts=seRNA),rowRanges = temp)

mouse_multiome_ArchR <- addGeneExpressionMatrix(input = mouse_multiome_ArchR,seRNA = seRNA,verbose = TRUE,force = TRUE)
getAvailableMatrices(mouse_multiome_ArchR)

temp <- getMatrixFromProject(ArchRProj = mouse_multiome_ArchR,useMatrix = 'GeneExpressionMatrix',verbose = TRUE)

#add meta data
mouse_multiome_ArchR$Gex_cell_type <- mouse_multiome_Seurat@meta.data[rownames(mouse_multiome_ArchR@cellColData),"cell_type"]
mouse_multiome_ArchR$Gex_macaque_cell_type <- mouse_multiome_Seurat@meta.data[rownames(mouse_multiome_ArchR@cellColData),"macaque_cell_type"]

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003", load = FALSE, overwrite = TRUE)

# QC plot -----------------------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003')
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

pdf(file = './res/step_15_fig_221003/TSS_enrichment_vs_nFrag_dotplot.pdf',width = 12,height = 10)
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

pdf(file = './res/step_15_fig_221003/TSS_enrichment_vlnplot.pdf',width = 4,height = 4)
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

pdf(file = './res/step_15_fig_221003/nFrag_vlnplot_.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_15_fig_221003/frag_length_.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_15_fig_221003/TSS_accessibility_by_sample.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003", load = FALSE, overwrite = TRUE)

# dim reduction -----------------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003')
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
  dimsToUse = 1:22, 
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
  dimsToUse = 1:22, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Gex_macaque_cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = 'cellColData', name = 'Gex_cell_type', embedding = 'UMAP')
p1+p2+p3+plot_layout(ncol = 3)

#cluster C11 seems weird
p1 <- my_dimplot(embedding = mouse_multiome_ArchR@embeddings$UMAP$df,meta_data = mouse_multiome_ArchR@cellColData,group.by = 'Clusters',cells.highlight = rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Clusters == 'C11'])
p2 <- my_dimplot(embedding = mouse_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = mouse_multiome_Seurat@meta.data,group.by = 'cell_type',cells.highlight = rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Clusters == 'C11'])
p1+p2+plot_layout(ncol = 2)

#filter C11
mouse_multiome_ArchR <- mouse_multiome_ArchR[!(mouse_multiome_ArchR$Clusters == 'C11')]

#redo dim reduction
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
  dimsToUse = 1:16, 
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
  dimsToUse = 1:16, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Gex_macaque_cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = 'cellColData', name = 'Gex_cell_type', embedding = 'UMAP')
p1+p2+p3+plot_layout(ncol = 3)

#consistancy
scibet::Confusion_heatmap(ori = mouse_multiome_ArchR$Clusters,prd = mouse_multiome_ArchR$Gex_cell_type)
scibet::Confusion_heatmap(ori = mouse_multiome_ArchR$Clusters,prd = mouse_multiome_ArchR$Gex_macaque_cell_type)

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003", load = FALSE, overwrite = TRUE)

# label transfer ----------------------------------------------------------
#load data
mouse_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003')
mouse_multiome_Seurat <- readRDS(file = '../processed_data/220802_summary/mouse_multiome_Seurat_220922.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(mouse_multiome_ArchR@cellColData)]

DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#constrained integration
p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Gex_macaque_cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(mouse_multiome_Seurat)[mouse_multiome_Seurat$macaque_cell_type %in% c('End','Mic','Per','VLMC')],
    ATAC = rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Clusters %in% c('C16','C17','C18')]
  ),
  In = SimpleList(
    RNA = colnames(mouse_multiome_Seurat)[mouse_multiome_Seurat$macaque_cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(mouse_multiome_ArchR@cellColData)[mouse_multiome_ArchR$Clusters %in% c('C1','C2')]
  ),
  Ex = SimpleList(
    RNA = colnames(mouse_multiome_Seurat)[mouse_multiome_Seurat$macaque_cell_type %in% c('Cycling','Ex-1','Ex-2','Ex-3','Ex-4','IP','RG-1')],
    ATAC = rownames(mouse_multiome_ArchR@cellColData)[!(mouse_multiome_ArchR$Clusters %in% c('C1','C2','C16','C17','C18'))]
  )
)

mouse_multiome_ArchR <- addImputeWeights(mouse_multiome_ArchR)
mouse_multiome_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = mouse_multiome_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = mouse_multiome_Seurat,
  addToArrow = TRUE, 
  groupList = group_list,
  groupRNA = "macaque_cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Gex_macaque_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

scibet::Confusion_heatmap(ori = mouse_multiome_ArchR$Gex_macaque_cell_type,prd = mouse_multiome_ArchR$predictedGroup)

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003", load = FALSE, overwrite = TRUE)

# call peak ---------------------------------------------------------------
#load data
mouse_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003')

p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Gex_macaque_cell_type", embedding = "UMAP")
p1+p2+p3+p4+plot_layout(ncol = 4)

#notice point:
#maxpeaks = 150k
#macs2 FDR = 0.05
#modify mincells and max cells while doing pseudobulk
#mincells: 40

#add pseudo bulk
temp <- paste(mouse_multiome_ArchR$donor,mouse_multiome_ArchR$Gex_macaque_cell_type,sep = '_')
table(temp) %>% max()
table(temp) %>% min()
mouse_multiome_ArchR <- addGroupCoverages(ArchRProj = mouse_multiome_ArchR, 
                                          groupBy = "Gex_macaque_cell_type", 
                                          minCells = 40, 
                                          maxCells = max(table(temp)), 
                                          maxFragments = 10^10, 
                                          minReplicates = 6, 
                                          maxReplicates = 6, 
                                          sampleRatio = 0.8)

temp <- as.data.frame(mouse_multiome_ArchR@projectMetadata$GroupCoverages$Gex_macaque_cell_type$coverageMetadata)

#call peak
pathToMacs2 <- '/data/User/sunym/env/MACS2/bin/macs2'

mouse_multiome_ArchR <- addReproduciblePeakSet(
  ArchRProj = mouse_multiome_ArchR, 
  groupBy = "Gex_macaque_cell_type", 
  pathToMacs2 = pathToMacs2,
  maxPeaks = 150000,
  cutOff = 0.05,
  force = TRUE
)

my_send_sms('call peak done!')

getPeakSet(mouse_multiome_ArchR)

#add peak matrix
mouse_multiome_ArchR <- addPeakMatrix(mouse_multiome_ArchR)
getAvailableMatrices(mouse_multiome_ArchR)

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003", load = FALSE, overwrite = TRUE)

# peak gene expression correlation ----------------------------------------
#load data
mouse_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003')

#add P2G
mouse_multiome_ArchR <- addPeak2GeneLinks(
  ArchRProj = mouse_multiome_ArchR,
  reducedDims = "IterativeLSI",
  useMatrix = 'GeneExpressionMatrix'
)

p2g <- getPeak2GeneLinks(
  ArchRProj = mouse_multiome_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "/data/User/sunym/project/Brain/processed_data/220802_summary/mouse_multiome_ArchR_221003", load = FALSE, overwrite = TRUE)

# save cell list ----------------------------------------------------------
cell_list <- rownames(mouse_multiome_ArchR@cellColData)
saveRDS(object = cell_list,file = './res/step_15_fig_221003/cell_list.rds')
