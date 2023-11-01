#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: preprocess of mouse multiome ATAC data                          ##
## Data: 2022.07.11                                                                ##
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

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmusculus.UCSC.mm10)
# 
# geneAnnotation <- rtracklayer::import('/data/User/sunym/project/Brain/data/reference/ensembl_gtf_for_mapping/Mus_musculus.GRCm38.102.gtf')
# geneAnnotation <- rtracklayer::as.data.frame(geneAnnotation)
# geneAnnotation$seqnames <- as.character(geneAnnotation$seqnames)
# geneAnnotation <- geneAnnotation[geneAnnotation$seqnames %in% c(as.character(1:19),'X','Y','MT'),]
# geneAnnotation$seqnames <- paste0('chr',geneAnnotation$seqnames)
# geneAnnotation <- GenomicRanges::makeGRangesFromDataFrame(geneAnnotation,keep.extra.columns = TRUE)
# 
# TxDb_object <- GenomicFeatures::makeTxDbFromGRanges(geneAnnotation)
# geneAnnotation <- createGeneAnnotation(TxDb = TxDb_object, OrgDb = org.Mm.eg.db)

# create arrow file -------------------------------------------------------
inputFiles <- list.files(path = './fragment_file/mouse_multiome_data/220711/')
temp <- sub(pattern = '_fragments.tsv.gz',replacement = '',x = inputFiles,fixed = TRUE)
inputFiles <- paste('./fragment_file/mouse_multiome_data/220711',inputFiles,sep = '/')
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

# create archrproject from arrow file -----------------------------------------------------
ArrowFiles <- list.files(path = './arrow_file/mouse_multiome_data/220711/')
ArrowFiles <- paste('./arrow_file/mouse_multiome_data/220711',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

mouse_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = './processed_data/mouse_multiome_ArchR_220711',
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

#add meta data
mouse_multiome_ArchR$species <- 'mouse'
mouse_multiome_ArchR$donor <- mouse_multiome_ArchR$Sample

#filter data
meta_data <- readRDS(file = '../res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_meta_data.rds')
table(rownames(meta_data) %in% rownames(mouse_multiome_ArchR@cellColData))
cell_list <- dplyr::intersect(x = rownames(meta_data),y = rownames(mouse_multiome_ArchR@cellColData))
mouse_multiome_ArchR <- mouse_multiome_ArchR[cell_list]

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

pdf(file = './res/step_10_fig_220711/TSS_enrichment_vs_nFrag_dotplot.pdf',width = 12,height = 5)
QC_E145_1+QC_E155_1+QC_E155_2+plot_layout(ncol = 3)
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

pdf(file = './res/step_10_fig_220711/TSS_enrichment_vlnplot.pdf',width = 4,height = 4)
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

pdf(file = './res/step_10_fig_220711/nFrag_vlnplot_.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_10_fig_220711/frag_length_.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_10_fig_220711/TSS_accessibility_by_sample.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220711/", load = FALSE, overwrite = TRUE)

# dim reduction -----------------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220711/')
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
