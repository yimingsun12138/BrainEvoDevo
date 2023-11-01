#####################################################################################
## Project: macaque fetal Brain multiome ATAC_seq                                  ##
## Script Purpose: preprocess of macaque multiome ATAC_seq                         ##
## Data: 2022.02.15                                                                ##
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
.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.3/',
            '/data/User/sunym/software/R_lib/R_4.1.3/'))
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
library(harmony,lib.loc = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.3/')
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(rtracklayer)
library(org.Mmu.eg.db)
library(clusterProfiler)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')
#initialize ArchR
addArchRThreads(threads = 5)

# genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmulatta.UCSC.rheMac10)
# 
# geneAnnotation <- rtracklayer::import('/data/User/sunym/project/Brain/data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf')
# geneAnnotation <- rtracklayer::as.data.frame(geneAnnotation)
# geneAnnotation$seqnames <- as.character(geneAnnotation$seqnames)
# geneAnnotation <- geneAnnotation[geneAnnotation$seqnames %in% c(as.character(1:20),'X','Y','MT'),]
# geneAnnotation$seqnames <- paste0('chr',geneAnnotation$seqnames)
# geneAnnotation <- GenomicRanges::makeGRangesFromDataFrame(geneAnnotation,keep.extra.columns = TRUE)
# 
# TxDb_object <- GenomicFeatures::makeTxDbFromGRanges(geneAnnotation)
# geneAnnotation <- createGeneAnnotation(TxDb = TxDb_object, OrgDb = org.Mmu.eg.db)

# create arrow file -------------------------------------------------------
inputFiles <- list.files(path = './fragment_file/macaque_multiome_data/')
temp <- base::lapply(X = inputFiles,FUN = function(x){
  return(strsplit(x = x,split = '_')[[1]][1])
})
temp <- unlist(temp)

inputFiles <- paste('./fragment_file/macaque_multiome_data',inputFiles,sep = '/')
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

# create archrproject -----------------------------------------------------
macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = './processed_data/macaque_multiome_ArchR_220216',
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

#add metadata
macaque_multiome_ArchR@cellColData$species <- 'macaque'
macaque_multiome_ArchR@cellColData$donor <- as.character(macaque_multiome_ArchR$Sample)
table(macaque_multiome_ArchR$donor)

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

pdf(file = './res/step_5_fig_220215/TSS_enrichment_vs_nFrag_dotplot_before_filter.pdf',width = 16,height = 5)
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

pdf(file = './res/step_5_fig_220215/TSS_enrichment_vlnplot_before_filter.pdf',width = 5,height = 5)
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

pdf(file = './res/step_5_fig_220215/nFrag_vlnplot_before_filter.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_5_fig_220215/frag_length_before_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_5_fig_220215/TSS_accessibility_by_sample_before_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#filter data
rownames(macaque_multiome_ArchR@cellColData)[1:3]
macaque_multiome_Seurat <- readRDS(file = '/data/User/sunym/project/Brain/res/step_15_fig_220208/macaque_multiome_meta_data.rds')
rownames(macaque_multiome_Seurat)[1:3]

cell_list <- rownames(macaque_multiome_Seurat)
cell_list <- base::lapply(X = cell_list,FUN = function(x){
  return(sub(pattern = '_',replacement = '#',x = x,fixed = TRUE))
})
cell_list <- unlist(cell_list)
table(cell_list %in% rownames(macaque_multiome_ArchR@cellColData))
macaque_multiome_ArchR <- macaque_multiome_ArchR[rownames(macaque_multiome_ArchR@cellColData) %in% cell_list]

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220216/", load = FALSE, overwrite = TRUE)

# process multiome ATAC_seq data ------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')

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

QC_A50A+QC_A50B+QC_A82A+QC_A82B+plot_layout(ncol = 4)

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

p

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

p

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

p + theme(legend.text = element_text(size = 10), legend.position = 'right')

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

p + theme(legend.text = element_text(size = 10), legend.position = 'right')

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
  dimsToUse = 1:30, 
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

#cluster C3 C9 C24 seems donor specific

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

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#marker gene
marker_gene <- c('PECAM1', #End
                 'PDGFRB', #Per
                 'CX3CR1', #Mic
                 'DNAH11','GJA1','SPATA17','LRRC71','SPAG17', #Ependymal
                 'SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3', #RG
                 'SOX10','OLIG2','EGFR', #OPC
                 'TOP2A','MKI67','CLSPN','AURKA', #Cyc
                 'EOMES','PPP1R17', #IP
                 'NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2', #Ex
                 'DLX5','GAD2','GAD1','DLX2', #In
                 'LHX6','SST', #InMGE
                 'SP8','NR2F2', #InCGE
                 'MEIS2','ETV1' #InPSB
)

macaque_multiome_ArchR <- addImputeWeights(macaque_multiome_ArchR,reducedDims = 'IterativeLSI')

p <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_multiome_ArchR)
)

#End
p1+p2+p$PECAM1+plot_layout(ncol = 3)
#cluster C1 C6

#Per
p1+p2+p$PDGFRB+plot_layout(ncol = 3)
#cluster C6

#Mic
p1+p2+p$CX3CR1+plot_layout(ncol = 3)
#cluster C1

#Ependymal
p1+p2+p$SPATA17+plot_layout(ncol = 3)
#not specific

#RG
p1+p2+p$VIM+plot_layout(ncol = 3)
#cluster C2 C3 C4 C5

#OPC
p1+p2+p$SOX10+plot_layout(ncol = 3)
p1+p2+p$OLIG2+plot_layout(ncol = 3)
#cluster C7 OPC
p1+p2+p$EGFR+plot_layout(ncol = 3)
#cluster C2 C3 RG-3

#Cyc
p1+p2+p$MKI67+plot_layout(ncol = 3)
#not specific

#IP
p1+p2+p$PPP1R17+plot_layout(ncol = 3)
#cluster C2 C12

#Ex
p1+p2+p$FEZF2+plot_layout(ncol = 3)
#cluster C11 C13 C14 C15 C16 C17

#In
p1+p2+p$GAD1+plot_layout(ncol = 3)
#cluster C19 - C27

#InMGE
p1+p2+p$LHX6+plot_layout(ncol = 3)
#cluster C19 C25 C23 C26 C24 C22 C28

#InCGE
p1+p2+p$SP8+plot_layout(ncol = 3)
#cluster C20 C21 C27

#InPSB
p1+p2+p$MEIS2+plot_layout(ncol = 3)
#not specific

#cluster filter
p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=TSSEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=nFrags)) + 
  geom_violin(mapping = aes(fill = Clusters))

p3 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=DoubletEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p1+p2+p3+plot_layout(ncol = 1)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#filter cluster C8 C9 C10 (double marker)
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(macaque_multiome_ArchR$Clusters %in% c('C8','C9','C10'))]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)
p1+p2+p3+plot_layout(ncol = 3)

# redo process ------------------------------------------------------------
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
  dimsToUse = 1:30, 
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

#cluster C17 C22 C10 C7 seems weird

#umap
macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:14, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

table(macaque_multiome_ArchR[macaque_multiome_ArchR$Clusters == 'C22']$donor)

#marker gene
marker_gene <- c('PECAM1', #End
                 'PDGFRB', #Per
                 'CX3CR1', #Mic
                 'DNAH11','GJA1','SPATA17','LRRC71','SPAG17', #Ependymal
                 'SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3', #RG
                 'SOX10','OLIG2','EGFR', #OPC
                 'TOP2A','MKI67','CLSPN','AURKA', #Cyc
                 'EOMES','PPP1R17', #IP
                 'NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2', #Ex
                 'DLX5','GAD2','GAD1','DLX2', #In
                 'LHX6','SST', #InMGE
                 'SP8','NR2F2', #InCGE
                 'MEIS2','ETV1' #InPSB
)

macaque_multiome_ArchR <- addImputeWeights(macaque_multiome_ArchR,reducedDims = 'IterativeLSI')

p <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_multiome_ArchR)
)

#In
p1+p2+p$GAD1+plot_layout(ncol = 3)

#InMGE
p1+p2+p$LHX6+plot_layout(ncol = 3)
p1+p2+p$SST+plot_layout(ncol = 3)

#InCGE
p1+p2+p$SP8+plot_layout(ncol = 3)
p1+p2+p$NR2F2+plot_layout(ncol = 3)

#filter cluster C19 C22 (doublet)
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(macaque_multiome_ArchR$Clusters %in% c('C19','C22'))]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)
p1+p2+p3+plot_layout(ncol = 3)

# redo process ------------------------------------------------------------
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
  dimsToUse = 1:30, 
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

#Cluster C5 seems weird

#umap
macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:14, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#cluster filter
p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=TSSEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=nFrags)) + 
  geom_violin(mapping = aes(fill = Clusters))

p3 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=DoubletEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p1+p2+p3+plot_layout(ncol = 1)

#filter cluster C9 (doublet)
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(macaque_multiome_ArchR$Clusters %in% c('C9'))]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#integrate with scRNA_seq
macaque_RNA_seurat <- readRDS(file = '/data/User/sunym/project/Brain/processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('Cyc-G2M','Cyc-S'))]

#recreate scRNA-seq data
temp <- macaque_RNA_seurat@assays$RNA@counts
meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_RNA','nFeature_RNA'))]
macaque_RNA_seurat <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
rm(list = c('temp','meta_data'))

#integrate
macaque_multiome_ArchR <- addImputeWeights(macaque_multiome_ArchR)
macaque_multiome_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = macaque_RNA_seurat,
  addToArrow = FALSE,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)

my_send_sms("intergration done!")

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")

p1+p2+plot_layout(ncol = 2)

#filter done! save data
macaque_multiome_ArchR@cellColData <- macaque_multiome_ArchR@cellColData[,!(colnames(macaque_multiome_ArchR@cellColData) %in% c('predictedCell_Un','predictedGroup_Un','predictedScore_Un'))]
macaque_multiome_ArchR@cellColData[1:3,]
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220216/", load = FALSE, overwrite = TRUE)

# final process -----------------------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')

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
  dimsToUse = 1:30, 
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

#cluster C22 C4 C16 seems weird
#umap
macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:14, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)
table(macaque_multiome_ArchR$Clusters)

#filter cluster C4
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(macaque_multiome_ArchR$Clusters %in% c('C4'))]

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220216/", load = FALSE, overwrite = TRUE)

# final process -----------------------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')

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
  dimsToUse = 1:12, 
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

#cluster C20 C7 seems weird

#umap
macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:12, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#low quality filter
p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=TSSEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=nFrags)) + 
  geom_violin(mapping = aes(fill = Clusters))

p3 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData),aes(x=Clusters,y=DoubletEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p1+p2+p3+plot_layout(ncol = 1)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220216/", load = FALSE, overwrite = TRUE)
meta_data <- macaque_multiome_ArchR@cellColData
meta_data <- as.data.frame(meta_data)
saveRDS(meta_data,file = './res/step_5_fig_220215/macaque_multiome_ArchR_meta_data.rds')
