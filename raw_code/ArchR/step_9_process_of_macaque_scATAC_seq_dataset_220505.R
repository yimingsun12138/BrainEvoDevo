#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: process of macaque scATAC_seq dataset                           ##
## Data: 2022.05.05                                                                ##
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
library(viridis)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load file ---------------------------------------------------------------
inputFiles <- list.files(path = './fragment_file/macaque_ATAC_data/')
temp <- base::lapply(X = inputFiles,FUN = function(x){
  x <- base::strsplit(x = x,split = '_')
  x <- x[[1]][1]
  return(x)
})
temp <- unlist(temp)
inputFiles <- paste('./fragment_file/macaque_ATAC_data',inputFiles,sep = '/')
names(inputFiles) <- temp
file.exists(inputFiles)

temp <- readRDS(file = './processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')

# create arrow file -------------------------------------------------------
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = getGeneAnnotation(temp),
  genomeAnnotation = getGenomeAnnotation(temp),
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

# create ArchR project ----------------------------------------------------
temp <- readRDS(file = './processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
ArrowFiles <- list.files(path = './arrow_file/macaque_ATAC_data/')
ArrowFiles <- paste('./arrow_file/macaque_ATAC_data',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

macaque_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = './processed_data/macaque_ATAC_ArchR_220506',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(temp),
  genomeAnnotation = getGenomeAnnotation(temp)
)

#add metadata
macaque_ATAC_ArchR$species <- 'macaque'
macaque_ATAC_ArchR$donor <- macaque_ATAC_ArchR$Sample
macaque_ATAC_ArchR$batch <- NA
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$donor %in% c('A68A','A68B','A84B','A84C'),"batch"] <- 'batch_1'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$donor %in% c('A50A','A82A','A82B'),"batch"] <- 'batch_2'
table(macaque_ATAC_ArchR$batch)

#filter data
temp <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')
table(rownames(temp@cellColData) %in% rownames(macaque_ATAC_ArchR@cellColData))
macaque_ATAC_ArchR <- macaque_ATAC_ArchR[rownames(macaque_ATAC_ArchR@cellColData) %in% rownames(temp@cellColData)]
macaque_ATAC_ArchR$cell_type <- NA
macaque_ATAC_ArchR@cellColData[rownames(temp@cellColData),"cell_type"] <- temp$cell_type
table(macaque_ATAC_ArchR$cell_type)

macaque_ATAC_ArchR <- macaque_ATAC_ArchR[!(macaque_ATAC_ArchR$cell_type %in% c('RG-3'))]

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)

# QC plot -----------------------------------------------------------------
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')
#QC plot
df <- macaque_ATAC_ArchR@cellColData[,c('nFrags','TSSEnrichment','Sample')]
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

pdf(file = './res/step_9_fig_220505/TSS_enrichment_vs_nFrag_dotplot.pdf',width = 16,height = 10)
QC_A50A+QC_A68A+QC_A68B+QC_A82A+QC_A82B+QC_A84B+QC_A84C+plot_layout(ncol = 4)
dev.off()

#tss enrichment
p <- plotGroups(ArchRProj = macaque_ATAC_ArchR, 
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

pdf(file = './res/step_9_fig_220505/TSS_enrichment_vlnplot.pdf',width = 5,height = 5)
p
dev.off()

#nFrag distribution
p <- plotGroups(ArchRProj = macaque_ATAC_ArchR, 
                groupBy = "Sample", 
                colorBy = "cellColData", 
                name = "log10(nFrags)",
                plotAs = "violin", 
                alpha = 0.4, 
                addBoxPlot = TRUE) + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12,face = 'bold')) + xlab('')

pdf(file = './res/step_9_fig_220505/nFrag_vlnplot.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_9_fig_220505/frag_length_by_sample.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_9_fig_220505/frag_length_by_batch.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_9_fig_220505/TSS_accessibility_by_sample.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_9_fig_220505/TSS_accessibility_by_batch.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)

# LSI dim reduction -------------------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')
#LSI
macaque_ATAC_ArchR <- addIterativeLSI(ArchRProj = macaque_ATAC_ArchR,
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
macaque_ATAC_ArchR <- addClusters(
  input = macaque_ATAC_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:16, 
  force = TRUE
)

#umap
macaque_ATAC_ArchR <- addUMAP(
  ArchRProj = macaque_ATAC_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:16, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "batch", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "donor", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

#the batch is so strong!


# harmony dim reduction ----------------------------------------------------
macaque_ATAC_ArchR <- addHarmony(ArchRProj = macaque_ATAC_ArchR,
                                 reducedDims = 'IterativeLSI',
                                 name = 'harmonyLSI',
                                 groupBy = 'batch',
                                 verbose = TRUE,
                                 force = TRUE)

#cluster
macaque_ATAC_ArchR <- addClusters(
  input = macaque_ATAC_ArchR,
  reducedDims = "harmonyLSI",
  method = "Seurat",
  name = "harmonyClusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:16, 
  force = TRUE
)

#umap
macaque_ATAC_ArchR <- addUMAP(
  ArchRProj = macaque_ATAC_ArchR, 
  reducedDims = "harmonyLSI", 
  name = "harmonyUMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:16, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "batch", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p1+p2+p3+plot_layout(ncol = 3)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)

# annotate based on clustering --------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')
macaque_ATAC_ArchR <- addImputeWeights(ArchRProj = macaque_ATAC_ArchR,reducedDims = 'IterativeLSI')

#marker express
marker_gene <- c('PECAM1', #End
                 'PDGFRB', #Per
                 'CX3CR1', #Mic
                 'DNAH11','GJA1','SPATA17','LRRC71','SPAG17', #Ependymal
                 'SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3', #RG
                 'SOX10','OLIG2','EGFR', #OPC
                 'TOP2A','MKI67','CLSPN','AURKA', #Cyc
                 'EOMES','PPP1R17', #IP
                 'NEUROD2','NEUROD6','TBR1','SATB2','FEZF2', #Ex
                 'CUX2','HS6ST3', #Ex-1
                 'NWD2','SNTG1', #Ex-2
                 'GRIK3','TRPM3', #Ex-3
                 'DLX5','GAD2','GAD1','DLX2', #In
                 'LHX6','SST', #InMGE
                 'SP8','NR2F2', #InCGE
                 'MEIS2','ETV1' #InPSB
)

p <- plotEmbedding(
  ArchRProj = macaque_ATAC_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "harmonyUMAP",
  imputeWeights = getImputeWeights(macaque_ATAC_ArchR)
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")

#End cluster C2
p1+p2+p$PECAM1+plot_layout(ncol = 3)

#Per cluster C2
p1+p2+p$PDGFRB+plot_layout(ncol = 3)

#Mic cluster C1
p1+p2+p$CX3CR1+plot_layout(ncol = 3)

#Ependymal not specific
p1+p2+p$DNAH11+plot_layout(ncol = 3)

#RG cluster C20 C21
p1+p2+p$SOX9+plot_layout(ncol = 3)
p1+p2+p$PAX6+plot_layout(ncol = 3)
p1+p2+p$VIM+plot_layout(ncol = 3)

#OPC cluster C18
p1+p2+p$SOX10+plot_layout(ncol = 3)
p1+p2+p$OLIG2+plot_layout(ncol = 3)
p1+p2+p$EGFR+plot_layout(ncol = 3)

#Cycling not specific
p1+p2+p$TOP2A+plot_layout(ncol = 3)
p1+p2+p$MKI67+plot_layout(ncol = 3)

#IP cluster C19
p1+p2+p$EOMES+plot_layout(ncol = 3)
p1+p2+p$PPP1R17+plot_layout(ncol = 3)

#Ex cluster C11 C12 C13 C14 C15 C16 C17
p1+p2+p$NEUROD2+plot_layout(ncol = 3)
p1+p2+p$NEUROD6+plot_layout(ncol = 3)

#Ex-1 cluster C12 C13 C15
p1+p2+p$NEUROD2+plot_layout(ncol = 3)
p1+p2+p$CUX2+plot_layout(ncol = 3)
p1+p2+p$HS6ST3+plot_layout(ncol = 3)

#Ex-2 cluster C14 C16 C17
p1+p2+p$NWD2+plot_layout(ncol = 3)
p1+p2+p$SNTG1+plot_layout(ncol = 3)
p1+p2+p$SATB2+plot_layout(ncol = 3)

#Ex-3 cluster C11
p1+p2+p$GRIK3+plot_layout(ncol = 3)
p1+p2+p$TRPM3+plot_layout(ncol = 3)

#InMGE cluster C6 C7 C8 C9 C10
p1+p2+p$LHX6+plot_layout(ncol = 3)
p1+p2+p$SST+plot_layout(ncol = 3)

#InCGE cluster C3 C4 C5
p1+p2+p$SP8+plot_layout(ncol = 3)
p1+p2+p$NR2F2+plot_layout(ncol = 3)

#annotate
p1+p2+plot_layout(ncol = 2)
macaque_ATAC_ArchR$cell_type <- NA
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C2'),"cell_type"] <- 'End/Per'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C1'),"cell_type"] <- 'Mic'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C18'),"cell_type"] <- 'OPC'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C20','C21'),"cell_type"] <- 'RG'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C19'),"cell_type"] <- 'IP'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C12','C13','C15'),"cell_type"] <- 'Ex-1'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C14','C16','C17'),"cell_type"] <- 'Ex-2'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C11'),"cell_type"] <- 'Ex-3'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C3','C4','C5'),"cell_type"] <- 'InCGE'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C6','C7','C8','C9','C10'),"cell_type"] <- 'InMGE'

#plot
p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p1+p2+plot_layout(ncol = 2)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)

# label transfer ----------------------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')
macaque_RNA_Seurat <- readRDS(file = '/data/User/sunym/project/Brain/processed_data/220504_summary/macaque_RNA_Seurat_220504.rds')

#check data
p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p1+p2+plot_layout(ncol = 2)

p1 <- DimPlot(macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_RNA_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#recreate macaque_RNA_Seurat
temp <- macaque_RNA_Seurat@meta.data
temp <- temp[,!(colnames(temp) %in% c('orig.ident','nCount_RNA','nFeature_RNA'))]
macaque_RNA_Seurat <- macaque_RNA_Seurat@assays$RNA@counts
macaque_RNA_Seurat <- CreateSeuratObject(counts = macaque_RNA_Seurat,project = 'temp',assay = 'RNA',meta.data = temp,min.cells = 0,min.features = 0)

#do constrained intehration
group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('End','Mic','OPC','Per')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C1','C2','C18')]
  ),
  Ex_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','IP','RG-1','RG-2','RG-3')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C20','C21','C19','C11','C12','C13','C14','C15','C16','C17')]
  ),
  In_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C3','C4','C5','C6','C7','C8','C9','C10')]
  )
)

macaque_ATAC_ArchR <- addImputeWeights(ArchRProj = macaque_ATAC_ArchR,reducedDims = 'IterativeLSI')
macaque_ATAC_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_ATAC_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = macaque_RNA_Seurat,
  addToArrow = TRUE, 
  groupList = group_list,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

my_send_sms('Seurat integration done!')

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "harmonyUMAP")
p1+p2+plot_layout(ncol = 2)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)

#compare cell type annotation
scibet::Confusion_heatmap(ori = macaque_ATAC_ArchR$cell_type,prd = macaque_ATAC_ArchR$predictedGroup)

# call peak ---------------------------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')
p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "harmonyUMAP")
p1+p2+plot_layout(ncol = 2)

#call peak
macaque_ATAC_ArchR <- addGroupCoverages(ArchRProj = macaque_ATAC_ArchR, 
                                        groupBy = "cell_type", 
                                        minCells = 40, 
                                        maxCells = 500, 
                                        minReplicates = 2, 
                                        maxReplicates = 5, 
                                        sampleRatio = 0.8)


pathToMacs2 <- '/data/User/sunym/env/MACS2/bin/macs2'

macaque_ATAC_ArchR <- addReproduciblePeakSet(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  pathToMacs2 = pathToMacs2,
  maxPeaks = 200000,
  genomeSize = 2.7e+09,
  force = TRUE
)

my_send_sms('call peak done!')

getPeakSet(ArchRProj = macaque_ATAC_ArchR)

#add peakmatrix
macaque_ATAC_ArchR <- addPeakMatrix(ArchRProj = macaque_ATAC_ArchR)
getAvailableMatrices(macaque_ATAC_ArchR)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)

# peak co-accessibility ---------------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')

#calculate co-accessibility
macaque_ATAC_ArchR <- addCoAccessibility(
  ArchRProj = macaque_ATAC_ArchR,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = macaque_ATAC_ArchR,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA

#plot
markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(macaque_ATAC_ArchR,resolution = 2000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG','OPC','Mic','End/Per')
)

#RG
grid::grid.newpage()
grid::grid.draw(p$SOX9)
grid::grid.newpage()
grid::grid.draw(p$PAX6)
grid::grid.newpage()
grid::grid.draw(p$VIM)
grid::grid.newpage()
grid::grid.draw(p$HOPX)

#IP
grid::grid.newpage()
grid::grid.draw(p$EOMES)
grid::grid.newpage()
grid::grid.draw(p$PPP1R17)

#OPC
grid::grid.newpage()
grid::grid.draw(p$OLIG2)
grid::grid.newpage()
grid::grid.draw(p$EGFR)

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

# peak gene express correlation -------------------------------------------
macaque_ATAC_ArchR <- addPeak2GeneLinks(
  ArchRProj = macaque_ATAC_ArchR,
  reducedDims = "IterativeLSI",
  useMatrix = 'GeneIntegrationMatrix'
)

p2g <- getPeak2GeneLinks(
  ArchRProj = macaque_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g

#plot
markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_ATAC_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG','OPC','Mic','End/Per')
)

#RG
grid::grid.newpage()
grid::grid.draw(p$SOX9)
grid::grid.newpage()
grid::grid.draw(p$PAX6)
grid::grid.newpage()
grid::grid.draw(p$VIM)
grid::grid.newpage()
grid::grid.draw(p$HOPX)

#IP
grid::grid.newpage()
grid::grid.draw(p$EOMES)
grid::grid.newpage()
grid::grid.draw(p$PPP1R17)

#OPC
grid::grid.newpage()
grid::grid.draw(p$OLIG2)
grid::grid.newpage()
grid::grid.draw(p$EGFR)

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
grid::grid.draw(p$GAD2)
grid::grid.newpage()
grid::grid.draw(p$DLX5)

p <- plotPeak2GeneHeatmap(ArchRProj = macaque_ATAC_ArchR, groupBy = "cell_type", k = 10)
ComplexHeatmap::draw(p, heatmap_legend_side = "bot", annotation_legend_side = "bot", width = unit(2,'inches'), heigh = unit(8,'inches'))

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)

# chromVar motif activity analysis ----------------------------------------

#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')

#add motif annotation
macaque_ATAC_ArchR <- addMotifAnnotations(ArchRProj = macaque_ATAC_ArchR,motifSet = 'cisbp',name = 'cisbpMotif',species = 'homo sapiens',force = TRUE)

#add background peaks
macaque_ATAC_ArchR <- addBgdPeaks(macaque_ATAC_ArchR)

#chromvar
macaque_ATAC_ArchR <- addDeviationsMatrix(ArchRProj = macaque_ATAC_ArchR,
                                          peakAnnotation = 'cisbpMotif',
                                          force = TRUE)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220506/", load = FALSE, overwrite = TRUE)
