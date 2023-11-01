#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: preprocess of macaque brain scATAC_seq                          ##
## Data: 2022.01.09                                                                ##
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
library(harmony)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(rtracklayer)
library(org.Mmu.eg.db)
library(clusterProfiler)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

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
inputFiles <- list.files(path = './fragment_file/')
temp <- base::lapply(X = inputFiles,FUN = function(x){
  return(strsplit(x = x,split = '_')[[1]][1])
})
temp <- unlist(temp)

inputFiles <- paste('./fragment_file',inputFiles,sep = '/')
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

# doublet inference -------------------------------------------------------
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# create archrproject -----------------------------------------------------
macaque_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = './processed_data/macaque_ATAC_ArchR_220109',
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

#add metadata
macaque_ATAC_ArchR@cellColData$species <- 'macaque'
macaque_ATAC_ArchR@cellColData$donor <- as.character(macaque_ATAC_ArchR$Sample)
macaque_ATAC_ArchR@cellColData$batch <- NA
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$donor %in% c('A68A','A68B','A84B','A84C'),"batch"] <- 'batch_1'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$donor %in% c('A50A','A82A','A82B'),"batch"] <- 'batch_2'
table(macaque_ATAC_ArchR$batch)

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

pdf(file = './res/step_1_fig_220109/TSS_enrichment_vs_nFrag_dotplot_before_doublet_filter.pdf',width = 16,height = 10)
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

pdf(file = './res/step_1_fig_220109/TSS_enrichment_vlnplot_before_doublet_filter.pdf',width = 5,height = 5)
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

pdf(file = './res/step_1_fig_220109/nFrag_vlnplot_before_doublet_filter.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_1_fig_220109/frag_length_by_sample_before_doublet_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_1_fig_220109/frag_length_by_batch_before_doublet_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_1_fig_220109/TSS_accessibility_by_sample_before_doublet_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_1_fig_220109/TSS_accessibility_by_batch_before_doublet_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220109/", load = FALSE, overwrite = TRUE)

# batch 2 -----------------------------------------------------------------
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')
batch_2_ArchR <- macaque_ATAC_ArchR[macaque_ATAC_ArchR$batch == 'batch_2',]
table(batch_2_ArchR$batch)

#LSI
batch_2_ArchR <- addIterativeLSI(ArchRProj = batch_2_ArchR,
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
batch_2_ArchR <- addClusters(
  input = batch_2_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

cM <- confusionMatrix(i = batch_2_ArchR$Clusters,j = batch_2_ArchR$Sample)
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
#umap
batch_2_ArchR <- addUMAP(
  ArchRProj = batch_2_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:15, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = batch_2_ArchR, 
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

batch_2_ArchR <- addImputeWeights(batch_2_ArchR,reducedDims = 'IterativeLSI')

p <- plotEmbedding(
  ArchRProj = batch_2_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(batch_2_ArchR)
)

#End
p1+p2+p$PECAM1+plot_layout(ncol = 3)
#cluster 12 1

#Per
p1+p2+p$PDGFRB+plot_layout(ncol = 3)
#cluster 1

#Mic
p1+p2+p$CX3CR1+plot_layout(ncol = 3)
#cluster 1

#Ependymal
p1+p2+p$SPAG17+plot_layout(ncol = 3)
#not specific

#RG
p1+p2+p$MOXD1+plot_layout(ncol = 3)
#cluster16 17

#OPC
p1+p2+p$EGFR+plot_layout(ncol = 3)
#cluster 15

#Cyc
p1+p2+p$AURKA+plot_layout(ncol = 3)
#not specific

#IP
p1+p2+p$PPP1R17+plot_layout(ncol = 3)
#cluster 11 10 18

#Ex
p1+p2+p$FEZF2+plot_layout(ncol = 3)
#cluster 12 13 14

#In
p1+p2+p$GAD2+plot_layout(ncol = 3)
#cluster 3 4 5 6 19 20 21 22

#InMGE
p1+p2+p$SST+plot_layout(ncol = 3)
#cluster 4 5 21 22

#InCGE
p1+p2+p$NR2F2+plot_layout(ncol = 3)
#cluster 3 6 19 20

#InPSB
p1+p2+p$ETV1+plot_layout(ncol = 3)
#cluster 2 3


#filter clusters
#cluster not included:7 8 9

cM <- My_Cell_Proportion(meta_data = batch_2_ArchR@cellColData,group.by = 'Sample',split.by = 'Clusters')
ggplot(data = cM, mapping = aes(x = Clusters,y = Proportion,fill = Sample))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of Sample contribution to each cluster',fill = 'Sample')+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()
#filter cluster 8, donor specific

p1 <- ggplot(data = as.data.frame(batch_2_ArchR@cellColData),aes(x=Clusters,y=TSSEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p2 <- ggplot(data = as.data.frame(batch_2_ArchR@cellColData),aes(x=Clusters,y=nFrags)) + 
  geom_violin(mapping = aes(fill = Clusters))

p3 <- ggplot(data = as.data.frame(batch_2_ArchR@cellColData),aes(x=Clusters,y=DoubletEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p1+p2+p3+plot_layout(ncol = 1)

#filter cluster 7 8 9
#cluster remain valid:2 3 16
batch_2_ArchR <- batch_2_ArchR[!(batch_2_ArchR$Clusters %in% c('C7','C8','C9')),]

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = batch_2_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#re-do cluster and umap
#LSI
batch_2_ArchR <- addIterativeLSI(ArchRProj = batch_2_ArchR,
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
batch_2_ArchR <- addClusters(
  input = batch_2_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "re_Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

#umap
batch_2_ArchR <- addUMAP(
  ArchRProj = batch_2_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:15, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

p1+p2+plot_layout(ncol = 2)

#integrate with scRNA-seq
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
batch_2_ArchR <- addImputeWeights(batch_2_ArchR)
batch_2_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = batch_2_ArchR, 
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

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")

p1+p2+plot_layout(ncol = 2)

#cluster 2 3 InPSB filter them
batch_2_ArchR <- batch_2_ArchR[!(batch_2_ArchR$Clusters %in% c('C2','C3')),]

#re-do cluster and umap
#LSI
batch_2_ArchR <- addIterativeLSI(ArchRProj = batch_2_ArchR,
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
batch_2_ArchR <- addClusters(
  input = batch_2_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "re_Clusters",
  resolution = 0.8, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

#umap
batch_2_ArchR <- addUMAP(
  ArchRProj = batch_2_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:20, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "re_Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = batch_2_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)
p1+p2+p3+plot_layout(ncol = 3)

#quality control
#cell type specific
cM <- My_Cell_Proportion(meta_data = batch_2_ArchR@cellColData,group.by = 'Sample',split.by = 're_Clusters')
ggplot(data = cM, mapping = aes(x = re_Clusters,y = Proportion,fill = Sample))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of Sample contribution to each cluster',fill = 'Sample')+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()
#depth
p1 <- ggplot(data = as.data.frame(batch_2_ArchR@cellColData),aes(x=re_Clusters,y=TSSEnrichment)) + 
  geom_violin(mapping = aes(fill = re_Clusters))

p2 <- ggplot(data = as.data.frame(batch_2_ArchR@cellColData),aes(x=re_Clusters,y=nFrags)) + 
  geom_violin(mapping = aes(fill = re_Clusters))

p3 <- ggplot(data = as.data.frame(batch_2_ArchR@cellColData),aes(x=re_Clusters,y=DoubletEnrichment)) + 
  geom_violin(mapping = aes(fill = re_Clusters))

p1+p2+p3+plot_layout(ncol = 1)

#marker express
batch_2_ArchR <- addImputeWeights(ArchRProj = batch_2_ArchR,reducedDims = 'IterativeLSI')

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

p <- plotEmbedding(
  ArchRProj = batch_2_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(batch_2_ArchR)
)

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "re_Clusters", embedding = "UMAP")

p1+p$NR2F2
#no cluster18

#cluster 18 marker
marker_gene <- getMarkerFeatures(
  ArchRProj = batch_2_ArchR, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "re_Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C18"
)

marker_gene <- getMarkers(marker_gene, cutOff = "FDR <= 0.05 & Log2FC >= 1")
p <- plotEmbedding(
  ArchRProj = batch_2_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene$C18$name[1:12], 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(batch_2_ArchR)
)

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "re_Clusters", embedding = "UMAP")
p1+p$CHST3
#no specific marker

#label transfer
macaque_RNA_seurat <- readRDS(file = '/data/User/sunym/project/Brain/processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('Cyc-G2M','Cyc-S','InPSB'))]

temp <- macaque_RNA_seurat@assays$RNA@counts
meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_RNA','nFeature_RNA'))]
macaque_RNA_seurat <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
rm(list = c('temp','meta_data'))

batch_2_ArchR <- addImputeWeights(batch_2_ArchR)
batch_2_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = batch_2_ArchR, 
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

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "re_Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)
table(batch_2_ArchR@cellColData[batch_2_ArchR$re_Clusters == 'C18',"predictedGroup_Un"])

#cluster 18 maybe doublet
batch_2_ArchR <- batch_2_ArchR[!(batch_2_ArchR$re_Clusters %in% c('C18')),]

#LSI
batch_2_ArchR <- addIterativeLSI(ArchRProj = batch_2_ArchR,
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
batch_2_ArchR <- addClusters(
  input = batch_2_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "re_Clusters",
  resolution = 0.8, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

#umap
batch_2_ArchR <- addUMAP(
  ArchRProj = batch_2_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:20, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "re_Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_2_ArchR, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)
meta_data <- batch_2_ArchR@cellColData
meta_data <- data.frame(meta_data)
meta_data$Clusters <- meta_data$re_Clusters
meta_data <- meta_data[,!(colnames(meta_data) %in% c('re_Clusters'))]
saveRDS(meta_data,file = './res/step_1_fig_220109/batch_2_meta_data_after_filter.rds')

# batch 1 -----------------------------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')
batch_1_ArchR <- macaque_ATAC_ArchR[macaque_ATAC_ArchR$batch == 'batch_1',]
table(batch_1_ArchR$batch)

#LSI
batch_1_ArchR <- addIterativeLSI(ArchRProj = batch_1_ArchR,
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
batch_1_ArchR <- addClusters(
  input = batch_1_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

cM <- My_Cell_Proportion(meta_data = batch_1_ArchR@cellColData,group.by = 'Sample','Clusters')
ggplot(data = cM, mapping = aes(x = Clusters,y = Proportion,fill = Sample))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of Sample contribution to each cluster',fill = 'Sample')+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()

#C12 C15 C4 C6 may be donor specific

#umap
batch_1_ArchR <- addUMAP(
  ArchRProj = batch_1_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:15, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = batch_1_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#marker express
batch_1_ArchR <- addImputeWeights(ArchRProj = batch_1_ArchR,reducedDims = 'IterativeLSI')

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

p <- plotEmbedding(
  ArchRProj = batch_1_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(batch_1_ArchR)
)

p1 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#End
p1+p2+p$PECAM1+plot_layout(ncol = 3)
#cluster 1 2 18 14

#Per
p1+p2+p$PDGFRB+plot_layout(ncol = 3)
#not specific

#Mic
p1+p2+p$CX3CR1+plot_layout(ncol = 3)
#cluster 18 2

#Ependymal
p1+p2+p$SPATA17+plot_layout(ncol = 3)
#not specific

#RG
p1+p2+p$HOPX+plot_layout(ncol = 3)
#cluster 10 11 12 13

#OPC
p1+p2+p$OLIG2+plot_layout(ncol = 3)
#cluster 12 13

#Cyc
p1+p2+p$MKI67+plot_layout(ncol = 3)
#not specific

#IP
p1+p2+p$PPP1R17+plot_layout(ncol = 3)
#cluster 10

#Ex
p1+p2+p$SATB2+plot_layout(ncol = 3)
#cluster 9 17 15 14 16

#In
p1+p2+p$GAD1+plot_layout(ncol = 3)
#cluster 3 4 5 6 7 8

#InMGE
p1+p2+p$LHX6+plot_layout(ncol = 3)
#cluster 7 8 

#InCGE
p1+p2+p$SP8+plot_layout(ncol = 3)
#cluster 3 4 5

#InPSB
p1+p2+p$MEIS2+plot_layout(ncol = 3)
#cluster 6

#cluster remain to be determined:cluster 16 18 19

#quality control
p1 <- ggplot(data = as.data.frame(batch_1_ArchR@cellColData),aes(x=Clusters,y=TSSEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p2 <- ggplot(data = as.data.frame(batch_1_ArchR@cellColData),aes(x=Clusters,y=nFrags)) + 
  geom_violin(mapping = aes(fill = Clusters))

p3 <- ggplot(data = as.data.frame(batch_1_ArchR@cellColData),aes(x=Clusters,y=DoubletEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p1+p2+p3+plot_layout(ncol = 1)

#cluster 18 low quality
batch_1_ArchR <- batch_1_ArchR[!(batch_1_ArchR$Clusters %in% c('C18')),]

#integrate with scRNA-seq
macaque_RNA_seurat <- readRDS(file = '/data/User/sunym/project/Brain/processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('Cyc-G2M','Cyc-S'))]

temp <- macaque_RNA_seurat@assays$RNA@counts
meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_RNA','nFeature_RNA'))]
macaque_RNA_seurat <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
rm(list = c('temp','meta_data'))

batch_1_ArchR <- addImputeWeights(batch_1_ArchR)
batch_1_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = batch_1_ArchR, 
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

p1 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

table(batch_1_ArchR@cellColData[batch_1_ArchR$Clusters == 'C19',"predictedGroup_Un"])
table(batch_1_ArchR@cellColData[batch_1_ArchR$Clusters == 'C16',"predictedGroup_Un"])

#cluster 6 InPSB
#cluster 16 19 doublet
batch_1_ArchR <- batch_1_ArchR[!(batch_1_ArchR$Clusters %in% c('C19','C6','C16')),]

#re-do
#LSI
batch_1_ArchR <- addIterativeLSI(ArchRProj = batch_1_ArchR,
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
batch_1_ArchR <- addClusters(
  input = batch_1_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

cM <- My_Cell_Proportion(meta_data = batch_1_ArchR@cellColData,group.by = 'Sample','Clusters')
ggplot(data = cM, mapping = aes(x = Clusters,y = Proportion,fill = Sample))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of Sample contribution to each cluster',fill = 'Sample')+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()

#cluster C14 C19 may be donor specific

#umap
batch_1_ArchR <- addUMAP(
  ArchRProj = batch_1_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:30, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = batch_1_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

#quality control
p1 <- ggplot(data = as.data.frame(batch_1_ArchR@cellColData),aes(x=Clusters,y=TSSEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p2 <- ggplot(data = as.data.frame(batch_1_ArchR@cellColData),aes(x=Clusters,y=nFrags)) + 
  geom_violin(mapping = aes(fill = Clusters))

p3 <- ggplot(data = as.data.frame(batch_1_ArchR@cellColData),aes(x=Clusters,y=DoubletEnrichment)) + 
  geom_violin(mapping = aes(fill = Clusters))

p1+p2+p3+plot_layout(ncol = 1)

#marker express
batch_1_ArchR <- addImputeWeights(ArchRProj = batch_1_ArchR,reducedDims = 'IterativeLSI')

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

p <- plotEmbedding(
  ArchRProj = batch_1_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(batch_1_ArchR)
)

p1 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#End
p1+p2+p$PECAM1+plot_layout(ncol = 3)
#cluster C1

#Per
p1+p2+p$PDGFRB+plot_layout(ncol = 3)
#not specific

#Mic
p1+p2+p$CX3CR1+plot_layout(ncol = 3)
#cluster C2

#Ependymal
p1+p2+p$SPAG17+plot_layout(ncol = 3)
#not specific

#RG
p1+p2+p$SOX9+plot_layout(ncol = 3)
#cluster C9 C10 C12

#OPC
p1+p2+p$EGFR+plot_layout(ncol = 3)
#cluster C9 RG-3 cluster C8 OPC

#Cyc
p1+p2+p$MKI67+plot_layout(ncol = 3)
#not specific

#IP
p1+p2+p$PPP1R17+plot_layout(ncol = 3)
#cluster C12

#Ex
p1+p2+p$FEZF2+plot_layout(ncol = 3)
#cluster C3 4 5 6 11

#In
p1+p2+p$DLX2+plot_layout(ncol = 3)
#cluster C7 13 14 15 16 17

#InCGE
p1+p2+p$NR2F2+plot_layout(ncol = 3)
#cluster C13 14 15

#InMGE
p1+p2+p$LHX6+plot_layout(ncol = 3)
#cluster 7 16 17

#integrate with scRNA-seq
macaque_RNA_seurat <- readRDS(file = '/data/User/sunym/project/Brain/processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('Cyc-G2M','Cyc-S'))]

temp <- macaque_RNA_seurat@assays$RNA@counts
meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_RNA','nFeature_RNA'))]
macaque_RNA_seurat <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
rm(list = c('temp','meta_data'))
gc()

batch_1_ArchR <- addImputeWeights(batch_1_ArchR)
batch_1_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = batch_1_ArchR, 
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

p1 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = batch_1_ArchR, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

table(batch_1_ArchR@cellColData[batch_1_ArchR$Clusters == 'C7',"predictedGroup_Un"])
meta_data <- batch_1_ArchR@cellColData
meta_data <- as.data.frame(meta_data)
saveRDS(meta_data,file = './res/step_1_fig_220109/batch_1_meta_data_after_filter.rds')
