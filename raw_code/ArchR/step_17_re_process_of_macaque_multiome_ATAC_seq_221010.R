#####################################################################################
## Project: macaque fetal Brain multiome ATAC_seq                                  ##
## Script Purpose: re_process of macaque multiome ATAC_seq                         ##
## Data: 2022.10.10                                                                ##
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
#load macaque multiome ArchR old data
temp <- loadArchRProject(path = '../processed_data/220718_summary/macaque_multiome_ArchR_220720/')

#load macaque multiome Seurat
macaque_multiome_Seurat <- readRDS(file = '../processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
table(rownames(temp@cellColData) %in% colnames(macaque_multiome_Seurat))
#come cells are filtered!

#load Arrow file
ArrowFiles <- list.files(path = './arrow_file/macaque_multiome_data/')
ArrowFiles <- paste('./arrow_file/macaque_multiome_data',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

#create ArchR project
macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

#add metadata
macaque_multiome_ArchR@cellColData$species <- 'macaque'
macaque_multiome_ArchR@cellColData$donor <- as.character(macaque_multiome_ArchR$Sample)
table(macaque_multiome_ArchR$donor)

#filter cells
table(colnames(macaque_multiome_Seurat) %in% rownames(macaque_multiome_ArchR@cellColData))
macaque_multiome_ArchR <- macaque_multiome_ArchR[colnames(macaque_multiome_Seurat)]

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/', load = FALSE, overwrite = TRUE)

# QC plot -----------------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')

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

pdf(file = './res/step_17_fig_221010/TSS_enrichment_vs_nFrag_dotplot_after_filter.pdf',width = 16,height = 5)
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

pdf(file = './res/step_17_fig_221010/TSS_enrichment_vlnplot_after_filter.pdf',width = 5,height = 5)
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

pdf(file = './res/step_17_fig_221010/nFrag_vlnplot_after_filter.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_17_fig_221010/frag_length_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_17_fig_221010/TSS_accessibility_by_sample_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/', load = FALSE, overwrite = TRUE)

# add gene expression matrix ----------------------------------------------
#load processed data
macaque_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')
macaque_multiome_Seurat <- readRDS(file = '../processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

#gene expression matrix
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
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/', load = FALSE, overwrite = TRUE)

# dim reduction -----------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')

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

#ndim 15 is ok!
#cluster
macaque_multiome_ArchR <- addClusters(
  input = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:14, 
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
  dimsToUse = 1:14, 
  force = TRUE
)

#add cell type
macaque_multiome_Seurat <- readRDS(file = '../processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
macaque_multiome_ArchR$Gex_cell_type <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]
macaque_multiome_ArchR$Gex_sub_cell_type <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub_cell_type"]
macaque_multiome_ArchR$cell_type <- macaque_multiome_ArchR$Gex_cell_type

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
pdf(file = './res/step_17_fig_221010/macaque_multiome_ArchR_cell_type_dimplot.pdf',width = 18,height = 8)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_sub_cell_type", embedding = "UMAP")
pdf(file = './res/step_17_fig_221010/macaque_multiome_ArchR_sub_cell_type_dimplot.pdf',width = 18,height = 8)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

scibet::Confusion_heatmap(ori = macaque_multiome_ArchR$Clusters,prd = macaque_multiome_ArchR$cell_type)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/', load = FALSE, overwrite = TRUE)

# label transfer ----------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')
macaque_multiome_Seurat <- readRDS(file = '../processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#constrained integration
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('End','Mic','OPC','Per','VLMC')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C1','C2','C10')]
  ),
  npc = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('Cycling','IP','RG-1','RG-2')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C11','C12','C14')]
  ),
  Ex = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C3','C4','C5','C6','C7','C8','C9','C13')]
  ),
  In = SimpleList(
    RNA = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C15','C16','C17','C18','C19','C20','C21')]
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

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/', load = FALSE, overwrite = TRUE)

# call peak ---------------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+p4+plot_layout(ncol = 4)

#notice point:
#maxpeaks = 150k
#macs2 FDR = 0.05
#modify mincells and max cells while doing pseudobulk
#mincells: 70

#add pseudo bulk
temp <- paste(macaque_multiome_ArchR$donor,macaque_multiome_ArchR$cell_type,sep = '_')
table(temp) %>% max()
table(temp) %>% min()
macaque_multiome_ArchR <- addGroupCoverages(ArchRProj = macaque_multiome_ArchR, 
                                            groupBy = "cell_type", 
                                            minCells = 70, 
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

# #sleep
# ii <- 1
# while(1){
#   cat(paste("round",ii),sep = "\n")
#   ii <- ii+1
#   Sys.sleep(30)
# }

getPeakSet(macaque_multiome_ArchR)

#add peak matrix
macaque_multiome_ArchR <- addPeakMatrix(macaque_multiome_ArchR)
getAvailableMatrices(macaque_multiome_ArchR)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/', load = FALSE, overwrite = TRUE)

# peak gene expression correlation ----------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')

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
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','OPC')
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
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/', load = FALSE, overwrite = TRUE)
