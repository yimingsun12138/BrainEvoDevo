#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: process of mouse multiome ATAC data                             ##
## Data: 2022.07.23                                                                ##
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

# create ArchRProject from Arrow file -------------------------------------
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

# add meta data -----------------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = '../processed_data/220718_summary/mouse_multiome_Seurat_220723.rds')

#add meta data
mouse_multiome_ArchR$species <- 'mouse'
mouse_multiome_ArchR$donor <- mouse_multiome_ArchR$Sample

#filter data
cell_list <- dplyr::intersect(colnames(mouse_multiome_Seurat),rownames(mouse_multiome_ArchR@cellColData))
table(mouse_multiome_Seurat@meta.data[!(colnames(mouse_multiome_Seurat) %in% cell_list),"cell_type"])
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

pdf(file = './res/step_12_fig_220723/TSS_enrichment_vs_nFrag_dotplot.pdf',width = 12,height = 5)
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

pdf(file = './res/step_12_fig_220723/TSS_enrichment_vlnplot.pdf',width = 4,height = 4)
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

pdf(file = './res/step_12_fig_220723/nFrag_vlnplot_.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_12_fig_220723/frag_length_.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = mouse_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_12_fig_220723/TSS_accessibility_by_sample.pdf',width = 7,height = 5)
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
  resolution = 1.5, 
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

mouse_multiome_ArchR$cell_type <- mouse_multiome_Seurat@meta.data[rownames(mouse_multiome_ArchR@cellColData),"cell_type"]
p1 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

scibet::Confusion_heatmap(ori = mouse_multiome_ArchR$Clusters,prd = mouse_multiome_ArchR$cell_type)

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220711/", load = FALSE, overwrite = TRUE)

# add gene expression matrix ----------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220711/')

#load seRNA
E145_1_seRNA <- import10xFeatureMatrix(input = '/data/User/sunym/project/Brain/data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E145_1')
E155_1_seRNA <- import10xFeatureMatrix(input = '/data/User/sunym/project/Brain/data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E155_1')
E155_2_seRNA <- import10xFeatureMatrix(input = '/data/User/sunym/project/Brain/data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix.h5',
                                       names = 'E155_2')

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
seRNA <- cbind(assay(E145_1_seRNA)[gene_list,],assay(E155_1_seRNA)[gene_list,],assay(E155_2_seRNA)[gene_list,])
seRNA <- seRNA[,rownames(mouse_multiome_ArchR@cellColData)]
seRNA <- SummarizedExperiment(assays = list(counts=seRNA),rowRanges = temp)

mouse_multiome_ArchR <- addGeneExpressionMatrix(input = mouse_multiome_ArchR,seRNA = seRNA,verbose = TRUE,force = TRUE)
getAvailableMatrices(mouse_multiome_ArchR)

#add meta data
mouse_multiome_ArchR$Gex_cell_type <- mouse_multiome_Seurat@meta.data[rownames(mouse_multiome_ArchR@cellColData),"cell_type"]
table(mouse_multiome_ArchR$cell_type == mouse_multiome_ArchR$Gex_cell_type)

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220711/", load = FALSE, overwrite = TRUE)

# call peak ---------------------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220711/')

#notice point:
#maxpeaks = 150k
#macs2 FDR = 0.05
#modify mincells and max cells while doing pseudobulk
#mincells: 80

#add pseudo bulk
temp <- paste(mouse_multiome_ArchR$donor,mouse_multiome_ArchR$cell_type,sep = '_')
table(temp) %>% max()
table(temp) %>% min()
mouse_multiome_ArchR <- addGroupCoverages(ArchRProj = mouse_multiome_ArchR, 
                                          groupBy = "cell_type", 
                                          minCells = 80, 
                                          maxCells = max(table(temp)), 
                                          maxFragments = 10^10, 
                                          minReplicates = 3, 
                                          maxReplicates = 3, 
                                          sampleRatio = 0.8)

temp <- as.data.frame(mouse_multiome_ArchR@projectMetadata$GroupCoverages$cell_type$coverageMetadata)

#call peak
pathToMacs2 <- '/data/User/sunym/env/MACS2/bin/macs2'

mouse_multiome_ArchR <- addReproduciblePeakSet(
  ArchRProj = mouse_multiome_ArchR, 
  groupBy = "cell_type", 
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
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220711/", load = FALSE, overwrite = TRUE)

# peak co-accessibility ---------------------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220711/')

#add co-accessibility
mouse_multiome_ArchR <- addCoAccessibility(
  ArchRProj = mouse_multiome_ArchR,
  reducedDims = "IterativeLSI"
)

#see the results
cA <- getCoAccessibility(
  ArchRProj = mouse_multiome_ArchR,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA

markerGenes  <- c(
  'Sox9','Pax6','Vim','Hopx', #RG
  'Olig2','Egfr', #OPC
  'Eomes','Ppp1r17', #IP
  'Neurod2','Neurod6','Satb2','Fezf2', #Ex
  'Dlx5','Gad2','Sst','Sp8','Lhx6' #In
)

p <- plotBrowserTrack(
  ArchRProj = mouse_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(mouse_multiome_ArchR,resolution = 2000),
  useGroups = c('InMGE','Ex-4','Ex-2','Ex-1','IP','RG')
)

#RG
grid::grid.newpage()
grid::grid.draw(p$Sox9)

grid::grid.newpage()
grid::grid.draw(p$Hopx)

grid::grid.newpage()
grid::grid.draw(p$Pax6)

grid::grid.newpage()
grid::grid.draw(p$Vim)

#OPC
grid::grid.newpage()
grid::grid.draw(p$Olig2)

grid::grid.newpage()
grid::grid.draw(p$Egfr)

#IP
grid::grid.newpage()
grid::grid.draw(p$Eomes)

grid::grid.newpage()
grid::grid.draw(p$Ppp1r17)

#Ex
grid::grid.newpage()
grid::grid.draw(p$Neurod2)

grid::grid.newpage()
grid::grid.draw(p$Neurod6)

grid::grid.newpage()
grid::grid.draw(p$Satb2)

grid::grid.newpage()
grid::grid.draw(p$Fezf2)

#In
grid::grid.newpage()
grid::grid.draw(p$Dlx5)

grid::grid.newpage()
grid::grid.draw(p$Gad2)

grid::grid.newpage()
grid::grid.draw(p$Sst)

grid::grid.newpage()
grid::grid.draw(p$Sp8)

grid::grid.newpage()
grid::grid.draw(p$Lhx6)

#save data
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220711/", load = FALSE, overwrite = TRUE)

# peak gene expression correlation ----------------------------------------
mouse_multiome_ArchR <- loadArchRProject(path = './processed_data/mouse_multiome_ArchR_220711/')

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

markerGenes  <- c(
  'Sox9','Pax6','Vim','Hopx', #RG
  'Olig2','Egfr', #OPC
  'Eomes','Ppp1r17', #IP
  'Neurod2','Neurod6','Satb2','Fezf2', #Ex
  'Dlx5','Gad2','Sst','Sp8','Lhx6' #In
)

p <- plotBrowserTrack(
  ArchRProj = mouse_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(mouse_multiome_ArchR,resolution = 1000),
  useGroups = c('InMGE','Ex-4','Ex-2','Ex-1','IP','RG')
)

#RG
grid::grid.newpage()
grid::grid.draw(p$Sox9)

grid::grid.newpage()
grid::grid.draw(p$Hopx)

grid::grid.newpage()
grid::grid.draw(p$Pax6)

grid::grid.newpage()
grid::grid.draw(p$Vim)

#OPC
grid::grid.newpage()
grid::grid.draw(p$Olig2)

grid::grid.newpage()
grid::grid.draw(p$Egfr)

#IP
grid::grid.newpage()
grid::grid.draw(p$Eomes)

grid::grid.newpage()
grid::grid.draw(p$Ppp1r17)

#Ex
grid::grid.newpage()
grid::grid.draw(p$Neurod2)

grid::grid.newpage()
grid::grid.draw(p$Neurod6)

grid::grid.newpage()
grid::grid.draw(p$Satb2)

grid::grid.newpage()
grid::grid.draw(p$Fezf2)

#In
grid::grid.newpage()
grid::grid.draw(p$Dlx5)

grid::grid.newpage()
grid::grid.draw(p$Gad2)

grid::grid.newpage()
grid::grid.draw(p$Sst)

grid::grid.newpage()
grid::grid.draw(p$Sp8)

grid::grid.newpage()
grid::grid.draw(p$Lhx6)

#save data
cell_list <- rownames(mouse_multiome_ArchR@cellColData)
saveRDS(object = cell_list,file = './res/step_12_fig_220723/cell_list.rds')
saveArchRProject(ArchRProj = mouse_multiome_ArchR, outputDirectory = "./processed_data/mouse_multiome_ArchR_220711/", load = FALSE, overwrite = TRUE)
