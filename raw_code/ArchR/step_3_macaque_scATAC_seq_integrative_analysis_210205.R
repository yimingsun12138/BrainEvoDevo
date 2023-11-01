#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: macaque scATAC_seq integrative analysis                         ##
## Data: 2022.02.05                                                                ##
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
library(circlize)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')

# marker gene validation --------------------------------------------------
macaque_ATAC_ArchR <- addImputeWeights(macaque_ATAC_ArchR)

#marker express
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
  ArchRProj = macaque_ATAC_ArchR, 
  colorBy = "GeneIntegrationMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_ATAC_ArchR)
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

#cluster 1 End Mic
p1+p2+p$PECAM1+plot_layout(ncol = 3)
p1+p2+p$CX3CR1+plot_layout(ncol = 3)
p1+p2+p$PDGFRB+plot_layout(ncol = 3)
#cluster 2 OPC
#cluster 6 7 RG-3
p1+p2+p$SOX10+plot_layout(ncol = 3)
p1+p2+p$OLIG2+plot_layout(ncol = 3)
p1+p2+p$EGFR+plot_layout(ncol = 3)
#cluster 8 9 10 11 12 13 14 15 In
p1+p2+p$DLX2+plot_layout(ncol = 3)
p1+p2+p$DLX5+plot_layout(ncol = 3)
p1+p2+p$GAD2+plot_layout(ncol = 3)
p1+p2+p$GAD1+plot_layout(ncol = 3)
#cluster 3 4 6 7 RG
p1+p2+p$SOX9+plot_layout(ncol = 3)
p1+p2+p$PAX6+plot_layout(ncol = 3)
p1+p2+p$VIM+plot_layout(ncol = 3)
#cluster 5 IP
p1+p2+p$EOMES+plot_layout(ncol = 3)
p1+p2+p$PPP1R17+plot_layout(ncol = 3)
#cluster 16 17 18 19 20 21 22 23 Ex
p1+p2+p$NEUROD2+plot_layout(ncol = 3)
p1+p2+p$NEUROD6+plot_layout(ncol = 3)
p1+p2+p$TBR1+plot_layout(ncol = 3)
p1+p2+p$SATB2+plot_layout(ncol = 3)
p1+p2+p$FEZF2+plot_layout(ncol = 3)


# call peak ---------------------------------------------------------------
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

getPeakSet(macaque_ATAC_ArchR)

#add peak matrix
macaque_ATAC_ArchR <- addPeakMatrix(macaque_ATAC_ArchR)
getAvailableMatrices(macaque_ATAC_ArchR)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220109/", load = FALSE, overwrite = TRUE)

# marker peak -------------------------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')
getAvailableMatrices(macaque_ATAC_ArchR)

#calculate marker peak
markersPeaks <- getMarkerFeatures(
  ArchRProj = macaque_ATAC_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE)
markerList

peakset <- getPeakSet(macaque_ATAC_ArchR)
names(peakset) <- paste(as.character(peakset@seqnames),as.character(peakset@ranges),sep = '-')
table(duplicated(names(peakset)))

#annotate marker peak
for (i in names(markerList)) {
  names(markerList[[i]]) <- paste(as.character(markerList[[i]]@seqnames),as.character(markerList[[i]]@ranges),sep = '-')
  if(sum(names(markerList[[i]]) %in% names(peakset)) < length(names(markerList[[i]]))){
    stop('peak name match error!')
  } else{
    markerList[[i]]$nearestGene <- peakset[names(markerList[[i]])]$nearestGene
  }
  print(paste(as.character(i),'done!',sep = ' '))
}

#get peak express matrix
express_matrix <- getMatrixFromProject(ArchRProj = macaque_ATAC_ArchR,useMatrix = 'PeakMatrix')

meta_data <- express_matrix@colData
meta_data <- as.data.frame(meta_data)
feature_list <- express_matrix@rowRanges
feature_list <- paste(as.character(feature_list@seqnames),as.character(feature_list@ranges),sep = '-')
feature_list[1:10]
table(duplicated(feature_list))

express_matrix <- express_matrix@assays@data$PeakMatrix
rownames(express_matrix) <- feature_list
temp <- Seurat::CreateSeuratObject(counts = express_matrix,project = 'temp',assay = 'RNA',meta.data = meta_data)
temp <- AverageExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp <- temp$RNA
temp[1:3,1:3]
express_matrix <- temp

markerList <- markerList[c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')]
feature_list <- base::lapply(names(markerList),FUN = function(x){
  temp_marker <- markerList[[x]]
  temp_marker <- paste(as.character(temp_marker@seqnames),as.character(temp_marker@ranges),sep = '-')
  return(temp_marker)
})
for (i in 1:length(feature_list)) {
  names(feature_list[[i]]) <- rep(names(markerList)[i],times = length(feature_list[[i]]))
}
feature_list <- unlist(feature_list)

table(feature_list %in% rownames(express_matrix))

express_matrix <- express_matrix[feature_list,]
express_matrix <- log1p(express_matrix)
express_matrix <- t(scale(t(express_matrix)))
express_matrix <- express_matrix[,names(markerList)]
express_matrix <- t(express_matrix)

col_fun = colorRamp2(c(-1,0,2), c('#1E78B4','white','#E11E26'))
pdf(file = './res/step_3_fig_220205/macaque_ATAC_marker_peak_heatmap.pdf',width = 12,height = 6)
Heatmap(matrix = express_matrix,cluster_columns = FALSE,cluster_rows = FALSE,show_row_names = TRUE,show_column_names = FALSE,
        name = 'z score',col = col_fun,width = unit(10,'inches'),height = unit(5,'inches'),
        column_split = factor(names(feature_list),levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')),
        border = TRUE,row_split = factor(c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic'),
                                         levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')))
dev.off()

rm(list = c('temp','express_matrix'))
gc()

#track plot
#RG
p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)['RG'],
  upstream = 50000,
  downstream = 50000,
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')
)

pdf(file = './res/step_3_fig_220205/macaque_ATAC_RG_VIM_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$VIM)
dev.off()

#OPC
p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('SOX10','OLIG2','EGFR'),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[c('RG-3','OPC')],
  upstream = 50000,
  downstream = 50000,
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')
)

pdf(file = './res/step_3_fig_220205/macaque_ATAC_OPC_EGFR_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$EGFR)
dev.off()

pdf(file = './res/step_3_fig_220205/macaque_ATAC_OPC_OLIG2_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$OLIG2)
dev.off()

#IP
p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('EOMES','PPP1R17'),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[c('IP')],
  upstream = 50000,
  downstream = 50000,
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')
)

pdf(file = './res/step_3_fig_220205/macaque_ATAC_IP_EOMES_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$EOMES)
dev.off()

#Ex
p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('NEUROD6','NEUROD2','SATB2','FEZF2'),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[c('Ex-1','Ex-2','Ex-3','Ex-4')],
  upstream = 50000,
  downstream = 50000,
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')
)

pdf(file = './res/step_3_fig_220205/macaque_ATAC_Ex_NEUROD6_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$NEUROD6)
dev.off()

pdf(file = './res/step_3_fig_220205/macaque_ATAC_Ex_FEZF2_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$FEZF2)
dev.off()

#In
p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('DLX5','GAD2','GAD1','DLX2','LHX6','SST','SP8','NR2F2','PROX1','HTR3A','CXCL14'),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[c('InCGE','InMGE')],
  upstream = 50000,
  downstream = 50000,
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG','End/Mic')
)

pdf(file = './res/step_3_fig_220205/macaque_ATAC_In_DLX5_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$DLX5)
dev.off()

pdf(file = './res/step_3_fig_220205/macaque_ATAC_In_LHX6_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$LHX6)
dev.off()

pdf(file = './res/step_3_fig_220205/macaque_ATAC_In_SP8_peak_trackplot.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$SP8)
dev.off()

# marker peak motif enrichment --------------------------------------------
#add motif annotation
macaque_ATAC_ArchR <- addMotifAnnotations(ArchRProj = macaque_ATAC_ArchR, motifSet = "cisbp", name = "Motif",species = 'homo sapiens')

#peak motif enrichment
markersPeaks <- getMarkerFeatures(
  ArchRProj = macaque_ATAC_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = macaque_ATAC_ArchR,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

# chromvar analysis -------------------------------------------------------
#add background peaks
macaque_ATAC_ArchR <- addBgdPeaks(macaque_ATAC_ArchR)

#calculate motif deviation
macaque_ATAC_ArchR <- addDeviationsMatrix(
  ArchRProj = macaque_ATAC_ArchR, 
  peakAnnotation = "Motif",
  force = TRUE,
  matrixName = 'MotifMatrix'
)

plotVarDev <- getVarDeviations(macaque_ATAC_ArchR, name = "MotifMatrix", plot = TRUE)
plotVarDev$data
plotVarDev

#find marker
express_matrix <- getMatrixFromProject(ArchRProj = macaque_ATAC_ArchR,useMatrix = 'MotifMatrix',useSeqnames = 'z')
temp <- express_matrix@assays@data$z
temp <- exp(temp)
meta_data <- as.data.frame(express_matrix@colData,stringsAsFactors = FALSE)

express_matrix <- CreateSeuratObject(counts = temp,project = 'temp',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
rm(list = c('temp','meta_data'))
gc()

Idents(express_matrix) <- 'cell_type'
markersMotif <- FindAllMarkers(object = express_matrix,assay = 'RNA',slot = 'counts',test.use = 'wilcox',only.pos = TRUE)

#feature plot
macaque_ATAC_ArchR <- addImputeWeights(macaque_ATAC_ArchR)

p <- plotEmbedding(
  ArchRProj = macaque_ATAC_ArchR, 
  colorBy = "MotifMatrix", 
  name = c('z:NR2E1_665','z:NKX21_476','z:EOMES_788','z:NEUROD1_63','z:NEUROD6_821','z:ZNF32_237','z:HSF4_624','z:DLX5_412'), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_ATAC_ArchR)
)

p <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 4),p))

# save data ---------------------------------------------------------------
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220109/", load = FALSE, overwrite = TRUE)

