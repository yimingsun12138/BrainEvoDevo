#####################################################################################
## Project: macaque fetal Brain multiome ATAC_seq                                  ##
## Script Purpose: re_process of macaque single cell ATAC_seq data                 ##
## Data: 2022.09.12                                                                ##
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

#load macaque_ATAC_ArchR
temp <- loadArchRProject(path = './processed_data/macaque_ATAC_ArchR_220506/')

#load Arrow file
ArrowFiles <- list.files(path = './arrow_file/macaque_ATAC_data/')
ArrowFiles <- paste('./arrow_file/macaque_ATAC_data',ArrowFiles,sep = '/')
file.exists(ArrowFiles)

#create ArchR project
macaque_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

#add meta data
macaque_ATAC_ArchR$species <- 'macaque'
macaque_ATAC_ArchR$donor <- macaque_ATAC_ArchR$Sample
table(macaque_ATAC_ArchR$Sample)
macaque_ATAC_ArchR$batch <- NA
macaque_ATAC_ArchR <- macaque_ATAC_ArchR[rownames(temp@cellColData)]
macaque_ATAC_ArchR$batch <- temp@cellColData[rownames(macaque_ATAC_ArchR@cellColData),"batch"]
macaque_ATAC_ArchR$cell_type <- temp@cellColData[rownames(macaque_ATAC_ArchR@cellColData),"cell_type"]
table(macaque_ATAC_ArchR$cell_type)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/', load = FALSE, overwrite = TRUE)

# QC plot -----------------------------------------------------------------

#load data
macaque_ATAC_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/')

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

pdf(file = './res/step_14_fig_220912/TSS_enrichment_vs_nFrag_dotplot_after_filter.pdf',width = 16,height = 10)
QC_A68A+QC_A68B+QC_A84B+QC_A84C+QC_A50A+QC_A82A+QC_A82B+plot_layout(ncol = 4)
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

pdf(file = './res/step_14_fig_220912/TSS_enrichment_vlnplot_after_filter.pdf',width = 5,height = 5)
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

pdf(file = './res/step_14_fig_220912/nFrag_vlnplot_after_filter.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution by sample
p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_14_fig_220912/frag_length_after_filter_by_sample.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#fragment length distribution by batch
p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_14_fig_220912/frag_length_after_filter_by_batch.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility by sample
p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_14_fig_220912/TSS_accessibility_by_sample_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility by batch
p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_14_fig_220912/TSS_accessibility_by_batch_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/', load = FALSE, overwrite = TRUE)

# LSI dim reduction -------------------------------------------------------
#load data
macaque_ATAC_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/')

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

#batch is so strong!

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
  dimsToUse = 1:25, 
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
  dimsToUse = 1:25, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "batch", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p1+p2+p3+plot_layout(ncol = 3)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/', load = FALSE, overwrite = TRUE)

# label transfer ----------------------------------------------------------
#load data
macaque_ATAC_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/')
macaque_RNA_Seurat <- readRDS(file = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')

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
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('End','Mic','OPC','Per','VLMC')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C18','C19')]
  ),
  Ex_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','IP','RG-1','RG-2')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C1','C2','C4','C5','C6','C7','C8','C20','C21')]
  ),
  SP_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('Ex-4')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C3')]
  ),
  In_neuron = SimpleList(
    RNA = colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_ATAC_ArchR@cellColData)[macaque_ATAC_ArchR$harmonyClusters %in% c('C9','C10','C11','C12','C13','C14','C15','C16','C17')]
  )
)

macaque_ATAC_ArchR <- addImputeWeights(ArchRProj = macaque_ATAC_ArchR,reducedDims = 'harmonyLSI')
macaque_ATAC_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_ATAC_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "harmonyLSI",
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
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p1+p2+p3+plot_layout(ncol = 3)

#Cluster 7 cell type 
table(macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters == 'C7',"cell_type"])
scibet::Confusion_heatmap(ori = macaque_ATAC_ArchR$harmonyClusters,prd = macaque_ATAC_ArchR$predictedGroup)

# annotate based on clustering --------------------------------------------
marker_gene <- c('PECAM1', #End
                 'PDGFRB', #Per
                 'CX3CR1', #Mic
                 'SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3', #RG
                 'SOX10','OLIG2','EGFR', #OPC
                 'EOMES','PPP1R17', #IP
                 'NEUROD2','NEUROD6','TBR1','SATB2','FEZF2', #Ex
                 'CUX2','HS6ST3','NRP1', #Ex-1
                 'NWD2','SNTG1', #Ex-2
                 'GRIK3','TRPM3', #Ex-3
                 'DLX5','GAD2','GAD1','DLX2', #In
                 'LHX6','SST', #InMGE
                 'SP8','NR2F2' #InCGE
)

p <- plotEmbedding(
  ArchRProj = macaque_ATAC_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "harmonyUMAP",
  imputeWeights = getImputeWeights(macaque_ATAC_ArchR)
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")

#End
p1+p2+p$PECAM1+plot_layout(ncol = 3)

#Per
p1+p2+p$PDGFRB+plot_layout(ncol = 3)

#Mic
p1+p2+p$CX3CR1+plot_layout(ncol = 3)

#RG
p1+p2+p$SOX9+plot_layout(ncol = 3)
p1+p2+p$PAX6+plot_layout(ncol = 3)
p1+p2+p$VIM+plot_layout(ncol = 3)

#OPC
p1+p2+p$OLIG2+plot_layout(ncol = 3)
p1+p2+p$EGFR+plot_layout(ncol = 3)

#IP
p1+p2+p$EOMES+plot_layout(ncol = 3)
p1+p2+p$PPP1R17+plot_layout(ncol = 3)

#Ex
p1+p2+p$SOX9+plot_layout(ncol = 3)
p1+p2+p$NEUROD2+plot_layout(ncol = 3)
p1+p2+p$NEUROD6+plot_layout(ncol = 3)
p1+p2+p$CUX2+plot_layout(ncol = 3)
p1+p2+p$HS6ST3+plot_layout(ncol = 3)
p1+p2+p$NWD2+plot_layout(ncol = 3)
p1+p2+p$NRP1+plot_layout(ncol = 3)

# re annotate -------------------------------------------------------------
#final cell type
macaque_ATAC_ArchR$cell_type <- macaque_ATAC_ArchR$predictedGroup
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$cell_type %in% c('End','Mic','OPC','Per','VLMC'),"cell_type"] <- 'End/Mic'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C18'),"cell_type"] <- 'End/Mic'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C19'),"cell_type"] <- 'OPC'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C1','C2'),"cell_type"] <- 'Ex-1'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C7') & macaque_ATAC_ArchR$cell_type %in% c('Ex-3'),"cell_type"] <- 'Ex-2'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C8'),"cell_type"] <- 'Ex-2'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C4','C5','C6'),"cell_type"] <- 'Ex-3'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C3'),"cell_type"] <- 'Ex-4'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C9','C10','C11','C13','C16','C17'),"cell_type"] <- 'InMGE'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C12','C14','C15'),"cell_type"] <- 'InCGE'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$cell_type %in% c('RG-1','RG-2'),"cell_type"] <- 'RG'

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "harmonyUMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "harmonyClusters", embedding = "harmonyUMAP")
p1+p2+p3+plot_layout(ncol = 3)

#cluster based cell type
macaque_ATAC_ArchR$cell_type_by_cluster <- macaque_ATAC_ArchR$harmonyClusters
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C18'),"cell_type_by_cluster"] <- 'End/Mic'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C19'),"cell_type_by_cluster"] <- 'OPC'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C21'),"cell_type_by_cluster"] <- 'RG'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C20'),"cell_type_by_cluster"] <- 'IP'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C1','C2'),"cell_type_by_cluster"] <- 'Ex-1'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C7','C8'),"cell_type_by_cluster"] <- 'Ex-2'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C4','C5','C6'),"cell_type_by_cluster"] <- 'Ex-3'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C3'),"cell_type_by_cluster"] <- 'Ex-4'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C9','C10','C11','C13','C16','C17'),"cell_type_by_cluster"] <- 'InMGE'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$harmonyClusters %in% c('C12','C14','C15'),"cell_type_by_cluster"] <- 'InCGE'

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "harmonyUMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type_by_cluster", embedding = "harmonyUMAP")
p1+p2+p3+plot_layout(ncol = 3)

scibet::Confusion_heatmap(ori = macaque_ATAC_ArchR$cell_type_by_cluster,prd = macaque_ATAC_ArchR$cell_type)
scibet::Confusion_heatmap(ori = macaque_ATAC_ArchR$cell_type_by_cluster,prd = macaque_ATAC_ArchR$predictedGroup)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/', load = FALSE, overwrite = TRUE)

# call peak ---------------------------------------------------------------
#load data
macaque_ATAC_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/')
p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "harmonyUMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "harmonyUMAP")
p1+p2+plot_layout(ncol = 2)

#notice point:
#maxpeaks = 150k
#macs2 FDR = 0.05
#modify mincells and max cells while doing pseudobulk
#mincells: 70

#add pseudo bulk
temp <- paste(macaque_ATAC_ArchR$donor,macaque_ATAC_ArchR$cell_type,sep = '_')
table(temp) %>% max()
table(temp) %>% min()
macaque_ATAC_ArchR <- addGroupCoverages(ArchRProj = macaque_ATAC_ArchR, 
                                        groupBy = "cell_type", 
                                        minCells = 70, 
                                        maxCells = max(table(temp)), 
                                        maxFragments = 10^10, 
                                        minReplicates = 7, 
                                        maxReplicates = 7, 
                                        sampleRatio = 0.8)

temp <- as.data.frame(macaque_ATAC_ArchR@projectMetadata$GroupCoverages$cell_type$coverageMetadata)

#call peak
pathToMacs2 <- '/data/User/sunym/env/MACS2/bin/macs2'

macaque_ATAC_ArchR <- addReproduciblePeakSet(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  pathToMacs2 = pathToMacs2,
  maxPeaks = 150000,
  genomeSize = 2.7e+09,
  cutOff = 0.05,
  force = TRUE
)

my_send_sms('call peak done!')

getPeakSet(macaque_ATAC_ArchR)

#add peak matrix
macaque_ATAC_ArchR <- addPeakMatrix(macaque_ATAC_ArchR)
getAvailableMatrices(macaque_ATAC_ArchR)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = '/data/User/sunym/project/Brain/processed_data/220802_summary/macaque_ATAC_ArchR_220912/', load = FALSE, overwrite = TRUE)
