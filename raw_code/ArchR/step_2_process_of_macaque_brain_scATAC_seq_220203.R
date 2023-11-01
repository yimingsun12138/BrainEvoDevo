#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: process of macaque brain scATAC_seq                             ##
## Data: 2022.02.03                                                                ##
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

#initialize ArchR
addArchRThreads(threads = 5)

#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')
batch_1_meta_data <- readRDS(file = './res/step_1_fig_220109/batch_1_meta_data_after_filter.rds')
batch_2_meta_data <- readRDS(file = './res/step_1_fig_220109/batch_2_meta_data_after_filter.rds')


# QC and doublet filter ---------------------------------------------------
cell_list <- append(rownames(batch_1_meta_data),rownames(batch_2_meta_data))
table(cell_list %in% rownames(macaque_ATAC_ArchR@cellColData))
dim(macaque_ATAC_ArchR@cellColData)
macaque_ATAC_ArchR <- macaque_ATAC_ArchR[cell_list]

macaque_ATAC_ArchR@cellColData[1:3,]
table(rownames(macaque_ATAC_ArchR@cellColData) == cell_list)

# QC plot -----------------------------------------------------------------
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

pdf(file = './res/step_2_fig_220203/TSS_enrichment_vs_nFrag_after_filter_dotplot.pdf',width = 16,height = 10)
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

pdf(file = './res/step_2_fig_220203/TSS_enrichment_after_filter_violin_plot.pdf',width = 5,height = 5)
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

pdf(file = './res/step_2_fig_220203/nFrag_distribution_after_filter_violin_plot.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_220203/fragment_length_by_sample_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_220203/fragment_length_by_batch_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_220203/TSS_accessibility_by_sample_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'batch',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_220203/TSS_accessibility_by_batch_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

# LSI dim reduction -------------------------------------------------------
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

# cluster and umap --------------------------------------------------------
macaque_ATAC_ArchR <- addClusters(
  input = macaque_ATAC_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:30, 
  force = TRUE
)

cM <- My_Cell_Proportion(meta_data = macaque_ATAC_ArchR@cellColData,group.by = 'Sample','Clusters')
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

cM <- My_Cell_Proportion(meta_data = macaque_ATAC_ArchR@cellColData,group.by = 'batch','Clusters')
ggplot(data = cM, mapping = aes(x = Clusters,y = Proportion,fill = batch))+
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
#batch seams strong

#umap
macaque_ATAC_ArchR <- addUMAP(
  ArchRProj = macaque_ATAC_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:10, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "batch", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "donor", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#batch is so strong, fuck!

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220109/", load = FALSE, overwrite = TRUE)

# annotation --------------------------------------------------------------
#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')
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
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_ATAC_ArchR)
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "batch", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#cluster 1 End Mic
p1+p2+p$PECAM1+plot_layout(ncol = 3)
p1+p2+p$CX3CR1+plot_layout(ncol = 3)
#cluster 2 OPC
#cluster 6 7 RG-3
p1+p2+p$SOX10+plot_layout(ncol = 3)
p1+p2+p$OLIG2+plot_layout(ncol = 3)
p1+p2+p$EGFR+plot_layout(ncol = 3)
#cluster 8 9 10 11 12 13 14 15 In
p1+p2+p$DLX5+plot_layout(ncol = 3)
#cluster 3 4 6 7 RG
p1+p2+p$SOX9+plot_layout(ncol = 3)
p1+p2+p$PAX6+plot_layout(ncol = 3)
p1+p2+p$VIM+plot_layout(ncol = 3)
#cluster 5 IP
p1+p2+p$EOMES+plot_layout(ncol = 3)
p1+p2+p$PPP1R17+plot_layout(ncol = 3)
#cluster 16 17 18 19 20 21 22 23 Ex

#do constrained integration
macaque_RNA_seurat <- readRDS(file = '/data/User/sunym/project/Brain/processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('Cyc-G2M','Cyc-S','InPSB'))]

temp <- macaque_RNA_seurat@assays$RNA@counts
meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_RNA','nFeature_RNA'))]
macaque_RNA_seurat <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
rm(list = c('temp','meta_data'))

group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('End','Per','OPC','Mic')]),
    ATAC = rownames(macaque_ATAC_ArchR[macaque_ATAC_ArchR$Clusters %in% c('C1','C2')]@cellColData)
  ),
  RG = SimpleList(
    RNA = colnames(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('RG-1','RG-2','RG-3')]),
    ATAC = rownames(macaque_ATAC_ArchR[macaque_ATAC_ArchR$Clusters %in% c('C3','C4','C6','C7')]@cellColData)
  ),
  In = SimpleList(
    RNA = colnames(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE')]),
    ATAC = rownames(macaque_ATAC_ArchR[macaque_ATAC_ArchR$Clusters %in% c('C8','C9','C10','C11','C12','C13','C14','C15')]@cellColData)
  ),
  Ex = SimpleList(
    RNA = colnames(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('IP','Ex-1','Ex-2','Ex-3','Ex-4')]),
    ATAC = rownames(macaque_ATAC_ArchR[macaque_ATAC_ArchR$Clusters %in% c('C5','C16','C17','C18','C19','C20','C21','C22','C23')]@cellColData)
  )
)

macaque_ATAC_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_ATAC_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = macaque_RNA_seurat,
  addToArrow = TRUE, 
  groupList = group_list,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

#final annotation
p <- plotEmbedding(
  ArchRProj = macaque_ATAC_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_ATAC_ArchR)
)

#InMGE
p1+p2+p$LHX6+plot_layout(ncol = 3)
p1+p2+p$SST+plot_layout(ncol = 3)
#cluster 9 12 13 14 15

#InCGE
p1+p2+p$SP8+plot_layout(ncol = 3)
p1+p2+p$NR2F2+plot_layout(ncol = 3)
#cluster 8 10 11

#RG cluster 3 4
#RG-3 cluster 6 7
#OPC cluster 2
#End/Mic cluster 1

macaque_ATAC_ArchR$cell_type <- as.character(macaque_ATAC_ArchR$predictedGroup)
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$Clusters %in% c('C9','C12','C13','C14','C15'),"cell_type"] <- 'InMGE'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$Clusters %in% c('C8','C10','C11'),"cell_type"] <- 'InCGE'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$Clusters %in% c('C3','C4'),"cell_type"] <- 'RG'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$Clusters %in% c('C6','C7'),"cell_type"] <- 'RG-3'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$Clusters %in% c('C2'),"cell_type"] <- 'OPC'
macaque_ATAC_ArchR@cellColData[macaque_ATAC_ArchR$Clusters %in% c('C1'),"cell_type"] <- 'End/Mic'

table(macaque_ATAC_ArchR$cell_type)

p1 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

#save data
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220109/", load = FALSE, overwrite = TRUE)
