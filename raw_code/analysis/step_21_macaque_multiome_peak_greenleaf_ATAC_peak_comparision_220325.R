#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque multiome peak greenleaf ATAC peak comparision           ##
## Data: 2022.03.25                                                                ##
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
setwd('/data/User/sunym/project/Brain/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(parallel)
library(Seurat)
library(ArchR)
library(rtracklayer)
library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(ComplexHeatmap)
library(parallel)
library(patchwork)
library(ggpubr)
library(ggvenn)
library(circlize)
library(ChIPpeakAnno)
library(GenomicFeatures)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# macaque peakset liftover ------------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_peakset <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
names(macaque_peakset) <- 1:length(macaque_peakset)
temp <- rtracklayer::as.data.frame(macaque_peakset)
names(macaque_peakset) <- paste(temp$seqnames,temp$start,temp$end,sep = '-')
table(duplicated(names(macaque_peakset)))
macaque_peakset$peakname <- names(macaque_peakset)

#lift over
# rtracklayer::export.bed(object = macaque_peakset,con = './res/step_21_fig_220325/macaque_multiome_peakset_2203251121.bed')
system(command = '/data/User/sunym/software/UCSC_LiftOver/liftOver /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_multiome_peakset_2203251121.bed /data/User/sunym/project/Brain/data/reference/rheMac10ToHg38.over.chain /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_multiome_peakset_lifted_2203251121.bed /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_multiome_peakset_unlifted_2203251121.bed -minMatch=0.9')
macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_21_fig_220325/macaque_multiome_peakset_lifted_2203251121.bed')
table(duplicated(macaque_peakset_lifted$name))

#mapped ratio
temp <- data.frame(group=c('mapped','unmapped'),
                   num=c(sum(names(macaque_peakset) %in% macaque_peakset_lifted$name),sum(!(names(macaque_peakset) %in% macaque_peakset_lifted$name))))
temp$label <- paste0(temp$group,' ',round(temp$num/sum(temp$num)*100),'%')

pdf(file = './res/step_21_fig_220325/macaque_multiome_liftover_mapped_proportion_pieplot.pdf',width = 3,height = 3.5)
ggpie(data = temp,x='num',label = 'label',lab.pos = 'out',fill='group') + 
  theme(legend.position = 'none') + 
  labs(title = 'macaque multiome peakset liftover') + 
  theme(plot.title = element_text(size = 12,face = 'bold',hjust = 0.5))
dev.off()

#lifted peak width
summary(macaque_peakset_lifted@ranges@width)
temp <- data.frame(width=log(macaque_peakset_lifted@ranges@width))

pdf(file = './res/step_21_fig_220325/macaque_multiome_liftover_mapped_width_distribution_plot.pdf',width = 6,height = 4)
ggplot(data = temp,aes(x=width)) + geom_density() + xlab('log(width)') + 
  geom_vline(aes(xintercept = 6.125),color = 'red',size = 0.5) + 
  geom_vline(aes(xintercept = 6.312),color = 'red',size = 0.5)
dev.off()

#cut at 550bp
macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_21_fig_220325/macaque_multiome_peakset_lifted_2203251121.bed')
macaque_peakset_lifted <- macaque_peakset_lifted[macaque_peakset_lifted@ranges@width < 550]
table(duplicated(macaque_peakset_lifted$name))

#mapped ratio
temp <- data.frame(group=c('mapped','unmapped'),
                   num=c(sum(names(macaque_peakset) %in% macaque_peakset_lifted$name),sum(!(names(macaque_peakset) %in% macaque_peakset_lifted$name))))
temp$label <- paste0(temp$group,' ',round(temp$num/sum(temp$num)*100),'%')

pdf(file = './res/step_21_fig_220325/macaque_multiome_liftover_after_width_filter_mapped_proportion_pieplot.pdf',width = 3,height = 3.5)
ggpie(data = temp,x='num',label = 'label',lab.pos = 'out',fill='group') + 
  theme(legend.position = 'none') + 
  labs(title = 'macaque multiome peakset liftover') + 
  theme(plot.title = element_text(size = 12,face = 'bold',hjust = 0.5))
dev.off()

#lifted peak width
pdf(file = './res/step_21_fig_220325/macaque_multiome_liftover_after_width_filter_width_histgram.pdf',width = 8,height = 5)
hist(macaque_peakset_lifted@ranges@width)
dev.off()

# merge with greenleaf peak -----------------------------------------------
#export append peakset
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220309/Save-ArchR-Project.rds')
human_peakset <- getPeakSet(Greenleaf_ATAC_ArchR)
human_peakset <- GRanges(seqnames = human_peakset@seqnames,ranges = human_peakset@ranges)
macaque_peakset_lifted <- GRanges(seqnames = macaque_peakset_lifted@seqnames,ranges = macaque_peakset_lifted@ranges)
macaque_peakset_lifted <- macaque_peakset_lifted[macaque_peakset_lifted@seqnames %in% unique(human_peakset@seqnames)]
temp <- append(human_peakset,macaque_peakset_lifted)
rtracklayer::export.bed(object = temp,con = './res/step_21_fig_220325/macaque_human_append_peakset_2203261253.bed')

#merge using bedtools
system(command = 'sort -k1,1 -k2,2n /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_human_append_peakset_2203261253.bed > /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_human_append_peakset_sorted_2203261257.bed')
system(command = '/data/User/sunym/software/bedtools/bedtools merge -d 0 -i /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_human_append_peakset_sorted_2203261257.bed > /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_human_merged_peakset_2203261305.bed')

macaque_human_merged_peakset <- rtracklayer::import.bed(con = './res/step_21_fig_220325/macaque_human_merged_peakset_2203261305.bed')
table(macaque_human_merged_peakset@seqnames)

summary(macaque_human_merged_peakset@ranges@width)
hist(macaque_human_merged_peakset@ranges@width)

#rename macaque_human_merged_peakset
temp <- macaque_human_merged_peakset
temp <- rtracklayer::as.data.frame(temp)
names(macaque_human_merged_peakset) <- paste(temp$seqnames,temp$start,temp$end,sep = '-')
table(duplicated(macaque_human_merged_peakset))

rtracklayer::export.bed(object = macaque_human_merged_peakset,con = './res/step_21_fig_220325/macaque_human_merged_peakset_renamed_2203261351.bed')

# re-liftover to macaque --------------------------------------------------
system(command = '/data/User/sunym/software/UCSC_LiftOver/liftOver /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_human_merged_peakset_renamed_2203261351.bed /data/User/sunym/project/Brain/data/reference/hg38ToRheMac10.over.chain /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_human_merged_peak_relifted_2203261355.bed /data/User/sunym/project/Brain/res/step_21_fig_220325/macaque_human_merged_peak_relifted_unmapped_2203261356.bed -minMatch=0.9')
macaque_human_merged_peakset_relifted <- rtracklayer::import.bed(con = './res/step_21_fig_220325/macaque_human_merged_peak_relifted_2203261355.bed')
macaque_human_merged_peakset_relifted <- macaque_human_merged_peakset_relifted[macaque_human_merged_peakset_relifted@seqnames %in% unique(getPeakSet(macaque_multiome_ArchR)@seqnames)]
unique(macaque_human_merged_peakset_relifted@seqnames)
table(duplicated(macaque_human_merged_peakset_relifted$name))

temp <- macaque_human_merged_peakset_relifted$name
macaque_human_merged_peakset <- macaque_human_merged_peakset[temp]
temp <- macaque_human_merged_peakset@ranges@width - macaque_human_merged_peakset_relifted@ranges@width
temp <- temp/macaque_human_merged_peakset@ranges@width
summary(temp)
table(temp < -0.1)

# cut ratio at 0.1
peak_name <- c(abs(temp) < 0.1)
peak_name <- names(macaque_human_merged_peakset)[peak_name]
macaque_human_merged_peakset <- macaque_human_merged_peakset[peak_name]
macaque_human_merged_peakset_relifted <- macaque_human_merged_peakset_relifted[macaque_human_merged_peakset_relifted$name %in% peak_name]
length(macaque_human_merged_peakset) == length(macaque_human_merged_peakset_relifted)

summary(macaque_human_merged_peakset@ranges@width)
summary(macaque_human_merged_peakset_relifted@ranges@width)

#save data
saveRDS(macaque_human_merged_peakset,file = './res/step_21_fig_220325/macaque_human_merged_peak_filted_2203261528.rds')
saveRDS(macaque_human_merged_peakset_relifted,file = './res/step_21_fig_220325/macaque_human_merged_peak_relifted_filted_2203261528.rds')

macaque_human_merged_peakset_relifted <- macaque_human_merged_peakset_relifted[countOverlaps(query = macaque_human_merged_peakset_relifted,subject = macaque_human_merged_peakset_relifted) == 1]
macaque_human_merged_peakset <- macaque_human_merged_peakset[names(macaque_human_merged_peakset) %in% macaque_human_merged_peakset_relifted$name]
length(macaque_human_merged_peakset) == length(macaque_human_merged_peakset_relifted)
saveRDS(macaque_human_merged_peakset,file = './res/step_21_fig_220325/macaque_human_merged_peak_filted_chr_filted_overlap_2203290057.rds')
saveRDS(macaque_human_merged_peakset_relifted,file = './res/step_21_fig_220325/macaque_human_merged_peak_relifted_filted_chr_filted_overlap_2203290058.rds')

# peak overlap ------------------------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220309/Save-ArchR-Project.rds')

macaque_human_merged_peakset <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_filted_chr_filted_overlap_2203290057.rds')
macaque_human_merged_peakset_relifted <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_relifted_filted_chr_filted_overlap_2203290058.rds')

#peak overlap
macaque_peakset <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
names(macaque_peakset) <- paste(macaque_peakset@seqnames,as.character(macaque_peakset@ranges),sep = '-')
human_peakset <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
names(human_peakset) <- paste(human_peakset@seqnames,as.character(human_peakset@ranges),sep = '-')

macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_21_fig_220325/macaque_multiome_peakset_lifted_2203251121.bed')
temp <- list(macaque=findOverlaps(query = macaque_peakset_lifted,subject = macaque_human_merged_peakset)@to,
             human=findOverlaps(query = human_peakset,subject = macaque_human_merged_peakset)@to)
pdf(file = './res/step_21_fig_220325/macaque_human_peakset_lifted_done_overlap_vennplot.pdf',width = 5,height = 5)
ggvenn::ggvenn(data = temp,c('macaque','human'))
dev.off()

# saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "/data/User/sunym/temp/macaque_multiome_ArchR_220329", load = FALSE, overwrite = TRUE)
# saveArchRProject(ArchRProj = Greenleaf_ATAC_ArchR, outputDirectory = "/data/User/sunym/temp/Greenleaf_ATAC_ArchR_220329", load = FALSE, overwrite = TRUE)

# create new peak matrix --------------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = '/data/User/sunym/temp/macaque_multiome_ArchR_220329/Save-ArchR-Project.rds')
Greenleaf_ATAC_ArchR <- readRDS(file = '/data/User/sunym/temp/Greenleaf_ATAC_ArchR_220329/Save-ArchR-Project.rds')

macaque_human_merged_peakset <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_filted_chr_filted_overlap_2203290057.rds')
macaque_human_merged_peakset_relifted <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_relifted_filted_chr_filted_overlap_2203290058.rds')

#add User peakset
macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,
                                     peakSet = macaque_human_merged_peakset_relifted,
                                     genomeAnnotation = getGenomeAnnotation(ArchRProj = macaque_multiome_ArchR),
                                     force = TRUE)
table(countOverlaps(query = getPeakSet(macaque_multiome_ArchR),subject = getPeakSet(macaque_multiome_ArchR)))
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)

Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,
                                   peakSet = macaque_human_merged_peakset,
                                   genomeAnnotation = getGenomeAnnotation(ArchRProj = Greenleaf_ATAC_ArchR),
                                   force = TRUE)
table(countOverlaps(query = getPeakSet(Greenleaf_ATAC_ArchR),subject = getPeakSet(Greenleaf_ATAC_ArchR)))
Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR,verbose = TRUE,force = TRUE)

my_send_sms('add matrix done!')

#create peak matrix
macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)
Greenleaf_peak_matrix <- getMatrixFromProject(ArchRProj = Greenleaf_ATAC_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

meta_data <- as.data.frame(macaque_peak_matrix@colData)
temp <- macaque_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- macaque_peak_matrix@rowRanges$name
macaque_peak_matrix <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

meta_data <- as.data.frame(Greenleaf_peak_matrix@colData)
temp <- Greenleaf_peak_matrix@assays@data$PeakMatrix
rownames(temp) <- paste(Greenleaf_peak_matrix@rowRanges@seqnames,as.character(Greenleaf_peak_matrix@rowRanges@ranges),sep = '-')
Greenleaf_peak_matrix <- CreateSeuratObject(counts = temp,project = 'human',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

gene_list <- dplyr::intersect(x = rownames(macaque_peak_matrix@assays$RNA@counts),y = rownames(Greenleaf_peak_matrix@assays$RNA@counts))

#create pseudo bulk peak matrix for macaque
rep_list <- base::lapply(names(macaque_multiome_ArchR@projectMetadata$GroupCoverages$cell_type$Params$cellGroups),FUN = function(x){
  temp <- names(macaque_multiome_ArchR@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[x]])
  return(paste(x,temp,sep = '#'))
})
rep_list <- unlist(rep_list)

temp <- do.call(what = cbind,args = base::lapply(X = rep_list,FUN = function(x){
  arg <- strsplit(x = x,split = '#')[[1]]
  cell_list <- macaque_multiome_ArchR@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[arg[1]]][[arg[2]]]
  temp <- rowSums(macaque_peak_matrix@assays$RNA@counts[,cell_list])
  temp <- data.frame(temp)
  colnames(temp) <- paste('macaque',x,sep = '#')
  rownames(temp) <- rownames(macaque_peak_matrix@assays$RNA@counts)
  return(temp)
}))
temp <- temp[gene_list,]
macaque_peak_matrix <- temp

#create pseudo bulk peak matrix for Greenleaf
rep_list <- base::lapply(names(Greenleaf_ATAC_ArchR@projectMetadata$GroupCoverages$cell_type$Params$cellGroups),FUN = function(x){
  temp <- names(Greenleaf_ATAC_ArchR@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[x]])
  return(paste(x,temp,sep = '#'))
})
rep_list <- unlist(rep_list)

temp <- do.call(what = cbind,args = base::lapply(X = rep_list,FUN = function(x){
  arg <- strsplit(x = x,split = '#')[[1]]
  cell_list <- Greenleaf_ATAC_ArchR@projectMetadata$GroupCoverages$cell_type$Params$cellGroups[[arg[1]]][[arg[2]]]
  temp <- rowSums(Greenleaf_peak_matrix@assays$RNA@counts[,cell_list])
  temp <- data.frame(temp)
  colnames(temp) <- paste('human',x,sep = '#')
  rownames(temp) <- rownames(Greenleaf_peak_matrix@assays$RNA@counts)
  return(temp)
}))
temp <- temp[gene_list,]
Greenleaf_peak_matrix <- temp

#merge
peak_matrix <- cbind(macaque_peak_matrix,Greenleaf_peak_matrix)
meta_data <- data.frame(sample=colnames(peak_matrix))
meta_data$species <- unlist(base::lapply(X = meta_data$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][1])
}))
meta_data$cell_type <- unlist(base::lapply(X = meta_data$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))
rownames(meta_data) <- meta_data$sample

peak_matrix <- CreateSeuratObject(counts = peak_matrix,project = 'temp',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
# saveRDS(object = peak_matrix,file = './res/step_21_fig_220325/macaque_human_merged_pseudo_bulk_peak_matrix_Seurat.rds')

# sample distance analysis ------------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220309/Save-ArchR-Project.rds')

macaque_human_merged_peakset <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_filted_chr_filted_overlap_2203290057.rds')
macaque_human_merged_peakset_relifted <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_relifted_filted_chr_filted_overlap_2203290058.rds')

peak_matrix <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_pseudo_bulk_peak_matrix_Seurat.rds')

#re annotate
peak_matrix <- peak_matrix[,!(peak_matrix$cell_type %in% c('Cycling','End','End/Per','Ependymal','Per','Mic','OPC'))]
Idents(peak_matrix) <- 'cell_type'

marker_list <- FindAllMarkers(object = peak_matrix,assay = 'RNA',slot = 'counts',test.use = 'DESeq2',logfc.threshold = 1,only.pos = TRUE,verbose = TRUE,min.cells.group = 2)
my_send_sms('marker done!')
marker_list_ori <- marker_list

temp <- unique(marker_list$gene)
VariableFeatures(peak_matrix) <- temp

#PCA
peak_matrix <- NormalizeData(object = peak_matrix,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
peak_matrix <- ScaleData(object = peak_matrix,features = VariableFeatures(peak_matrix),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE)
peak_matrix <- RunPCA(object = peak_matrix,assay = 'RNA',features = VariableFeatures(peak_matrix),npcs = 50,verbose = TRUE)

temp <- peak_matrix@meta.data
temp <- cbind(temp,peak_matrix@reductions$pca@cell.embeddings)
ggplot(temp,aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + theme(aspect.ratio = 1)

#UMAP
peak_matrix <- RunUMAP(object = peak_matrix,dims = 1:8,reduction = 'pca')
DimPlot(object = peak_matrix,group.by = 'cell_type',label = TRUE,repel = TRUE)

temp <- peak_matrix@meta.data
temp <- cbind(temp,peak_matrix@reductions$umap@cell.embeddings)

pdf(file = './res/step_21_fig_220325/macaque_human_merged_bulk_peak_matrix_UMAP.pdf',width = 8,height = 7)
ggplot(temp,aes(x=UMAP_1,y=UMAP_2)) + 
  geom_point(aes(color=cell_type,shape=species),size = 4,alpha = 0.5) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'UMAP plot')
dev.off()

#cluster analysis
mat <- peak_matrix[,peak_matrix$cell_type %in% c('Ex-3','Ex-2','Ex-1','IP','RG','RG-1','RG-2','InCGE','InMGE')]
marker_list <- marker_list[order(marker_list$p_val_adj,decreasing = FALSE),]
temp <- base::lapply(X = c('Ex-3','Ex-2','Ex-1','IP','RG','RG-1','RG-2','InCGE','InMGE'),FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(dim(temp)[1] < 2000){
    return(temp$gene)
  }else{
    return(temp$gene[1:2000])
  }
})
temp <- unique(unlist(temp))
group <- kmeans(mat@assays$RNA@scale.data[temp,],centers = 9)$cluster
col_fun = colorRamp2(c(-2, 0, 4), c("green", "white", "red"))
col_anno <- HeatmapAnnotation(cell_type = mat$cell_type,which = 'column')

pdf(file = './res/step_21_fig_220325/macaque_human_merged_ulk_peak_matrix_heatmap.pdf',width = 14,height = 8)
Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data[temp,]),
        show_column_names = TRUE,show_row_names = FALSE,
        cluster_rows = cluster_within_group(mat = t(as.matrix(mat@assays$RNA@scale.data[temp,])),factor = group),
        cluster_columns = cluster_within_group(mat = as.matrix(mat@assays$RNA@scale.data[temp,]),factor = mat$cell_type),
        col = col_fun,height = unit(4,'inches'),width = unit(12,'inches'),name = 'z-score',top_annotation = col_anno)
dev.off()

# peak annotation ---------------------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220309/Save-ArchR-Project.rds')

macaque_human_merged_peakset <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_filted_chr_filted_overlap_2203290057.rds')
macaque_human_merged_peakset_relifted <- readRDS(file = './res/step_21_fig_220325/macaque_human_merged_peak_relifted_filted_chr_filted_overlap_2203290058.rds')

macaque_anno <- makeTxDbFromGFF(file = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
human_anno <- makeTxDbFromGFF(file = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')

#macaque peak
macaque_peakset <- getPeakSet(macaque_multiome_ArchR)
macaque_peak_anno <- genomicElementDistribution(macaque_peakset, 
                                                TxDb = macaque_anno,
                                                promoterRegion=c(upstream=2000, downstream=500),
                                                geneDownstream=c(upstream=0, downstream=2000),
                                                promoterLevel=list(
                                                  breaks = c(-2000, -1000, -500, 0, 500),
                                                  labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                             "upstream <500b", "TSS - 500b"),
                                                  colors = c("#FFE5CC", "#FFCA99", 
                                                             "#FFAD65", "#FF8E32")))

pdf(file = './res/step_21_fig_220325/macaque_peak_region_anno.pdf',width = 8,height = 6)
macaque_peak_anno$plot
dev.off()

#human peak
human_peakset <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
human_peak_anno <- genomicElementDistribution(human_peakset, 
                                              TxDb = human_anno,
                                              promoterRegion=c(upstream=2000, downstream=500),
                                              geneDownstream=c(upstream=0, downstream=2000),
                                              promoterLevel=list(
                                                breaks = c(-2000, -1000, -500, 0, 500),
                                                labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                           "upstream <500b", "TSS - 500b"),
                                                colors = c("#FFE5CC", "#FFCA99", 
                                                           "#FFAD65", "#FF8E32")))

pdf(file = './res/step_21_fig_220325/human_peak_region_anno.pdf',width = 8,height = 6)
human_peak_anno$plot
dev.off()

#overlapped peak
macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_21_fig_220325/macaque_multiome_peakset_lifted_2203251121.bed')
temp <- dplyr::intersect(x = findOverlaps(query = macaque_peakset_lifted,subject = macaque_human_merged_peakset)@to,
                         y = findOverlaps(query = human_peakset,subject = macaque_human_merged_peakset)@to)
temp <- macaque_human_merged_peakset[temp]
overlapped_peak_anno <- genomicElementDistribution(temp, 
                                                   TxDb = human_anno,
                                                   promoterRegion=c(upstream=2000, downstream=500),
                                                   geneDownstream=c(upstream=0, downstream=2000),
                                                   promoterLevel=list(
                                                     breaks = c(-2000, -1000, -500, 0, 500),
                                                     labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                                "upstream <500b", "TSS - 500b"),
                                                     colors = c("#FFE5CC", "#FFCA99", 
                                                                "#FFAD65", "#FF8E32")))

pdf(file = './res/step_21_fig_220325/overlapped_peak_region_anno.pdf',width = 8,height = 6)
overlapped_peak_anno$plot
dev.off()

#merge plot
macaque_peak_anno$peaks$peakset <- 'macaque'
human_peak_anno$peaks$peakset <- 'human'
overlapped_peak_anno$peaks$peakset <- 'overlapped'
temp <- append(macaque_peak_anno$peaks,human_peak_anno$peaks)
temp <- append(temp,overlapped_peak_anno$peaks)

get_anno <- function(peak_anno,group.by,split.by = 'peakset'){
  names(peak_anno) <- 1:length(peak_anno)
  peak_anno <- rtracklayer::as.data.frame(peak_anno)
  peak_anno <- My_Cell_Proportion(meta_data = peak_anno,group.by = group.by,split.by = 'peakset')
  colnames(peak_anno)[2] <- 'type'
  peak_anno$level <- group.by
  peak_anno <- peak_anno[peak_anno$type != 'undefined',]
  return(peak_anno)
}


geneLevel <- get_anno(peak_anno = temp,group.by = 'geneLevel')
ExonIntron <- get_anno(peak_anno = temp,group.by = 'ExonIntron')
Exons <- get_anno(peak_anno = temp,group.by = 'Exons')
temp <- rbind(geneLevel,ExonIntron,Exons)

#change factor
temp <- rbind(geneLevel,ExonIntron,Exons)
temp$level <- factor(temp$level,levels = c('geneLevel','ExonIntron','Exons'))
temp$peakset <- factor(temp$peakset,levels = c('overlapped','human','macaque'))
temp$type <- factor(temp$type,levels = c('otherExon','CDS','utr5','utr3',
                                         'intergenic','intron','exon',
                                         'distalIntergenic','geneBody','geneDownstream','promoter'))

pdf(file = './res/step_21_fig_220325/macaque_human_overlapped_peakset_region_anno_barplot.pdf',width = 8,height = 4)
ggplot(data = temp,aes(x=peakset,y=Proportion,fill=type)) + 
  geom_bar(stat="identity",position="stack",width = 0.8) + 
  coord_flip() + facet_grid(rows = 'level') + theme_cowplot() + 
  scale_fill_manual(values = c(ggsci::pal_npg()(10),'#008280FF')) + 
  theme(aspect.ratio = 0.25,
        axis.title = element_text(face = 'bold'))
dev.off()