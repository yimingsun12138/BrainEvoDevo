#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: ITGA2 RNA_Scope fluorescent strength analysis                   ##
## Data: 2023.06.30                                                                ##
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
setwd('/content/data/sunym/project/Brain/')
.libPaths('/content/data/sunym/software/R_lib/R_4.3.1/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(dplyr)
library(Seurat)
library(ggplot2)
library(scibet)
library(Matrix)
library(tidyverse)
library(cowplot)
library(viridis)
library(ComplexHeatmap)
library(parallel)
library(ggsignif)
library(RColorBrewer)
library(ggsci)
library(scales)
library(patchwork)
library(ggpointdensity)
library(latex2exp)
library(ArchR)
library(scales)
library(circlize)
library(ggpubr)
library(ggtext)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(BSgenome.Mmusculus.UCSC.mm10)
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(universalmotif)
library(topGO)
library(future.apply)
library(transPlotR)
library(aplot)
library(Biostrings)
library(riverplot)
library(CellChat)
library(OpenAI4R)
library(paletteer)
library(ggpattern)
library(ggrepel)
library(openxlsx)
library(readxl)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0613')

# load RNA scope data -----------------------------------------------------
RNA_Scope_stat <- './res/step_123_fig_230630/ITGA2.xlsx'
sheet_names <- excel_sheets(path = RNA_Scope_stat)
for (i in sheet_names) {
  temp <- i
  temp <- read_excel(path = RNA_Scope_stat,sheet = temp)
  assign(x = i,value = as.data.frame(temp))
}

#check data
table(human_DAPI$scale == human_ITGA2$scale)
table(human_DAPI$region == human_ITGA2$region)
table(human_DAPI$scale == human_PAX6$scale)
table(human_DAPI$region == human_PAX6$region)

table(monkey_DAPI$scale == monkey_ITGA2$scale)
table(monkey_DAPI$region == monkey_ITGA2$region)
table(monkey_DAPI$scale == monkey_PAX6$scale)
table(monkey_DAPI$region == monkey_PAX6$region)

# calculate normalized fluorescence intensity ------------------------------
human_DAPI$species <- 'human'
human_DAPI$ITGA2 <- human_ITGA2$IntDen / human_DAPI$IntDen
human_DAPI$PAX6 <- human_PAX6$IntDen / human_DAPI$IntDen

monkey_DAPI$species <- 'monkey'
monkey_DAPI$ITGA2 <- monkey_ITGA2$IntDen / monkey_DAPI$IntDen
monkey_DAPI$PAX6 <- monkey_PAX6$IntDen / monkey_DAPI$IntDen

temp <- rbind(human_DAPI,monkey_DAPI)

bar_table <- base::do.call(what = rbind,args = base::lapply(X = c('human','monkey'),FUN = function(x){
  temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('ITGA2','PAX6'),FUN = function(y){
    temp_temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('VZ','OSVZ'),FUN = function(z){
      mean_value <- mean(temp[which(temp$species == x & temp$region == z),y])
      sd_value <- sd(temp[which(temp$species == x & temp$region == z),y])
      return(data.frame(mean = mean_value,sd = sd_value,species = x,gene = y,region = z))
    }))
  }))
}))

bar_table$region <- factor(bar_table$region,levels = c('VZ','OSVZ'))

pdf(file = './res/step_123_fig_230630/ITGA2_normalized_fluorescence_intensity_barplot.pdf',width = 4.5,height = 4)
ggplot(data = bar_table,mapping = aes(x = species,y = mean,fill = gene)) + 
  facet_wrap(facets = ~ region,ncol = 2) + 
  geom_bar(stat = 'identity',position = position_dodge(width = 0.7),width = 0.6) + 
  geom_errorbar(mapping = aes(ymin = mean - 0.25*sd,ymax = mean + 0.25*sd),position = position_dodge(width = 0.7),width = 0.3,linewidth = 0.2) + 
  theme_bw() + 
  scale_fill_manual(values = c('#5fb257','#af2318')) + 
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = 'right') + 
  ylab('normalized fluorescence intensity (a.u.)')
dev.off()

p1 <- ggplot(data = bar_table,mapping = aes(x = species,y = mean,fill = gene)) + 
  facet_wrap(facets = ~ region,ncol = 2) + 
  geom_bar(stat = 'identity',position = position_dodge(width = 0.7),width = 0.6) + 
  geom_errorbar(mapping = aes(ymin = mean - 0.25*sd,ymax = mean + 0.25*sd),position = position_dodge(width = 0.7),width = 0.3,linewidth = 0.2) + 
  theme_bw() + 
  scale_fill_manual(values = c('#5fb257','#af2318')) + 
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = 'right') + 
  ylab('normalized fluorescence intensity (a.u.)')

# merge with GPX3 data ----------------------------------------------------
RNA_Scope_stat <- './res/step_123_fig_230630/GPX3.xlsx'
sheet_names <- excel_sheets(path = RNA_Scope_stat)
for (i in sheet_names) {
  temp <- i
  temp <- read_excel(path = RNA_Scope_stat,sheet = temp)
  assign(x = i,value = as.data.frame(temp))
}

#check data
table(human_DAPI$scale == human_GPX3$scale)
table(human_DAPI$region == human_GPX3$region)
table(human_DAPI$scale == human_PAX6$scale)
table(human_DAPI$region == human_PAX6$region)

table(monkey_DAPI$scale == monkey_GPX3$scale)
table(monkey_DAPI$region == monkey_GPX3$region)
table(monkey_DAPI$scale == monkey_PAX6$scale)
table(monkey_DAPI$region == monkey_PAX6$region)

#calculate normalized fluorescence intensity
human_DAPI$species <- 'human'
human_DAPI$GPX3 <- human_GPX3$IntDen / human_DAPI$IntDen
human_DAPI$PAX6 <- human_PAX6$IntDen / human_DAPI$IntDen

monkey_DAPI$species <- 'monkey'
monkey_DAPI$GPX3 <- monkey_GPX3$IntDen / monkey_DAPI$IntDen
monkey_DAPI$PAX6 <- monkey_PAX6$IntDen / monkey_DAPI$IntDen

temp <- rbind(human_DAPI,monkey_DAPI)

bar_table <- base::do.call(what = rbind,args = base::lapply(X = c('human','monkey'),FUN = function(x){
  temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('GPX3','PAX6'),FUN = function(y){
    temp_temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('VZ','OSVZ'),FUN = function(z){
      mean_value <- mean(temp[which(temp$species == x & temp$region == z),y])
      sd_value <- sd(temp[which(temp$species == x & temp$region == z),y])
      return(data.frame(mean = mean_value,sd = sd_value,species = x,gene = y,region = z))
    }))
  }))
}))

bar_table$region <- factor(bar_table$region,levels = c('VZ','OSVZ'))

p2 <- ggplot(data = bar_table,mapping = aes(x = species,y = mean,fill = gene)) + 
  facet_wrap(facets = ~ region,ncol = 2) + 
  geom_bar(stat = 'identity',position = position_dodge(width = 0.7),width = 0.6) + 
  geom_errorbar(mapping = aes(ymin = mean - 0.25*sd,ymax = mean + 0.25*sd),position = position_dodge(width = 0.7),width = 0.3,linewidth = 0.2) + 
  theme_bw() + 
  scale_fill_manual(values = c('#5fb257','#af2318')) + 
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = 'right') + 
  ylab('normalized fluorescence intensity (a.u.)')

#merge plot
p1 <- p1 + coord_cartesian(ylim = c(0,0.5))
p2 <- p2 + coord_cartesian(ylim = c(0,0.5)) + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

pdf(file = './res/step_123_fig_230630/ITGA2_GPX3_normalized_fluorescence_intensity.pdf',width = 8,height = 4)
p1+p2+plot_layout(ncol = 2)
dev.off()