#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: GPX3 RNA_Scope dot number analysis                              ##
## Data: 2023.06.23                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.0/')
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

# load RNA-Scope data -----------------------------------------------------
RNA_Scope_stat <- './res/step_122_fig_230623/GPX3.xlsx'
sheet_names <- excel_sheets(path = RNA_Scope_stat)

for (i in sheet_names) {
  temp <- i
  temp <- read_excel(path = RNA_Scope_stat,sheet = temp)
  assign(x = i,value = as.data.frame(temp))
}

# dot number per DAPI region ----------------------------------------------
human_GPX3$species <- 'human'
human_GPX3$gene <- 'GPX3'
human_GPX3$DAPI <- human_DAPI$Area * human_DAPI$percent_Area / 100

human_PAX6$species <- 'human'
human_PAX6$gene <- 'PAX6'
human_PAX6$DAPI <- human_DAPI$Area * human_DAPI$percent_Area / 100

monkey_GPX3$species <- 'monkey'
monkey_GPX3$gene <- 'GPX3'
monkey_GPX3$DAPI <- monkey_DAPI$Area * monkey_DAPI$percent_Area / 100

monkey_PAX6$species <- 'monkey'
monkey_PAX6$gene <- 'PAX6'
monkey_PAX6$DAPI <- monkey_DAPI$Area * monkey_DAPI$percent_Area / 100

temp <- rbind(human_GPX3,human_PAX6,monkey_GPX3,monkey_PAX6)
temp[which(temp$scale == 20),"DAPI"] <- temp[which(temp$scale == 20),"DAPI"] * (1000/20)^2
temp[which(temp$scale == 63),"DAPI"] <- temp[which(temp$scale == 63),"DAPI"] * (1000/63)^2
temp$Count_Per_DAPI <- temp$Count / temp$DAPI

bar_table <- base::do.call(what = rbind,args = base::lapply(X = c('human','monkey'),FUN = function(x){
  temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('GPX3','PAX6'),FUN = function(y){
    temp_temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('VZ','OSVZ'),FUN = function(z){
      mean_value <- mean(temp[which(temp$species == x & temp$gene == y & temp$region == z),"Count_Per_DAPI"])
      sd_value <- sd(temp[which(temp$species == x & temp$gene == y & temp$region == z),"Count_Per_DAPI"])
      return(data.frame(mean = mean_value,sd = sd_value,species = x,gene = y,region = z))
    }))
  }))
}))

bar_table$region <- factor(bar_table$region,levels = c('VZ','OSVZ'))

pdf(file = './res/step_122_fig_230623/GPX3_dot_count_per_DAPI_barplot.pdf',width = 4.5,height = 4)
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
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'right') + 
  ylab('fluorescent dot count per DAPI region')
dev.off()

# merge ITGA2 and GPX3 ----------------------------------------------------
p_ITGA2 <- p_ITGA2 + 
  coord_cartesian(ylim = c(0,0.0025))

p_GPX3 <- p_GPX3 + 
  coord_cartesian(ylim = c(0,0.0025)) + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pdf(file = './res/step_122_fig_230623/ITGA2_GPX3_dot_number_per_DAPI_region.pdf',width = 8,height = 4)
p_ITGA2 + p_GPX3 + plot_layout(ncol = 2)
dev.off()