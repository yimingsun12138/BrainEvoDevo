# BrainEvoDevo

[Comparative single-cell multiome identifies evolutionary changes in neural progenitor cells during primate brain development](https://doi.org/10.1016/j.devcel.2024.10.005)

Update: 2025/02/16

## Abstract

Understanding the cellular and genetic mechanisms driving human-specific features of cortical development remains a challenge. We generated a cell-type resolved atlas of transcriptome and chromatin accessibility in the developing macaque and mouse prefrontal cortex (PFC). Comparing with published human data, our findings demonstrate that although the cortex cellular composition is overall conserved across species, progenitor cells show significant evolutionary divergence in cellular properties. Specifically, human neural progenitors exhibit extensive transcriptional rewiring in growth factor and extracellular matrix (ECM) pathways. Expression of the human-specific progenitor marker ITGA2 in the fetal mouse cortex increases the progenitor proliferation and the proportion of upper-layer neurons. These transcriptional divergences are primarily driven by altered activity in the distal regulatory elements. The chromatin regions with human-gained accessibility are enriched with human-specific sequence changes and polymorphisms linked to intelligence and neuropsychiatric disorders. Our results identify evolutionary changes in neural progenitors and putative gene regulatory mechanisms shaping primate brain evolution.

![Project Design](https://s3.ap-east-1.amazonaws.com/mrdoge-s3-bucket/share/web_resource/BrainEvoDevo/project_design.png "Project Design")

## Data Visualization

To facilitate user interaction with our data, a Shiny app has been encapsulated within a Docker container. Users can install Docker on their computing devices and execute the provided code to run the Shiny app, thereby enabling access to http://127.0.0.1:3838 via a web browser for data browsing and interaction. It is important to note that a minimum of 16 GB of free memory is required to run this container. Additionally, after the code has finished executing, a brief waiting period is advisable to allow the container to complete data loading and processing.

```
wget https://s3.ap-east-1.amazonaws.com/mrdoge-s3-bucket/share/web_resource/BrainEvoDevo/docker/BrainEvoDevo_Shiny_app.tar
docker load -i BrainEvoDevo_Shiny_app.tar
docker run -d -p 3838:3838 brainevodevo_shiny_app
```

## Data Download

The raw data can be downloaded from [GSE241429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241429).

The processed data can be downloaded by the following URL:

* [Macaque snMultiome-seq RNA modality Seurat object](https://s3.ap-east-1.amazonaws.com/mrdoge-s3-bucket/share/web_resource/BrainEvoDevo/macaque_multiome_Seurat.rds)
* [Macaque snMultiome-seq ATAC fragment files](https://s3.ap-east-1.amazonaws.com/mrdoge-s3-bucket/share/web_resource/BrainEvoDevo/macaque_snMultiome_atac_fragments.tar.gz)
* [Mouse snMultiome-seq RNA modality Seurat object](https://s3.ap-east-1.amazonaws.com/mrdoge-s3-bucket/share/web_resource/BrainEvoDevo/mouse_multiome_Seurat.rds)
* [Mouse snMultiome-seq ATAC fragment files](https://s3.ap-east-1.amazonaws.com/mrdoge-s3-bucket/share/web_resource/BrainEvoDevo/mouse_snMultiome_atac_fragments.tar.gz)
* [Species integrated transcriptome Seurat object](https://s3.ap-east-1.amazonaws.com/mrdoge-s3-bucket/share/web_resource/BrainEvoDevo/Brain_integrated_RNA_Seurat.rds)

## Code

All raw code used in this project are deposited in the [`/raw_code`](https://github.com/yimingsun12138/BrainEvoDevo/tree/main/raw_code) folder.

## Citation

If you use the data in your research, please cite:

Liu, Y., Luo, X., Sun, Y., Chen, K., Hu, T., You, B., ... & Su, B. (2024). Comparative single-cell multiome identifies evolutionary changes in neural progenitor cells during primate brain development. Developmental Cell.
