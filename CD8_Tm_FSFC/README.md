
# Demo and instructions to use

## 1 System requirements

- R version 4.3.0 (2023-04-21)
- Platform: x86_64-pc-linux-gnu (64-bit)
- Running under: Ubuntu 22.04.4 LTS
- A moderately sized workstation with at least 8Gb, better 16Gb RAM can run the scripts in this folder.

---

### R package dependencies 

- DT_0.33
- ReactomePA_1.43.0
- DOSE_3.25.0
- GSEABase_1.62.0
- graph_1.77.3
- annotate_1.77.0
- XML_3.99-0.16.1
- biomaRt_2.55.4
- gridExtra_2.3
- scales_1.3.0
- tximport_1.28.0
- DESeq2_1.39.8
- SummarizedExperiment_1.29.1
- MatrixGenerics_1.11.1
- matrixStats_1.3.0
- GenomicRanges_1.51.4
- GenomeInfoDb_1.35.17
- plyr_1.8.9
- reshape2_1.4.4
- openxlsx_4.2.5.2
- survminer_0.4.9
- survival_3.7-0
- rstatix_0.7.2
- writexl_1.5.0
- igraph_2.0.3
- here_1.0.1
- pathview_1.39.0
- circlize_0.4.16
- ComplexHeatmap_2.15.3
- org.Hs.eg.db_3.17.0
- AnnotationDbi_1.64.1
- IRanges_2.33.1
- S4Vectors_0.37.7
- Biobase_2.59.0
- BiocGenerics_0.45.3
- limma_3.55.10
- clusterProfiler_4.7.1
- ggsci_3.1.0
- cowplot_1.1.3
- fgsea_1.25.2
- data.table_1.15.4
- janitor_2.2.0
- ggpubr_0.6.0
- magrittr_2.0.3
- lubridate_1.9.3
- forcats_1.0.0
- stringr_1.5.1
- dplyr_1.1.4
- purrr_1.0.2
- readr_2.1.5
- tidyr_1.3.1
- tibble_3.2.1
- ggplot2_3.5.1
- tidyverse_2.0.0

- devtools
- BiocManager
- mclust
- scales
- tibble
- multimode
- viridis
- plotly
- data.table
- tidyverse
- ggrepel
- rmarkdown
- RColorBrewer
- dplyr
- reshape2
- sBIC
- Rtsne
- ggplot2
- diptest
- LaplacesDemon
- Biobase
- diptest
- mousetrap
- readxl
- magrittr
- ggridges
- emdist
- plyr

### Bioconductor package dependencies 

- flowCore
- FlowSOM
- flowAI
- PeacoQC
- flowWorkspace
- uwot
- pheatmap
- tidyr
- ggpubr
- ggcyto
- flowVS
- preprocessCore
- flowTrans
- harmony

### Github package dependencies 

- CytoExploreRData
- CytoExploreR
- Spectre
- cyCombine 

## 2. Installation guide

- R and RStudio can be downloaded from https://posit.co/download/rstudio-desktop/
- The same URL hold install instructions for every operation system 
- package dependencies can be resolved by installing packages with install.packages(<Package-name>) directly, via BioConductor (instructions: https://www.bioconductor.org/install/) or through devtools (https://github.com/r-lib/devtools)

## 3. Demo and instructions to use

- CD8_Tm_FSFC relevant scripts and data are available at https://doi.org/10.5281/zenodo.12201093.
- Each CD8_Tm_FSFC relevant script and dependent data can be sourced and produce output including figures and tables that were produced as results of the manuscript.
- Comments are included to describe the intention of the script. Once all scripts are sourced, functions can be used independently provided the necessary parameters or following the main.R, providing parameters through the exemplary Parameters.R script.
- If included, seeds need to be left untouched to reproduce results of the manuscript.
- The expected runtime ranges from 30 min to several hours depending on file sizes, machine resources and parameters.
