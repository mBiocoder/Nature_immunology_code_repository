
# Demo and instructions to use

## 1 System requirements

- R version 4.3.0 (2023-04-21)
- Platform: x86_64-pc-linux-gnu (64-bit)
- Running under: Ubuntu 22.04.4 LTS
- A moderately sized workstation with at least 8Gb, better 16Gb RAM can run the scripts in this folder.

---

### R package dependencies 

- devtools_2.4.5
- BiocManager_1.30.21
- mclust_6.0.0
- scales_1.2.1
- tibble_3.2.1
- multimode_1.3.6
- viridis_0.6.2
- plotly_4.10.1
- data.table_1.14.8
- tidyverse_2.0.0
- ggrepel_0.9.3
- rmarkdown_2.24
- RColorBrewer_1.1.3
- dplyr_1.1.2
- reshape2_1.4.4
- sBIC_1.2.5
- Rtsne_0.15
- ggplot2_3.5.1
- diptest_0.76-0
- LaplacesDemon_16.1.6
- Biobase_2.60.0
- mousetrap_3.2.2
- readxl_1.4.3
- magrittr_2.0.3
- ggridges_0.5.4
- emdist_0.3-2
- plyr_1.8.8

### Bioconductor package dependencies 

- flowCore_2.10.0
- FlowSOM_2.4.0
- flowAI_1.24.0
- PeacoQC_1.18.0
- flowWorkspace_4.8.0
- uwot_0.1.14
- pheatmap_1.0.12
- tidyr_1.3.0
- ggpubr_0.6.0
- ggcyto_1.26.0
- flowVS_1.48.0
- preprocessCore_1.60.2
- flowTrans_1.48.0
- harmony_0.1.1

### Github package dependencies 

- CytoExploreRData_1.1.0
- CytoExploreR_1.1.0
- Spectre_0.6.0
- cyCombine_0.2.2 

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
