#=============================================================================================================================#
# Title: Package installation for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

### ======================= INSTALL & LOAD PACKAGES ======================== ###

## =========================== General packages ============================= ##

# List of required packages that can be installed with install.packages()
packages = c("devtools","BiocManager","mclust","scales","tibble", "multimode", 
             "viridis","plotly","data.table",  "tidyverse", "ggrepel",
             "rmarkdown", "RColorBrewer", "dplyr", "reshape2", "sBIC","Rtsne", 
             "ggplot2","diptest", "LaplacesDemon", "Biobase", "diptest", 
             "mousetrap", "readxl", "magrittr","ggridges", 'emdist', "plyr")

# List of general required packages that are not already installed
req_packages = packages[!(packages %in% installed.packages()[,"Package"])]
if(length(req_packages)){install.packages(req_packages, dependencies = TRUE)} # install the rest
invisible(sapply(req_packages, library, character.only = TRUE)) # Load packages

## ========================== Bioconductor packages ========================= ##

# Repeat the same procedure for Bioconductor packages
biopackages = c("flowCore", "FlowSOM", "flowAI", "PeacoQC","flowWorkspace", 
                "uwot", "pheatmap","tidyr", "ggpubr", "ggcyto", "flowVS",
                "preprocessCore", "flowTrans", "preprocessCore","harmony")
req_biopackages = biopackages[!(biopackages %in% installed.packages()[,"Package"])]
if(length(req_biopackages)){BiocManager::install(req_biopackages)}
invisible(sapply(req_biopackages, library, character.only = TRUE))

## =========================== Github repositories ========================== ##
# setRepositories(ind = c(1:6, 8))
# Repeat the same procedure for Bioconductor packages
githubrepo = c("DillonHammill/CytoExploreRData","DillonHammill/CytoExploreR",
               "biosurf/cyCombine","ImmuneDynamics/Spectre")
gitpackages = c("CytoExploreRData", "CytoExploreR", "Spectre", "cyCombine")
req_github = githubrepo[!(githubrepo %in% installed.packages()[,"Package"])]
if(length(req_github)){devtools::install_github(req_github)}
invisible(sapply(gitpackages, library, character.only = TRUE))

# options(repos = BiocManager::repositories()) # required for publishing 
