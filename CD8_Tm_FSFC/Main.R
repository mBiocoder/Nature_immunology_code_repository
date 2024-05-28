#=============================================================================================================================#
# Title: Main file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

# Install and load packages
source("Packages.R", chdir = FALSE)

# increase memory and use cache of larger files
knitr::opts_chunk$set(cache = TRUE, warning = FALSE, message = FALSE, cache.lazy = FALSE)

# ================================= INPUT ======================================

# Load file with specific parameter selection
parameter_file <- "Parameters_Salt.R"

# Load all files through the specified parameters
source("Load.R", chdir = TRUE)

# ============================= PREPROCESSING =================================

# Preprocess the files
source("Preprocess.R", chdir = TRUE)

# In case of starting with already preprocessed files, load the preprocessed files
source("Load_pp.R")

# =============================== ANALYSE =====================================

# Reduce the dimensionality of the preprocessed data and plot the embedding
source("UMAP.R")

# Cluster the preprocessed data in FlowSOM and plot corresponding figures
source("FlowSOM.R")




