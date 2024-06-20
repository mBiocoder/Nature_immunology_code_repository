library(reticulate)
library(sceasy)
use_condaenv("r-reticulate", required=TRUE)
loompy <- reticulate::import('loompy')
library(Seurat)

rds_data <- readRDS("../data/pancreas/pk_all.rds")
