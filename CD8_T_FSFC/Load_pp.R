#=============================================================================================================================#
# Title: Loading for preprocessed files for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

# ============================ Read preprocessed ===============================

# create an empty list to store the flowFrame objects
all_preprocessed <- list()
pp_files <- list.files(path=pp_output_path, pattern=".fcs$")
print(pp_files)
pp_fcs_files <- paste(pp_output_path,pp_files, sep="")

# loop through each FCS file
if(length(pp_files)<1){
  print("The preprocessed output folder is empty, please check whether you have 
        preprocessed your files correctly. Otherwise preprocess the files by 
        running Preprocess.R before loading them.")
}else{
  for (i in 1:length(pp_fcs_files)) {
    
    # read in the FCS file without transforming nor using a maximal truncation range
    flow_frame <- read.FCS(pp_fcs_files[i], 
                           transformation = FALSE,
                           truncate_max_range = FALSE)
    
    # add the flowFrame object to the list
    all_preprocessed[[i]] <- flow_frame
  }
}


# ====================== Count cells after Preprocessing =======================
tot_pp <- c()
for (i in 1:length(all_preprocessed)){
  ncells <- nrow(all_preprocessed[[i]])
  print(paste("After preprocessing,", sort(sample_names)[i], "had",
              ncells, "cells.", sep=" "))
  tot_pp[i] <- ncells
}
