#=============================================================================================================================#
# Title: Loading file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

### ============================ LOAD INPUT DATA =========================== ###
# Import functions
source_files <- c(parameter_file,"Compare.R", "Create_ffr.R", "Homogenize.R", 
                  "Downsampling.R") # "Density_plots.R", UMAP.R
lapply(source_files, source, chdir = TRUE)

## ====================== Create output directory paths ===================== ##

results_path <- if (is.na(!file.info(project_name)$isdir)) {
  file.path(getwd(), project_name)} else {project_name }

rw_files_path <- paste0(results_path, "/raw_files/")
plot_path <-paste0(results_path,"/raw_plots/")
downsample_path <-paste0(results_path,"/downsampled_files/")

dir_create(c(results_path, downsample_path, rw_files_path, plot_path))

## ========================= Read & Check FCS files ==========================##

# List all files in input directory ending with .fcs
input_fcs_files <- list.files(path = input_file_path, 
                              pattern = ".fcs$", 
                              full.names = TRUE)

# Extract the file names from the paths
input_fcs_names <-gsub(input_file_path, "", input_fcs_files)

# Check that there are .fcs files in the input_file_directory
cat("Reading FCS-files from", input_file_path,"\n")
if (length(input_fcs_files) == 0) {
  cat("There are no FCS files in the specified input file directory,",
      "please select another folder to proceed.\n")
} else {
  # Read all fcs files into a list of flowFrame objects
  flow_frames <- c()
  for (i in 1:length(input_fcs_files)) {
    flow_frames[[i]] <- load_cytoframe_from_fcs(input_fcs_files[i],
                                           transformation = FALSE, 
                                           truncate_max_range = truncation,
                                           min.limit = min_limit)
    # Print number of cells before preprocessing
    cat("Sample", sample_names[i],"=", input_fcs_names[i], "had a total of", 
        nrow(flow_frames[[i]]), "cells.\n")
    # Extract expression matrix and plot density profile per input file
    if(extract_rw_csv){
      cat("Extracting expression matrix of input sample", sample_names[i],"into a CSV file.\n \n")
      write.csv(flow_frames[[i]]@exprs, paste0(rw_files_path, "/", sample_names[i],"_bf_ro.csv"))
    }
  }
  # Check that there are the same number of sample_names as flow_frames
  if (length(flow_frames) != length(sample_names)) {
    message("There should be the same number of file names as FCS-files in the",
            "specified directory.\n \n")
  }
  # Iterate over all pairwise combinations of flow frames and check if their markernames are equal, adding them to an existing subgroup if they are or
  # create a new one if they aren't. Finally, print the subgroups and markernames within each subgroup.
  if (length(flow_frames) > 1) {
    if(check_panels){
      panel_dict <- compare_panels(flow_frames)
      cat(sprintf("Among the %d FCS-files, there are a total of %d different", 
                  length(flow_frames),length(panel_dict)), 
          "panels. \n You can access them",
          sprintf("by typing <panel_dict[1:%d]>.\n",length(panel_dict)))
    }
  } else {print("You only inputted one file, there is nothing to compare. \n \n")}
  
  # Create metadata.csv for all files
  metdat <- data.frame(filename = paste0( sample_names,"_af_ro.fcs")
 #                      ,batch_names, condition_names
                       )
  metdat_filepath <- paste0(rw_files_path, "reordered_metadata.csv")
  write.csv(write.csv(metdat, file = metdat_filepath,row.names = FALSE))
  
  # Homogenize all flowframes if necessary
  if(homogenize == "before"){
    flow_frames <- homogenize_ffl(flow_frames,
                                  sample_names, 
                                  extract_ro_fcs,
                                  extract_ro_csv,
                                  markers_of_interest,
                                  markername_replacing,
                                  rw_files_path)
  }
  
  # Downsampling
  ss_ffr <- subsample_flowframes(flow_frames, 
                                 target_ncell_subsampling, 
                                 subsampling_method, 
                                 downsample_path,
                                 sample_names)
}
