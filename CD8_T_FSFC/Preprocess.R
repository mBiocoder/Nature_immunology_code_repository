#=============================================================================================================================#
# Title: Main preprocessing file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

# PREPROCESSING

# Load required functions
source(parameter_file, chdir=TRUE)
source("Preprocessing_functions.R", chdir = TRUE)

# Create folder to store output preprocessed files & another for quality control results
pp_output_path <-paste0(results_path, "/preprocessed_files/")
qc_output_path <- switch(qc_type,
                         flowAI = file.path(results_path, "quality_control"),
                         PeacoQC = file.path(results_path, "PeacoQC_results", "fcs_files"),
                         results_path)
dir_create(c(pp_output_path,qc_output_path))

# ========================= Preprocessing ======================================

pp_flowframes <- list()

function_order <- if (pp_order == "QT") {
  list(quality_control,transform_flowframe)
} else if (pp_order == "TQ") {
  list(transform_flowframe,quality_control)
} else {
  stop("Please give preprocessing order as either QT or TQ.")
}

for (pos in 1:length(flow_frames)) {
  cat("\n")
  # List of marker positions
  marker_indices_of_interest <- explore_data(flow_frames[[pos]])$mcol #Check if true
  # Execute functions in the specified order
  for (func in function_order) {
    pp_flowframes[[pos]] <- do.call(func, list(flow_frames[[pos]], sample_names[pos]))
  }
  
  if (quant_norm){
    cat("Quantile normalizing", sample_names[pos], "\n")
    pp_flowframes[[pos]] <- quantile_norm(pp_flowframes[[pos]],
                                          marker_indices_of_interest)
  }
  if (scaling){
    cat("Scaling", sample_names[pos], "\n")
    pp_flowframes[[pos]] <- z_score_scaling(pp_flowframes[[pos]],
                                            marker_indices_of_interest)
  }
  if (save_pp_fcs){
    # cyto_channels_restrict(pp_flowframes[[pos]])
    cat("Saving preprocessed FCS file of sample", sample_names[pos],". \n")
    flowCore::write.FCS(pp_flowframes[[pos]], 
                        paste0(pp_output_path, "pp_", 
                               sample_names[pos],".fcs"))
  }
  # marker_data <- pp_flowframes[[pos]]@exprs[,marker_indices_of_interest]
  # colnames(marker_data)<- sort(markernames(pp_flowframes[[pos]][,marker_indices_of_interest]))
  if (save_pp_csv){  
    cat("Saving preprocessed expression matrix of file", sample_names[pos],"as a CSV file. \n\n")
    write.csv(pp_flowframes[[pos]]@exprs, 
              paste0(pp_output_path, "pp_",
                     sample_names[pos], ".csv"))
  }
}
if (homogenize == "after"){
  hg_pp_flowframes <- homogenize_ffl(pp_flowframes, 
                                     sample_names, 
                                     extract_ro_fcs,
                                     extract_ro_csv,
                                     markers_of_interest,
                                     markername_replacing,
                                     pp_output_path)
}
# Plot per sample
if(plot_density_per_samplenames){
  overlap_per_marker(pp_flowframes, paste0(plot_path,"pp_density_per_sample/"),
                     sample_names,target_ncell_subsampling)
}
# Plot per batch
if(plot_density_per_batch){
  overlap_per_marker(pp_flowframes, paste0(plot_path,"pp_batches/"),
                     batch_names,target_ncell_subsampling)
}
# Plot per condition
if(plot_density_per_condition){
  overlap_per_marker(pp_flowframes, paste0(plot_path,"pp_conditions/"),
                     condition_names,target_ncell_subsampling)
}
