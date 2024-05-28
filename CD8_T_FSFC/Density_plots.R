#=============================================================================================================================#
# Title: Density plot  file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

plotMarkerDensity <- function(step_name, flowframe, output_folder) {
  # Create the output folder if it does not exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  # Get the marker names
  marker_names <- markernames(flowframe)
  
  # Loop over the markers
  for (i in 1:length(marker_names)) {
    # Get the marker expression data
    marker_data <- exprs(flowframe)[,i]
    
    # Create a density plot for the marker
    p <- ggplot2::ggplot(data = data.frame(x=marker_data), mapping=ggplot2::aes(x=x)) +
      ggplot2::geom_density() +
      ggplot2::ggtitle(marker_names[i])
    
    # Get the marker name or feature name
    if (!is.na(marker_names[i])) {
      file_name <- paste(step_name, marker_names[i], ".png", sep = "")
    } else {
      file_name <- paste(step_name, colnames(flowframe)[i], ".png", sep = "")
    }
    # Save the density plot
    file_path <- file.path(output_folder, file_name)
    ggplot2::ggsave(file_path, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  }
}

overlap_per_marker <- function(flowfr, output_folder, groupnames,downsample_size){
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  cat("\nPlotting the cell density per marker for the following groups:\n", unique(groupnames), "\n")
  for(i in 1:length(markernames(flowfr[[1]]))){
    marker_i <- c()
    group_i <- c()
    for(j in 1:length(flowfr)){
      marker_ij <- exprs(flowfr[[j]])[,i]
      group_ij <- rep(groupnames[j], length(marker_ij))
      
      marker_i <- append(marker_i, marker_ij)
      group_i <- append(group_i, group_ij)
    }
    marker_name_i <-toString(markernames(flowfr[[j]])[i])
    all_markers_i <- structure(list(marker_i,group_i), 
                               .Names = c(marker_name_i,"Group"),
                               class = "data.frame", 
                               row.names = c(NA, -length(group_i)))
    p <- ggplot2::ggplot() + 
      ggplot2::geom_density(data=all_markers_i, ggplot2::aes(x=marker_i, group=group_i,
                                           color=group_i), adjust=2) +
      ggplot2::xlab(marker_name_i) +
      ggplot2::ylab("Cell Density")+
      ggplot2::theme_classic()
    file_path <- file.path(output_folder, paste0("overlap_density_", 
                                                marker_name_i, ".png"))
    ggplot2::ggsave(file_path, plot = p, width = 10, height = 6, units = "in", dpi = 300)
    
  }
}
