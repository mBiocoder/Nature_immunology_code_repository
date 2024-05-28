#=============================================================================================================================#
# Title: Downsampling file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

consistent_subsampling <- function(flow_frames, target_nrow) {
  smallest_nrow <- min(sapply(flow_frames, nrow))
  actual_target_nrow <- ifelse(smallest_nrow < target_nrow, smallest_nrow, target_nrow)
  
  subsampling_list <- lapply(flow_frames, function(ff) {
    sample_indices <- sample(1:nrow(ff), size = actual_target_nrow, replace = FALSE)
    sample_indices
  })
  return(subsampling_list)
}

maximal_subsampling <- function(flow_frames, target_nrow) {
  largest_nrow <- max(sapply(flow_frames, nrow))
  actual_target_nrow <- ifelse(target_nrow > largest_nrow, largest_nrow, target_nrow)
  subsampling_list <- lapply(flow_frames, function(ff) {
    if (target_nrow < nrow(ff)) {
      size = target_nrow
    } else {
      size = nrow(ff)
    }
    sample_indices <- sample(1:nrow(ff), size = size, replace = FALSE)
  })
  return(subsampling_list)
}


subsample_flowframes <- function(list_ffr, target_nrow, subsampling_method, outputpath,sample_names){
  cat("Generate subsampling indices with", subsampling_method, "Subsampling method. \n")
  # Perform subsampling based on the chosen method
  subsampled_indices <- switch(
    subsampling_method,
    Consistent = consistent_subsampling(list_ffr, target_nrow),
    Maximal = maximal_subsampling(list_ffr, target_nrow),
    stop("Invalid subsampling method. Please choose either 'Consistent' or 'Maximal'.")
  )
  
  # Determine the maximum length among all columns
  max_length <- max(sapply(subsampled_indices, length), na.rm = TRUE)
  
  # Fill in missing values with NA to make all columns of the same length
  na_indices <- lapply(subsampled_indices, function(col) {
    if (length(col) < max_length) {
      c(col, rep(NA, max_length - length(col)))
    } else {col}
  })
  
  ss_filename <- paste0(outputpath, "subsampling_indices.csv")
  cat("Saving subsampling indices in", ss_filename," \n")
  
  index_df <- data.frame(na_indices)
  names(index_df) <- sample_names
  write.csv(index_df, file = ss_filename, row.names = FALSE)
  
  # Apply the subsampled indices to the flow frames
  cat("Subsampling flowframes with", subsampling_method, "Subsampling method. \n \n")
  subsampled_flow_frames <- list()  # Initialize an empty list to store results
  for (i in seq_along(flow_frames)) {
    subsampled_flow_frames[[i]] <- flow_frames[[i]][subsampled_indices[[i]], ]
    if(save_downsampling_fcs){flowCore::write.FCS(subsampled_flow_frames[[i]], 
                                        paste0(outputpath, "/", sample_names[i],
                                               "_downsampled.fcs"))
    }
  }
  return(list(ffr = subsampled_flow_frames, indices = subsampled_indices))
}


