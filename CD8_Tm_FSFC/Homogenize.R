#=============================================================================================================================#
# Title: Homogenize input file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

translate_vector <- function(vector, dictionary){
  translated_vector <- sapply(vector, function(letter) {
    if (letter %in% names(dictionary)) {
      return(dictionary[[letter]])
    } else {
      return("Unknown")  # Handle cases where the letter is not in the dictionary
    }
  })
}

# Replace markernames 
replace_markernames <- function(flowframe, replace_map){
  marker_list <-markernames(flowframe)
  marker_pos <- explore_data(flowframe)$mcol
  for (i in 1:length(marker_list)) {
    if (marker_list[i] %in% names(replace_map)) {
      new_marker <- replace_map[[marker_list[i]]]
      print(paste0("Marker ",marker_list[i], " was replaced for ", new_marker))
      names(new_marker) <- featureNames(flowframe)[marker_pos[i]]
      markernames(flowframe)[i] <- new_marker
    }
  }
  return(flowframe)
}

reorder_exprmtx <-function(ffr, marker_name_list){
  exprs_matr <- data.frame(matrix(NA, nrow = nrow(exprs(ffr)),
                                              ncol = length(marker_name_list)))
  colnames(exprs_matr)<-marker_name_list
  for(markernam in 1:length(marker_name_list)){
    column_index <- which(explore_data(ffr)$markers 
                          == marker_name_list[markernam])
    exprs_matr[,markernam] <- exprs(ffr)[, column_index]
  }
  # reordered_ffr[[ffr]] <-flowFrame(exprs = as.matrix(exprs_matr))
  reordered_mtx <- as.matrix(exprs_matr)
  return(reordered_mtx)
}

replace_reorder <- function(ffr, markers, marker_replacing, sample_name, cellID){
  # Replace markernames
  rp_ff <- replace_markernames(ffr, marker_replacing)
  
  # Create an expression matrix with reordered markernames
  exprs_mtx <- reorder_exprmtx(rp_ff, markers)
  
  # Integrate the new expression matrix into a new flowframe
  new_ff <- from_mtx_to_ffr(exprs_mtx, markers, sample_name, cellID)
  return(new_ff)
}

reorder_flowframes <-function(ffr_list, marker_name_list, replace_mapping){
  reordered_ffr <- vector(mode = "list", length = length(ffr_list))
  for(ffr in 1:length(ffr_list)){
    print(paste0("Reordering flowframe ",ffr))
    reordered_ffr[[ffr]] <- replace_reorder(ffr_list[[ffr]], marker_name_list, 
                                            replace_mapping)
  }
  return(reordered_ffr)
}

homogenize_ffl <- function(ff_list, sample_nams, expression_fcs, expression_csv,
                           marker_name_list, replace_mapping, output_path){
  #Generate cellID for each flowframe
  cell_IDs <- generate_cellID(ff_list)
  
  # Homogenize all flowframes if necessary
  homogene_ffL <- list()
  if(length(ff_list)<= 1){
    cat("Flowframe homogenization skipped automatically")
  } else{
    for(i in 1:length(ff_list)){
      
      if(i==1){cat("\nHomogenizing all flowframes to have the same panel according to the markers of interest:\n")}
      cat("Replacing markernames of", sample_nams[i], ":\n")
      homogene_ffL[[i]] <-replace_reorder(ff_list[[i]],marker_name_list,
                                          replace_mapping, sample_nams[i],
                                          cell_IDs[[i]])
      
      if(expression_fcs){
        cat("Saving reordered and renamed", sample_nams[i],"flowframe as a FCS file.\n")
        write.FCS(homogene_ffL[[i]], paste0(output_path, "/", sample_nams[i],"_af_ro.fcs"))
      }
      if(expression_csv){
        cat("Saving reordered and renamed ", sample_nams[i]," flowframe as a CSV file.\n \n")
        write.csv(homogene_ffL[[i]]@exprs, paste0(output_path, "/", sample_nams[i], "_af_ro.csv"))
      }
    } 
  }
  return(homogene_ffL)
}


