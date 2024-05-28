#=============================================================================================================================#
# Title: Flowframe creation file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

# Credit https://gist.github.com/yannabraham/c1f9de9b23fb94105ca5
extract_exprmtx_list <- function(list_of_ffr){
  expression_mtx_list <- lapply(list_of_ffr, function(flowframe) {
    as.data.frame(exprs(flowframe))
    })
  return(expression_mtx_list)
}

generate_cellID <- function(list_of_ffr){
  expression_mtx_list <- extract_exprmtx_list(list_of_ffr)
    
  # Initialize a variable to keep track of the cell ID
  total_cells <- 0
  cell_ID <- c()
  # Loop through the list of dataframes
  for (i in seq_along(expression_mtx_list)) {
    # Get the current dataframe
    current_df <- expression_mtx_list[[i]]
    
    # Add a "cell_ID" column with unique numbering
    cell_ID[[i]] <- c((total_cells + 1):(total_cells + nrow(current_df)))
    
    # Update the total rows count
    total_cells <- total_cells + nrow(current_df)
    
    # Update the dataframe in the list
    expression_mtx_list[[i]] <- current_df
  }
  
  # Print the modified list of dataframes
  return(cell_ID)
}

from_mtx_to_ffr <- function(data_mtx, ordered_markernames, sample_name, cellID){
  print(paste0("Generating a reordered and renamed flowframe for ", sample_name,sep = ""))
  # Appending additional columns to the expression matrix
  rownames(data_mtx) <- cellID
  # Replace markernames of matrix for those provided
  #dimnames(data_mtx)[[2]] <- ordered_markernames
  dimnames(data_mtx)[[2]] <-c(ordered_markernames)
  meta <- data.frame(name=dimnames(data_mtx)[[2]], desc=dimnames(data_mtx)[[2]])
  meta$range <- apply(apply(data_mtx,2,range),2,diff)
  meta$minRange <- apply(data_mtx,2,min)
  meta$maxRange <- apply(data_mtx,2,max)
  
  #meta <- data.frame(name=dimnames(data_mtx)[[2]], desc=dimnames(data_mtx)[[2]])
  #meta$minRange <- apply(data_mtx[,1:length(ordered_markernames)],2,min)
  #meta$maxRange <- apply(data_mtx[,1:length(ordered_markernames)],2,max)
  #meta$range <- as.numeric(meta$maxRange) - as.numeric(meta$minRange)
  write.csv(x = meta,file = paste0(rw_files_path,"/",sample_name,"_metadata.csv"))
  ff <- new("flowFrame",
            exprs=data_mtx,
            parameters=Biobase::AnnotatedDataFrame(meta)
  )
  cat("New flowframe created. \n")
  markernames(ff)<-setNames(ordered_markernames, ordered_markernames)
  return(ff)
}

from_ffl_to_fSet <- function(ffl){
  flowSet(ffl)
}

