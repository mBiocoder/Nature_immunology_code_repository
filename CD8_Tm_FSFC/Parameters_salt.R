#=============================================================================================================================#
# Title: Parameter file of the full spectrum flow cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: figure 2k-m, extended data 9a-b
#=============================================================================================================================#

# Input directory path (either absolute or relative to the current working directory path)
input_file_path <- "data/MM20230606  CF-CD8 Tm salt treatment live cell export/"

# Output directory path(either absolute or relative to the current working directory path)
project_name <- "Salt_new_0613/"


sample_names <- c("aB2_h", "aB4_h", "aB6_h", "aC2_h", "aC4_h", "aC6_h","aD2_l", "aD4_l", "aD6_l",
                  "bE2_l", "aE2_l", "bE3_l", "bE4_h", "aE4_l", "bE5_h", "aE6_l",
                  "bF2_l", "bF3_l", "bF4_h", "bF5_h", "bG2_l", "bG3_l", "bG4_h", "bG5_h",
                  "bH2_l", "bH3_l", "bH4_h", "bH5_h")

# List technical batches for batch correction or check
batch_names <- substr(sample_names, 1, 1)

# List biologically interesting conditions
condition_names <- substr(sample_names, nchar(sample_names), nchar(sample_names)) 

high_salt <- which(condition_names == "h")
low_salt <- which(condition_names == "l")

# Dictionary of possible markernames stored in the input FCS files (left column) 
# and the actual markernames they should be replaced with in the output (right)
markername_replacing <- list(
  `PD1 (CD279)` = 'PD1',
  `Fas (CD95)` = 'CD95',
  `CCR7 (CD197)` = 'CCR7',
  `HLA-DR` = 'HLADR'
  # Add more mappings as needed
)

# Name a subset of markers of interest from the renamed markers
markers_of_interest <- c("CXCR6",  "CD38",   "CD127",  "HLA-DR", "PD1",  "CD5", 
                         "LAG3",   "CD57",  "CD49a",  "CCR6",   "CCR7",  
                         "CXCR3",  "CD27",   "CD56",   "CTLA4",  "CD69",   "CD28", "CXCR5",
                         "CD103",  "CD62L",  "CD25" ,  "CD95",   "CD31",   "CCR4",   "ICOS")



#============================== Loading parameters ===================================
truncation = FALSE # Truncate the highest outliers according to the upper $THRESHOLD value
min_limit = NULL # Lower truncation limit

extract_rw_csv = FALSE # Takes up to 5min per file
check_panels = FALSE # Compare pannels with one another

# Plot the density plots of the raw files
plot_density_per_samplenames = FALSE # per marker for all sample names 
plot_density_per_batch = FALSE # per marker for all subgroup1 members
plot_density_per_condition = FALSE # per marker for all subgroup2 members

## Subsampling
downsampling = TRUE
# Target number of cells to subsample for dimensionality reduction
target_ncell_subsampling = NULL
# Subsampling type to choose from "Consistent"(equel to the minimum) or "Maximal"(target if larger, total if smaller)
subsampling_method = "Consistent"
save_downsampling_fcs = FALSE
plot_umap_raw_per_file = FALSE

# Homogenize all flow_frames =  remove all columns that are not "of interest" and
# order those columns to match the order of the variable markers_of_interest
# CAVE: This step currently cannot be run before quality control 
homogenize = "before" # accepts either ["before" or "after" or "None"]
extract_ro_csv = FALSE # Save the reordered expression matrix in a CSV (May take a long time)
extract_ro_fcs = TRUE # Saves reordered FCS file (May take a long time)

#======================== Preprocessing parameters =============================

### Choose transformation type: either "arcsinh" "logicle" or "biexponential"
transform_type <- "arcsinh"
cofactor_default <- 150

# Default list of cofactors = all markers will be transformed with the same cofactor
list_of_cofactors <- rep(cofactor_default, times = length(markers_of_interest))

# Choose quality control algorithm: "None","flowAI", "PeacoQC", (not yet implemented: "flowCut" or "flowClean")
qc_type <- "PeacoQC"

# Choose preprocessing order, either first quality control, then transform "QT" 
# or the other way around "TQ"
pp_order <- "TQ"

quant_norm = FALSE
scaling = FALSE
flowSOM_scaling = FALSE

save_pp_fcs = TRUE
save_pp_csv = TRUE

normalize = FALSE
norm_covar = "batch"

