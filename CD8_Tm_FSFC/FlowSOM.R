#=============================================================================================================================#
# Title: FlowSOM file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: Extended data 9b
#=============================================================================================================================#

# Output directory
out_dir <- paste(output_directory, "results_", project_name, "/flowSOM/", sep="")

# Number of selected metaclusters
n_meta <- 22

# Generate flowSOM
fsom <- FlowSOM::FlowSOM(input = pp_flowSet,#out_pp_files
                         scale = flowSOM_scaling,
                         colsToUse = clustering_markers_of_interest,
                         seed = seed,
                         nClus = n_meta,
                         xdim = SOM_x, 
                         ydim = SOM_y,
                         rlen = 200,
                         compensate=FALSE,
                         transform = FALSE)
saveRDS(fsom, paste0(out_dir, "fsom_", name,".rds"))

# manual_labels <- readRDS("/home/igarcia/R/03_Salt_data/results_Salt_0607/flowSOM/fsom_as_ds_ns_2.rds")
# fsom <- manual_labels

FlowSOM::FlowSOMmary(fsom = fsom,plotFile = paste0(out_dir, "fsommary_",name,".pdf"))
types <- c("counts", "percentages", "MFIs")
MFIs <- markernames_of_interest
ds_file_names <- sub(".*15_(.*).fcs", "\\1", pp_flowSet)


# Get the features
features <- FlowSOM::GetFeatures(fsom = fsom,
                                 files = out_ds_files,
                                 filenames = ds_file_names,
                                 type = types,
                                 MFI = MFIs)

# Define the groups and feature you would want to compare.
# feature <- "cluster_percentages"#"cluster_percentages"
grouplist <- list("low_salt" = ds_file_names[low_salt],
                  "high_salt" = ds_file_names[high_salt])
stat <- "fold changes"

# Compare the 2 groups of interest
stats <- FlowSOM::GroupStats(features = features[["cluster_percentages"]],
                             groups = grouplist[1:2])
metastats <- FlowSOM::GroupStats(features = features[["metacluster_percentages"]],
                                 groups = grouplist[1:2])

# Show the findings of last step on the trees
# Define the plotting variables
stat_levels <- c(paste0(names(grouplist)[2], " underrepresented compared to ",
                        names(grouplist)[1]),
                 paste0(names(grouplist)[1], " underrepresented compared to ",
                        names(grouplist)[2]),
                 "--")
colors <- c("blue", "red", "white")

udn_levels <- c("under","-", "over")
udn_stat <- stats[stat,]
udn_stat <- factor(ifelse(udn_stat < -2.5, udn_levels[1],
                          ifelse(udn_stat > 2.5, udn_levels[3],
                                 udn_levels[2])), 
                   levels = udn_levels)
udn_stat[is.na(udn_stat)] <- udn_levels[2]

list(udn_levels)

table(fsom$cluster)

# Show statistical findings on FlowSOM trees
cluster_stat <- stats[stat,]
cluster_stat <- factor(ifelse(cluster_stat < -2.5, stat_levels[1],
                              ifelse(cluster_stat > 2.5, stat_levels[2],
                                     stat_levels[3])), 
                       levels = stat_levels)


cluster_stat[is.na(cluster_stat)] <- stat_levels[3]
gr_1 <- FlowSOM::PlotStars(fsom = fsom, title = names(grouplist)[1], 
                           nodeSizes = stats[paste0("medians ", names(grouplist)[1]),], 
                           backgroundValues = cluster_stat,
                           backgroundColors = colors,
                           maxNodeSize = 2,
                           markers=clustering_markers_of_interest,
                           list_insteadof_ggarrange = TRUE)
gr_2 <- FlowSOM::PlotStars(fsom = fsom, title = names(grouplist)[2], 
                           nodeSizes = stats[paste0("medians ", names(grouplist)[2]),],
                           backgroundValues = cluster_stat,
                           backgroundColors = colors,
                           maxNodeSize = 2,
                           markers=clustering_markers_of_interest,
                           list_insteadof_ggarrange = TRUE)
ggpubr::ggarrange(plotlist = list(gr_1$tree, gr_2$tree, gr_2$starLegend, 
                                  gr_2$backgroundLegend), 
                  heights = c(3,1))
ggsave(paste0(out_dir, "fsom_star_",name,".pdf"), width = 20, height = 10)
