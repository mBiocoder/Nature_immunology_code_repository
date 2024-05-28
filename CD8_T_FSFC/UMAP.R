#=============================================================================================================================#
# Title: UMAP generation file for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: Extended data 9a 
#=============================================================================================================================#

# source missing libraries and files
library(viridis)
source("UMAP_functions.R", chdir=TRUE)

list_of_ffr <- pp_flowframes
samples <- sample_names
condition <- batch_names

# Create a dictionary mapping values from vector1 to vector2
cond_dictionary <- setNames(condition, samples)

# Extracting the number of rows from each flowframe's expression matrix
num_rows <- sapply(list_of_ffr, function(flowframe) nrow(exprs(flowframe)))

# Replicating each letter based on the number of rows in the corresponding flowframe
sample_procedence <- unlist(mapply(rep, samples, num_rows))

# Condition procedence using condition dictionary
condition_procedence <- cond_dictionary[sample_procedence]


df <-concatenate_exprs_mtx(list_of_ffr)
na_rows <- which(apply(df, 1, anyNA))  # Rows with at least one NA
na_cols <- which(apply(df, 2, anyNA))  # Columns with at least one NA
sample_procedence <- sample_procedence[-na_rows]
condition_procedence  <- condition_procedence[-na_rows]
df <- na.omit(df)
write.csv(df, "df.csv", row.names = FALSE)

umap_df <- uwot::umap(X=df[,1:33],n_neighbors = 15, min_dist = 0.1)
write.csv(umap_df, "umap_df.csv", row.names = FALSE)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umapdf_with_labels <- cbind(umap_df, sample_procedence, condition_procedence)
write.csv(umapdf_with_labels, "umapdf_with_labels.csv", row.names = FALSE)

plot <- ggplot(data.frame(umapdf_with_labels), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2), 
                                                   color = sample_procedence)) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Sample") +
  #scale_color_manual(values = c("deepskyblue", "brown1"))+
  guides(colour = guide_legend(override.aes = list(size=5),  ncol = 4))+
  theme_bw()

# Save the plot as a PNG file
ggsave(paste0(final_plots,"sample_procedence.png"), plot,width = 6.05, height = 2.9, units = "in")

plot <- ggplot(data.frame(umapdf_with_labels), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2), 
                   color = condition_procedence)) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Treatment") +
  #scale_color_manual(values = c("deepskyblue", "brown1"))+
  guides(fill=guide_legend(ncol=1, byrow=TRUE))+
  guides(colour = guide_legend(override.aes = list(size=5),  ncol = 1))+
  theme_bw()

# Save the plot as a PNG file
ggsave(paste0(final_plots,"condition_procedence.png"), plot,width = 4.4, height = 2.9, units = "in")

file_precedence <- rep(sample_names,each=downsample_size)
umap_vector(umap_df, file_precedence,"file_procedence", rainbow(16)) 


umap_vector(umap_df, salt_concentration,"Salt condition", c("aquamarine", "lightsalmon"))
tsne_vector(dimred_res$layout[,c(1,2)], salt_concentration,"Salt condition", c("aquamarine", "lightsalmon"))
umap_vec(tsne_dim$Y, salt_concentration,"Salt condition", c("aquamarine", "lightsalmon"))


umap_vector(umap_df, udn_vector,"Cluster group ", c("Blue", "Grey", "Red"))

metacluster_vector <- GetMetaclusters(fsom)
metacluster_origin <- metacluster_vector[subsampling_idx]
umap_vector(umap_df, metacluster_origin,"metacluster origin", rainbow(22))
print(table(udn_vector))
