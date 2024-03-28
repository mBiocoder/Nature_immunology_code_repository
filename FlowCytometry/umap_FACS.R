# Install required packages
install.packages("umap")
BiocManager::install("flowCore")

# Load required libraries
library(umap)
library(flowCore)
library(ggplot2)

# Get the names of all .fcs files in the current directory
fcs_files <- list.files(pattern = "\\.fcs$", full.names = TRUE)
print(fcs_files)

# Read FCS files
high_NaCl_files <- list.files(pattern = "hi NaCl", full.names = TRUE)
low_NaCl_files <- list.files(pattern = "low NaCl", full.names = TRUE)

high_NaCl_data <- lapply(high_NaCl_files, read.FCS)
low_NaCl_data <- lapply(low_NaCl_files, read.FCS)

# Read and convert FCS files to data frames
high_NaCl_data <- lapply(high_NaCl_files, function(file) {
  fcs <- read.FCS(file)
  exprs <- exprs(fcs)
  as.data.frame(exprs)
})

low_NaCl_data <- lapply(low_NaCl_files, function(file) {
  fcs <- read.FCS(file)
  exprs <- exprs(fcs)
  as.data.frame(exprs)
})

# Combine data
high_NaCl_data <- do.call("rbind", high_NaCl_data)
low_NaCl_data <- do.call("rbind", low_NaCl_data)
combined_data <- rbind(high_NaCl_data, low_NaCl_data)


# Perform UMAP dimensionality reduction with adjusted parameters
umap_result <- umap(combined_data, n_neighbors = 30, min_dist = 0.1, spread = 2.0)


# Create a data frame with UMAP embeddings and condition labels
umap_df <- data.frame(UMAP1 = umap_result$layout[,1],
                      UMAP2 = umap_result$layout[,2],
                      Condition = c(rep("High NaCl", nrow(high_NaCl_data)), 
                                    rep("Low NaCl", nrow(low_NaCl_data))))

# Plot UMAP colored by condition
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 0.3) +
  scale_color_manual(values = c("brown1", "#30D5C8")) + 
  theme_minimal() +
  labs(title = "UMAP colored by High and Low NaCl condition")
p

# Save the plot as a PDF file with specified height and width
ggsave("umap_facs_condition_modified.pdf", plot = p, height = 10, width = 10)

