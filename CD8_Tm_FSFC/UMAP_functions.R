#=============================================================================================================================#
# Title: UMAP functions for the Full Spectrum Flow Cytometry CD8 high salt and low salt preprocessing and analysis
# Author: Ignacio Garcia Ribelles
# Figures: extended data 9a
#=============================================================================================================================#


concatenate_exprs_mtx <- function(flowframe_list){
  expression_matrices <- lapply(flowframe_list, exprs)
  concat_ffr <- do.call(rbind, expression_matrices)
  return(concat_ffr)
}

concatenate_ff <- function(list_ffr){
  # Append a column with corresponding numbering per origin
  numbered_ffs <- lapply(seq_along(list_ffr), function(i) {
    df <- list_ffr[[i]]@exprs
    df$sample_origin <- i - 1  # Subtract 1 to start numbering from 0
  })
  # Concatenate data frames row-wise
  combined_df <- do.call(rbind, numbered_ffs)
}


umap_marker <- function(umap, matrix, file_name, df_col, marker_name){
  as.data.frame(umap)%>%
    ggplot2::ggplot(ggplot2::aes(x = "UMAP1", 
                                 y = "UMAP2"))+
    ggplot2::theme_bw()+
    ggplot2::geom_point(ggplot2::aes(color = matrix[,df_col]), shape = 20, size = .1) +
    scale_color_gradientn(colors = c("#0000FFFF","#008080","#00FF00FF","#80ff00","#FFFF00FF", "#ff8000","#FF0000FF"))+
    ggplot2::labs(x = "UMAP1",
                  y = "UMAP2",
                  subtitle = file_name,
                  color = paste0(marker_name, "          ", sep=""))
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep=""), 
                  dpi=400)
}
umap_vector <- function(umap1,umap2, vector, file_name, palette){
  ggplot2::ggplot(data = data.frame(x = umap1, 
                                    y = umap2),
                  ggplot2::aes(x = umap1, 
                               y = umap2))+
    ggplot2::theme_bw()+
    ggplot2::geom_point(ggplot2::aes(color = c(vector)), shape = 1, size = 0.8, stroke = 0.1) +
    #ggplot2::scale_colour_manual(values = palette)+
    ggplot2::labs(x = "UMAP1",
                  y = "UMAP2",
                  subtitle = file_name,
                  color =file_name) 
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep=""))
}

tsne_vector <- function(tsne, vector, file_name, palette){
  as.data.frame(tsne)%>%
    ggplot(aes(x = dimred_1, 
               y = dimred_2))+
    theme_bw()+
    geom_point(aes(color = vector), shape = 1, size = 0.8, stroke = 0.2) +
    scale_colour_manual(values = palette)+
    labs(x = "tSNE1",
         y = "tSNE2",
         subtitle = file_name,
         color =file_name) 
  ggsave(paste(plot_path, file_name, ".pdf", sep=""), 
         width= 7, height=5, dpi=800)
}
umap_vec <- function(umap, vector, file_name, palette){
  as.data.frame(umap)%>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2))+
    theme_bw()+
    geom_point(aes(color = vector), shape = 1, size = 0.8, stroke = 0.1) +
    scale_colour_manual(values = palette)+
    labs(x = "UMAP1",
         y = "UMAP2",
         subtitle = file_name,
         color =file_name) 
  ggsave(paste(plot_path, file_name, ".pdf", sep=""), 
         width= 7, height=5, dpi=800)
}


uwot_plot <- function(df, uwot_embeding, file_name, file_location, df_col){
  cat("Computing uwot_embedding for file", file_name )
  uwot_plot <- as.data.frame(uwot_embeding)%>%
    ggplot(aes(x = UMAP1, 
               y = UMAP2))+
    geom_point(aes(color = df_col), shape = 20, size = .1) +
    #scale_color_gradientn(colors = c("#0000FFFF","#00baff","#00FF00FF","#b1ff00","#FFFF00FF", "#ff9f00","#FF0000FF"))+
    scale_color_gradientn(colors = c("#0000FFFF","#008080","#00FF00FF","#80ff00","#FFFF00FF", "#ff8000","#FF0000FF"))+
    labs(x = "UMAP1",
         y = "UMAP2",
         subtitle = file_name) 
  ggsave(paste(file_location, file_name, ".png", sep=""), 
         width= 5, height=3, dpi=300)
}
