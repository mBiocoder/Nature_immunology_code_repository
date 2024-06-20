#=============================================================================================================================#
# Title: scRNAseq CD8 Tn high salt and low salt preprocessing and analysis
# Author: Mahima Arunkumar
# Figures: Extended data 11 a-d
#=============================================================================================================================#


#=============================================================================================================================#
# Prepare R session
#=============================================================================================================================#

library(dplyr)
library(ggpubr)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(DoubletFinder)
library(ggplot2)
library(RColorBrewer)
library(ArchR)
library(leiden)
library(EnhancedVolcano)
library(DESeq2)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(cowplot)
library(openxlsx)
library(ggsignif)
library(presto)
library(ComplexHeatmap)
library(enrichplot)
library(rstatix)
library(vegan)
library(compositions)
library(reshape2)

################################# Preprocessing ###############################################################################
# CD8 T cells
# 0170: CD8 Tm (CD8+CD45RA-) low salt
# 0171: CD8 Tm + NaCl (CD8+CD45RA-) high salt
# 0172: CD8 Tn (CD8+CD45RA+) low salt 
# 0173: CD8 Tn + NaCl (CD8+CD45RA+) high salt
# 0174: CD8 (CD8+)

theme_custom <- theme(axis.text.x = element_text(size = 12.8, color = "black"),
                      axis.text.y = element_text(size = 12.8, color = "black"), 
                      axis.title.x = element_text(size = 16, margin = margin(6,0,0,0)), 
                      axis.title.y = element_text(size = 16, margin = margin(0,8,0,0)),
                      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(5,0,10,0)),
                      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(0,0,10,0)), 
                      axis.ticks.x = element_line(size = 0.4, colour = "black"), 
                      axis.ticks.y = element_line(size = 0.4, colour = "black"),
                      axis.line = element_blank(), #axis.line = element_line(size = 0.4, colour="black"),
                      panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
                      panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
                      panel.border = element_rect(colour = "black", size = 0.8, fill = NA),
                      panel.background = element_rect(colour = NA, fill = NA),
                      strip.text = element_text(size = 16, face = "bold", margin = margin(0,0,10,0)), 
                      strip.background = element_blank(), 
                      legend.key = element_blank(), legend.key.size = unit(0.5,"cm"),
                      legend.text = element_text(size = 12.8), 
                      legend.title = element_text(size = 12.8),
                      plot.background = element_rect(colour = NA, fill = NA), aspect.ratio = 1)


NVAR = 8000   # number of variable features to be computed
var_cutoff = 0.8 # to choose number of PCA components

# Load high and low salt data and create Seurat object
seurat_0172 <- Read10X(data.dir = "raw_data/0172/") #low salt
seurat_0173 <- Read10X(data.dir = "raw_data/0173/") # high salt
seurat_0172 <- CreateSeuratObject(counts = seurat_0172$`Gene Expression`, min.cells = 0, min.features = 0, 
                                  project = "low salt",assay="RNA")
seurat_0173 <- CreateSeuratObject(counts = seurat_0173$`Gene Expression`, min.cells = 0, min.features = 0, 
                                  project = "high salt",assay="RNA")

# Add conditions
seurat_0172$condition = "Low NaCl"
seurat_0173$condition = "High NaCl"

combined.object <- merge(seurat_0172, c(seurat_0173) , 
                         add.cell.ids = c("Low NaCl", "High NaCl"), 
                         project = "salt")

combined.object[["percent.mt"]] <- PercentageFeatureSet(combined.object, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(combined.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filtering
combined.object <- subset(combined.object, 
                          nFeature_RNA < 9000 &
                            nCount_RNA < 80000 & percent.mt < 20)

# Normalize and scale the data
combined.object <- NormalizeData(combined.object, normalization.method = "LogNormalize", scale.factor = 1e4)
combined.object <- FindVariableFeatures(combined.object, selection.method = "vst", nfeatures = NVAR)
combined.object <- ScaleData(combined.object,vars.to.regress=c("nCount_RNA", "nFeature_RNA"))
combined.object <- RunPCA(combined.object)

# Choosing dimension of the PCA
pca = combined.object[["pca"]]
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
npca <- which(cumsum(varExplained)>var_cutoff)[1]

#combined.object$salt_condition <- ifelse(test = combined.object$condition == "low salt", yes = "Low NaCl", no = "High NaCl")
combined.object <- RunUMAP(combined.object, reduction = "pca", dims = 1:npca,
                           seed.use = 1, reduction.name = "umap",umap.method="uwot")
p2<-DimPlot(combined.object, reduction = "umap", group.by = "condition",pt.size=0.1) + theme_custom + ggtitle("UMAP")
p2
ggsave("./figures/UMAP_condition.pdf",p2)

combined.object <- FindNeighbors(combined.object, reduction = "pca", k.param = 20, dims = 1:npca, do.plot = F)
combined.object<- FindClusters(combined.object, random.seed = 0)
col.ls <- ArchRPalettes[1]
col.ls <- as.vector(col.ls$stallion)
p3 <- DimPlot(combined.object, reduction = "umap", group.by = "RNA_snn_res.0.8",pt.size=0.1,label=TRUE,cols=col.ls) + theme_custom + ggtitle("leiden clusters")
p3
ggsave("./figures/UMAP_leiden_clustering_condition.pdf",p3)

# Identify the dead-cells
mito_genes <- rownames(combined.object)[grep(pattern = "^MT-",rownames(combined.object))]

p4<-DotPlot(object = combined.object, features = mito_genes,col.min = -1,col.max = 1) + 
  theme_custom + theme(axis.text.x = element_text(size = 8.8, color = "black",angle = 90)) + 
  scale_x_discrete(name ="Mitochondrial genes")+scale_y_discrete(name ="Leiden clusters")
p4

# Doublet detection and removal
sweep.res <- paramSweep(combined.object, PCs = 1:10, sct = FALSEE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2,xlab="pK",ylab="BCmetric") + theme_custom
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
annotations <- combined.object@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(optimal.pk * nrow(combined.object@meta.data)) 
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

combined.object <- doubletFinder_v3(seu = combined.object, 
                                 PCs = 1:npca, 
                                 pK = optimal.pk,
                                 nExp = nExp.poi.adj)
DF.name = colnames(combined.object@meta.data)[grepl("DF.classification", colnames(combined.object@meta.data))]

p8 <- DimPlot(combined.object, reduction = "umap",group.by = DF.name,pt.size=0.1) + theme_custom + ggtitle("Doublets")
p8
ggsave("./figures/doubletfinder_UMAP.pdf",p8, units = "cm", width = 15, height = 15)

p9 <- DimPlot(combined.object, reduction = "tsne",group.by = DF.name,pt.size=0.1) + theme_custom + ggtitle("Doublets")
p9
ggsave("./figures/doubletfinder_tSNE.pdf",p9, units = "cm", width = 15, height = 15)

combined.object = combined.object[, combined.object@meta.data[, DF.name] == "Singlet"]

# Save and load object
saveRDS(combined.object,"./output_files/combined_object_no_doublets")
combined.object <- readRDS("./combined_object_no_doublets")

# Clustering 
combined.object <- FindNeighbors(combined.object, reduction = "pca", k.param = 20, do.plot = F)
combined.object<- FindClusters(combined.object, random.seed = 0, resolution = 0.8)
combined.object.final <- RunUMAP(combined.object, reduction = "pca", dims = 1:24,
                                        seed.use = 1, 
                                        reduction.name = "umap",umap.method="uwot")

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6)) #concat two color palettes for more than 12 colors
p10<-DimPlot(combined.object.final, reduction = "umap", 
            group.by = "orig.ident",pt.size=0.1) + theme_custom + ggtitle("UMAP")
p10
ggsave("./figures/UMAP.pdf",p10, units = "cm", width = 15, height = 15)

p11<-DimPlot(combined.object.final, reduction = "tsne", 
             group.by = "orig.ident",pt.size=0.1) + theme_custom + ggtitle("tSNE")
p11
ggsave("./figures/tSNE.pdf",p11, units = "cm", width = 15, height = 15)

#Plot leiden clusters
p12 <- DimPlot(combined.object.final, reduction = "umap", pt.size=0.1, label = TRUE ) + theme_custom + ggtitle("UMAP")
p12
ggsave("./figures/UMAP_leiden.pdf",p12, units = "cm", width = 15, height = 15)

p13 <- DimPlot(combined.object.final, reduction = "tsne", pt.size=0.1, cols = mycolors) + theme_custom + ggtitle("tSNE")
p13
ggsave("./figures/tSNE_leiden.pdf",p13, units = "cm", width = 15, height = 15)

#Save and load object
saveRDS(combined.object.final,"./output_files/")
#combined.object.final <- readRDS("./output_files/MA_combined_object_after_clustering")
combined.object.final <- readRDS("./output_files/combined_object_no_doublets")


####################################### General ####################################################

# Deterine number of cell counts
salt_conditions <- combined.object.final$condition

# Count the number of cells in each condition
cell_count_by_condition <- table(salt_conditions)

# Print or access the counts
print(cell_count_by_condition) # high salt: 2987  low salt: 5442


################################################################################################################################

# Stacked bar plot (high salt/low salt percentage per cluster)
cell_clusters = dplyr::select(combined.object.final@meta.data,seurat_clusters,orig.ident)
df_rel_freq=table(cell_clusters)/as.vector(table(combined.object.final$orig.ident))*100
par(mfrow=c(1, 1), mar=c(5, 4, 5, 6))
cell_numbers = as.vector(table(combined.object.final$seurat_clusters))
no_of_clusters <- length(unique(combined.object.final$seurat_clusters))
ticks <- sprintf("%s",seq(1:no_of_clusters))
a1<-as.data.frame(df_rel_freq)
a1_low <- a1[a1$orig.ident=="low salt",]
a1_high <- a1[a1$orig.ident=="high salt",]
df= cbind(a1_high$Freq,a1_low$Freq)

p3<-barplot(t(df_rel_freq),legend = c("High NaCl","Low NaCl"),col=c("brown1","#30D5C8"),ylab="Fraction of cells (%)", xlab="Leiden clusters",ylim=c(0,70),args.legend = list(bty = "n",x = "topright"),bg = 'transparent', space= 0.5, width = c(0.3,0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)) 
text(p3,apply(df, 1, sum)+1.5 , labels = cell_numbers)
ggsave("./figures/stacked_barplot.pdf",text(p3,apply(df, 1, sum)+1.5 , labels = cell_numbers), units = "cm", width = 15, height = 15)
dev.off()


#Heatmap per leiden cluster
markers<- presto::wilcoxauc(combined.object.final, 'seurat_clusters', assay = 'data')
markers<- top_markers(markers, n = 10)
markers <- markers[, c("rank", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")] #reorder columns

all_markers <- markers %>%
  dplyr::select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>% .[!is.na(.)]

mat<- combined.object.final[["RNA"]]@data[all_markers, ] %>% as.matrix()
mat<- t(scale(t(mat)))
cluster_anno<- combined.object.final@meta.data$seurat_clusters
quantile(mat, c(0.1, 0.95)) # what's the value range in the matrix
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

cheatmap <- Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 9),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 6),
        column_title_rot = 90,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
        show_column_names = FALSE,
        use_raster = FALSE,
        raster_quality = 4)
cheatmap
ggsave("./figures/Heatmap_per_leiden_all_top_10_genes_per_cluster.pdf",cheatmap, units = "cm", width = 15, height = 20)


####################################### Differential gene expression #########################################################

# Differential Gene Expression analysis using FindMarkers 
Idents(object = combined.object.final ) <- "condition"
markers_high_vs_low <- FindMarkers(combined.object.final, ident.1 = "High NaCl", ident.2 = "Low NaCl")

sum(markers_high_vs_low$avg_log2FC < 0) # downregulated genes (403)
sum(markers_high_vs_low$avg_log2FC > 0) # upregulated genes (593)

#Write DEG table to file
write.table(markers_high_vs_low, file="./output_files/DEG_high_salt_vs_low_salt.tsv", quote=FALSE, sep='\t', col.names = NA)

#Write in Excel files (one for upreg. and one for downreg.)
markers_high_vs_low <- tibble::rownames_to_column(markers_high_vs_low, "gene")
markers_high_vs_low <- markers_high_vs_low %>% rowwise() %>% mutate(differential_expression = if_else(avg_log2FC > 0,'upregulated','downregulated'))

upregulated_all_DEGs <- markers_high_vs_low[which(markers_high_vs_low$differential_expression == "upregulated"),] 
downregulated_all_DEGs <- markers_high_vs_low[which(markers_high_vs_low$differential_expression == "downregulated"),]

write.table(upregulated_all_DEGs, './output_files/Upregulated_all_DEGs.tsv')
write.table(downregulated_all_DEGs, './output_files/Downregulated_all_DEGs.tsv')

symbols <- mapIds(org.Hs.eg.db, keys = rownames(markers_high_vs_low), keytype = "SYMBOL", column="ENSEMBL")

#Heatmap per group (high salt/low salt)

# Filter markers by p-value and log-fold change
up_markers <- markers_high_vs_low[markers_high_vs_low$p_val_adj < 0.05 & markers_high_vs_low$avg_log2FC > 0, ] %>% arrange(desc(avg_log2FC)) %>% slice(1:20) #647
down_markers <- markers_high_vs_low[markers_high_vs_low$p_val_adj < 0.05 & markers_high_vs_low$avg_log2FC < 0, ] %>% arrange(avg_log2FC) %>% slice(1:20) #290

# Heatmaps (up -and downregulated genes)
p16_a <- DoHeatmap(combined.object.final, group.by = "orig.ident", features = rownames(up_markers))  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p16_a
ggsave("./figures/Heatmap_per_salt_identity_upregulated.pdf",p16_a, units = "cm", width = 10, height = 15)

p16_b <- DoHeatmap(combined.object.final, group.by = "orig.ident", slot= "counts", features = rownames(down_markers))  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p16_b
ggsave("./figures/Heatmap_per_salt_identity_downregulated.pdf",p16_b, units = "cm", width = 10, height = 15)


all_markers <- c(rownames(up_markers), rownames(down_markers))
p16_c <- DoHeatmap(combined.object.final, group.by = "orig.ident", slot= "counts", features = all_markers)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
p16_c
ggsave("./figures/Heatmap_per_salt_identity_up_and_downregulated.pdf",p16_c, units = "cm", width = 10, height = 15)


################################## Enrichment ##################################################

# Enrichment analysis high vs. low salt
DEG_genes = subset(markers_high_vs_low, p_val_adj < 0.05) #(989)
original_gene_list <- DEG_genes$avg_log2FC
names(original_gene_list) <- rownames(DEG_genes)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
write.csv(names(gene_list),file="./output_files/genelist_sorted_for_enrichment.csv")

#All genes from the list
gene_list_up = gene_list[gene_list > 0]
gene_list_dn = gene_list[gene_list < 0]

#upregulated significant genes
enrich_GO_up_res <- enrichGO(gene = names(gene_list_up), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                      ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)

#downregulated significant genes
enrich_GO_dn_res <- enrichGO(gene = names(gene_list_dn), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                      ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff= 0.05, readable = T)

#Write enrichment tables to file
write.table(enrich_GO_up_res, file="./output_files/enrich_GO_up_res.tsv", quote=FALSE, sep='\t', col.names = NA)
write.table(enrich_GO_dn_res, file="./output_files/enrich_GO_dn_res.tsv", quote=FALSE, sep='\t', col.names = NA)

#Enrichment dotplot for only upregulated genes
selected_pathways <- c("regulation of T cell activation", "regulation of cell-cell adhesion" , "mononuclear cell differentiation" , "extrinsic apoptotic signaling pathway", "regulation of antigen receptor-mediated signaling pathway",
                       "positive regulation of nitric oxide metabolic process" , "cellular response to interleukin-1", "positive regulation of cellular amide metabolic process" , "positive regulation of leukocyte activation",
                       "negative regulation of phosphate metabolic process", "regulation of carbohydrate metabolic process", "regulation of vitamin metabolic process", "pyrimidine ribonucleoside metabolic process", "nitric oxide metabolic process", "regulation of generation of precursor metabolites and energy")

dp_up <- dotplot(enrich_GO_up_res, showCategory=selected_pathways) + ggtitle("Overrepresentation: GO term \n (Biological process) for upregulated genes")
ggsave("./figures/Enrichment_for_upreg_signif_genes_BP_GO.pdf",dp_up, units = "cm", width = 15, height = 20)

#Enrichment dotplot for only downregulated genes
dp_dn <- dotplot(enrich_GO_dn_res, showCategory=15) + ggtitle("Overrepresentation: GO term \n (Biological process) for downregulated genes")
ggsave("./figures/Enrichment_for_downreg_signif_genes_BP_GO.pdf",dp_dn, units = "cm", width = 15, height = 20)

gse <- gseGO(geneList=gene_list, 
             ont ="BP",
             nPerm = 10000,
             keyType = "SYMBOL",
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "BH")

pdot <- dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign,scales="free",
                                                                  space="free_y") + theme(axis.text.y = element_text(size = 4, color = "black"),
                                                                                          axis.text.x = element_text(size = 3, color = "black"),axis.title.x = element_text(size = 10, margin = margin(6,0,0,0)),
                                                                                          panel.background = element_rect(fill='transparent'),                                                                                       plot.background = element_rect(fill='transparent', color=NA),) + ggtitle("") + theme(plot.title = element_text(hjust = 0.5))

# Enrichment analysis for cluster 2,3,10,11 vs rest
Idents(object = combined.object.final ) <- "seurat_clusters"
cluster_markers_high_vs_low <- FindMarkers(combined.object.final, ident.1 = c("2","3","10","11"), ident.2 = c("1", "4", "5", "6", "7", "8", "9", "12", "13", "14"))
cluster_DEG_genes = subset(cluster_markers_high_vs_low, p_val_adj < 0.05) #(781)
cluster_original_gene_list <- cluster_DEG_genes$avg_log2FC
names(cluster_original_gene_list) <- rownames(cluster_DEG_genes)
gene_list<-na.omit(cluster_original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
write.csv(names(gene_list),file="./output_files/genelist_sorted_for_enrichment_for_clusters.csv")

gene_list_up = gene_list[gene_list > 0]
gene_list_dn = gene_list[gene_list < 0]
sum(cluster_DEG_genes$avg_log2FC < 0) # downregulated genes (103)
sum(cluster_DEG_genes$avg_log2FC > 0) # upregulated genes (678)

#upregulated significant genes
enrich_GO_up_res <- enrichGO(gene = names(gene_list_up), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                             ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = T)

#downregulated significant genes
enrich_GO_dn_res <- enrichGO(gene = names(gene_list_dn), keyType = 'SYMBOL', OrgDb = org.Hs.eg.db, 
                             ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff= 0.05, readable = T)

#Write enrichment tables to file
write.table(enrich_GO_up_res, file="./output_files/enrich_clusters_GO_up_res.tsv", quote=FALSE, sep='\t', col.names = NA)
write.table(enrich_GO_dn_res, file="./output_files/enrich_clusters_GO_dn_res.tsv", quote=FALSE, sep='\t', col.names = NA)

#Enrichment dotplot for only upregulated genes
p21 <- dotplot(enrich_GO_up_res, showCategory=15) + ggtitle("Overrepresentation: GO term (Biological process) \n for upregulated genes ([2,3,10,11] vs. rest)")
ggsave("./figures/Enrichment_clusters_for_upreg_signif_genes_BP_GO.pdf",p21, units = "cm", width = 17, height = 20)

#Enrichment dotplot for only downregulated genes
p22 <- dotplot(enrich_GO_dn_res, showCategory=15) + ggtitle("Overrepresentation: GO term (Biological process) \n for downregulated genes ([2,3,10,11] vs. rest)")
ggsave("./figures/Enrichment_clusters_for_downreg_signif_genes_BP_GO.pdf",p22, units = "cm", width = 17, height = 20)


#Heatmap for DEG for clusters [2,3,10,11] vs. rest
# Check if cells belong to high salt cluster or low salt cluster
combined.object.final$salt_label <- ifelse(test = combined.object.final$seurat_clusters %in% c(2,3,10,11), yes = "high salt cluster", no = "low salt cluster")

clusters_up_markers <- cluster_markers_high_vs_low[cluster_markers_high_vs_low$p_val_adj < 0.05 & cluster_markers_high_vs_low$avg_log2FC > 0, ] %>% arrange(desc(avg_log2FC)) %>% slice(1:20)
clusters_down_markers <- cluster_markers_high_vs_low[cluster_markers_high_vs_low$p_val_adj < 0.05 & cluster_markers_high_vs_low$avg_log2FC < 0, ] %>% arrange(avg_log2FC) %>% slice(1:20) 

clusterheat_a <- DoHeatmap(combined.object.final, group.by = "salt_label", features = rownames(clusters_up_markers), size = 3)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
clusterheat_a
ggsave("./figures/Heatmap_per_salt_cluster_identity_upregulated.pdf",clusterheat_a, units = "cm", width = 15, height = 15)

clusterheat_b <- DoHeatmap(combined.object.final, group.by = "salt_label", features = rownames(clusters_down_markers), size = 3)  + scale_fill_gradientn(colors = c("blue", "white", "red"))
clusterheat_b
ggsave("./figures/Heatmap_per_salt_cluster_identity_downregulated.pdf",clusterheat_b, units = "cm", width = 15, height = 15)


######################## Compute stats for boxplots ##################################

#combined.object.final$salt_condition <- combined.object.final$condition

# List of genes to plot
genes <- c("ORAI1") #Extended data 11c

# Function to extract statistics from a boxplot
extract_stats <- function(data, group, feature) {
  group_data <- data %>% filter(!!sym(group) == !!feature)
  n <- nrow(group_data)
  min_val <- min(group_data$expression)
  max_val <- max(group_data$expression)
  q1 <- quantile(group_data$expression, 0.25)
  median <- quantile(group_data$expression, 0.5)
  q3 <- quantile(group_data$expression, 0.75)
  iqr <- IQR(group_data$expression)
  lower_whisker <- max(min(group_data$expression), q1 - 1.5 * iqr)
  upper_whisker <- min(max(group_data$expression), q3 + 1.5 * iqr)
  return(c(n, min_val, max_val, q1, median, q3, lower_whisker, upper_whisker))
}

# Loop over each gene and create plots and statistics
for (gene in genes) {
  combined.object.final <- AddModuleScore(combined.object.final,
                                          features = list(gene),
                                          name=paste0("module_", gene))
  
  # Extract statistics for both conditions
  high_salt_stats <- extract_stats(combined.object.final@meta.data %>% mutate(expression = get(paste0("module_", gene, "1"))), "salt_condition", "High NaCl")
  low_salt_stats <- extract_stats(combined.object.final@meta.data %>% mutate(expression = get(paste0("module_", gene, "1"))), "salt_condition", "Low NaCl")
  
  # Perform Wilcoxon rank-sum test
  wilcox_test <- combined.object.final@meta.data %>%
    mutate(expression = get(paste0("module_", gene, "1"))) %>%
    wilcox_test(expression ~ salt_condition, alternative = "greater")
  
  # Get the p-value
  p_value <- wilcox_test$p
  
  # Prepare the output data
  output_data <- data.frame(
    n = c(high_salt_stats[1], low_salt_stats[1]),
    type = c("High NaCl", "Low NaCl"),
    statistical_test = "one-sided",
    min = c(high_salt_stats[2], low_salt_stats[2]),
    max = c(high_salt_stats[3], low_salt_stats[3]),
    q1 = c(high_salt_stats[4], low_salt_stats[4]),
    median = c(high_salt_stats[5], low_salt_stats[5]),
    q3 = c(high_salt_stats[6], low_salt_stats[6]),
    lower_whisker = c(high_salt_stats[7], low_salt_stats[7]),
    upper_whisker = c(high_salt_stats[8], low_salt_stats[8]),
    p_value = p_value
  )
  
  # Write the output to a CSV file
  write.csv(output_data, paste0("statistics_summary_", gene, ".csv"), row.names = FALSE)
  
  # Generate and save the plots
  p <- plot_grid(
    FeaturePlot(object = combined.object.final, features = paste0("module_", gene, "1"), 
                pt.size = 0.4, cols = c("blue", "red"), 
                combine = TRUE) + ggtitle(paste("Module Score", gene)) +
      theme(plot.title = element_text(size = 10)) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
    VlnPlot(object = combined.object.final, features = paste0("module_", gene, "1"), 
            group.by = "salt_condition") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +  
      ggtitle(paste("Module Score", gene)) +
      theme(plot.title = element_text(size = 10)) + stat_compare_means(method.args = list(alternative = "two.sided")) + 
      geom_boxplot(width=0.8, size=0.5)
  )
  
  ggsave(paste0("./Feature_and_Violinplot_", gene, "_high_vs_low_salt.pdf"), p, units = "cm", width = 27, height = 15)
}


# Extract statistics function
extract_stats <- function(data, group, feature) {
  group_data <- data %>% filter(!!sym(group) == !!feature)
  n <- nrow(group_data)
  min_val <- min(group_data$expression)
  max_val <- max(group_data$expression)
  q1 <- quantile(group_data$expression, 0.25)
  median <- quantile(group_data$expression, 0.5)
  q3 <- quantile(group_data$expression, 0.75)
  iqr <- IQR(group_data$expression)
  lower_whisker <- max(min(group_data$expression), q1 - 1.5 * iqr)
  upper_whisker <- min(max(group_data$expression), q3 + 1.5 * iqr)
  return(c(n, min_val, max_val, q1, median, q3, lower_whisker, upper_whisker))
}

# List of pathways and their respective genes
pathways <- list(
  activation_GO_genes = c("ZP4","SMARCB1","IL6ST","CD46","RASGRP1","TESPA1","CD274","CD74","THY1","CD160","DNAJA3","NFKBID","TNFSF14","LILRB1","PRKCQ","PIK3CA","SELENOK","RARA","HLX","HLA-DQB2","JAK3","ICOS","ARID1B","RAG1","SLC7A1","CCR7","LILRB4","HSPD1","SMARCD1","BCL10","ABL1","ABL2","DPP4","KLHL25","DHPS","CD28","NCK2","PTPN11","PNP","SOCS5","DOCK8","ICOSLG","MALT1","TMIGD2","CD1D","PCK1","TNFSF13B","GLI3","GLI2","CD86","PTPN22","VNN1","ZBTB1","HLA-DRB3","IHH","HES1","ZMIZ1","IL6","VAV1","HHLA2","PYCARD","CD27","RASAL3","TRAF6","IGFBP2","KITLG","BTN2A2","CD80","CARD11","BCL6","RPS3","NKAP","CCDC88B","CD55","CCL19","SHH","SHB","KLRK1","FYN","LCK","MDK","PHF10","KAT5","CLECL1P","IL7R","SMARCD2","BAD","ZBTB16","PRKCZ","TNFRSF14","AKT1","EFNB3","IGF2","GPAM","ZAP70","ZP3","NFKBIZ","NCKAP1L","WNT10B","IL2RA","EPO","SMARCC1","DUSP10","SYK","CCR2","CD5","IL15","XCL1","EFNB2","ACTL6A","RHOA","HLA-DPB1","HLA-A","PBRM1","HLA-DQA2","HLA-DRB1","HLA-DQB1","SMARCA4","SMARCA2","B2M","MAP3K8","TNFSF9","LGALS9","CSK","HSPH1","NOD2","CD4","IL21","HLA-DRB5","STAT5B","HAVCR2","CORO1A","AGER","SIRPB1","AP3B1","CBFB","IL1B","IL1A","IFNG","FOXP3","LILRB2","CCL21","AIF1","HLA-DOA","IL2RG","RHOH","CD47","HLA-DMB","HLA-DMA","LEP","IL1RL2","RUNX3","PTPRC","IL12B","IL12A","CD70","TGFBR2","RUNX1","CD3E","CD83","PIK3R6","PPP3CA","SPTA1","NLRP3","EFNB1","ADA","RIPK2","SOX4","VSIR","LYN","YES1","VCAM1","CD209","ITPKB","PTPN6","CR1","PDCD1LG2","TYK2","EGR3","BRD7","SASH3","IL23A","XBP1","HLA-G","CD81","GATA3","FCHO1","SLAMF1","ARID1A","JAK2","IL4","HLA-DPA1","SOX13","CD6","FLOT2","IGF1","EBI3","AMBRA1","SMARCE1","HLA-DOB","HLA-DRB4","HLA-E","IL2","TNFSF11","ANXA1","ACTL6B","NCK1","SIRPG","HLA-DQA1","HLA-DRA","SART1","TNFSF4","VTCN1","IL7","IL18","CD24","SMARCC2","BMI1","IL36B","CAV1","IL4I1","CD276","TNFRSF13C","FADD","TFRC","CD40LG","HMGB1","IL4R","ZBTB7B","SRC","SOX12","ARID2","AP3D1","LGALS1","ACTB","SMARCD3","CYRIB","ADAM8","SIRPA","SOCS1","CCL5","CCL2","LEF1","FOXO3","IL23R","IL12RB1"),
  nfat_genes = c("CIB1", "CIB1", "PPP3CC", "PPP3CC", "SPPL3", "PTBP1", "PPP3CA", "PPP3CA", "LACRT", "SLC9A1", "AKAP6", "PPP3R2", "PPP3R2", "CABIN1", "AKAP5", "TMEM62", "PPP3CB", "PPP3CB", "PPP3R1", "PPP3R1", "TNF", "CALHM1", "LMCD1", "ERBB3", "IGF1", "CMLC1", "CMLC1"),
  TCR_sig_GO_genes = c("RELA","CCR7","BCL10","CD81","TESPA1","RAB29","CARD11","RPS3","TRAT1","IKBKG","LCK","CD226","PRKD2","NECTIN2","ADA","KCNN4"),
  MAPK_GO_genes = c("MYD88","MAPK10","ARHGEF6","LGALS9","ERCC6","SH2D3C","WNT5A","NOD2","MAP4K3","CASR","IRAK1","AGT","MAPK11","MAPK8IP2","LRRK2","MAPK1","ZFP36L1","MINK1","DUSP10","MAP4K1","ZFP36","TNF","CDC42EP5","DUSP9","PAFAH1B1","IL1B","MAP3K5","NOX1","IKBKB","MAP3K12","MAP4K2","MAP2K3","MAPK8IP1","PBK","ATF2","MAP2K4","MAPK9","MAPK8","TNFSF11","MAP2K7","TNFRSF19","TLR4","TLR3","STRADB","IRAK4","MAPK13","STK3","PTGER4","CRYAB","SH2D3A","RIPK2","NPHS1","GPS1","TLR7","SMAD3","MAP3K20","CARD9","MAP3K11","MAPK14","CCM2","MAP2K6","MAP3K7","MAPK3","DAXX","MAPKAPK2","HACD3","NFKB1","GPS2","CRKL","TRIB1","MAP3K10","MAP3K13","NOD1","TAOK2"),
  cytotoxicity_GO_genes = c("HLA-C", "CD1D", "CD1E", "SLC22A13", "HLA-F", "MR1", "P2RX7", "CD1A", "ULBP3", "XCL1", "HFE", "ULBP1", "ULBP2", "HLA-A", "HLA-DRB1", "RAET1E", "B2M", "NECTIN2", "TAP2", "PTPRC", "IL12B", "IL12A", "RAET1L", "IL23A", "HLA-G", "CD1C", "CD1B", "MICA", "MICB", "HLA-E", "HLA-DRA", "HLA-H", "HLA-B", "PVR", "STX7", "FADD", "RAET1G", "CYRIB", "IL23R", "IL12RB1"),
  effector_GO_genes = c("GZMB", "GZMK", "IFNG", "EOMES", "ITGAD", "ITGAX", "ITGB7", "CXCR3", "CCR5"),
  TRM_GO_genes = c("XIST", "UBC", "LGALS3", "MT-CO2", "VIM", "ANKRD28", "RGS1", "RGCC", "HSPA1B", "MT-ND4", "HSP90AB1", "PPP1R15A")
)


process_pathway <- function(pathway_genes, pathway_name) {
  combined.object.final <- AddModuleScore(combined.object.final,
                                          features = list(pathway_genes),
                                          name = paste0("module_", pathway_name))
  
  # Extract statistics for both conditions
  high_salt_stats <- extract_stats(combined.object.final@meta.data %>% mutate(expression = get(paste0("module_", pathway_name, "1"))), "salt_condition", "High NaCl")
  low_salt_stats <- extract_stats(combined.object.final@meta.data %>% mutate(expression = get(paste0("module_", pathway_name, "1"))), "salt_condition", "Low NaCl")
  
  # Perform Wilcoxon rank-sum test
  wilcox_test <- combined.object.final@meta.data %>%
    mutate(expression = get(paste0("module_", pathway_name, "1"))) %>%
    wilcox_test(expression ~ salt_condition, alternative = "greater")
  
  # Get the p-value
  p_value <- wilcox_test$p
  
  # Boxplot
  p <- combined.object.final@meta.data %>%
    mutate(expression = get(paste0("module_", pathway_name, "1"))) %>%
    ggplot(aes(x = salt_condition, y = expression, fill = salt_condition)) +
    geom_boxplot() +
    ggtitle(paste(pathway_name, "Pathway Activity Score")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  
  list(high_salt_stats = high_salt_stats, low_salt_stats = low_salt_stats, p_value = p_value, plot = p, pathway_name = pathway_name)
}

# Process all pathways
results <- lapply(names(pathways), function(pathway_name) {
  process_pathway(pathways[[pathway_name]], pathway_name)
})

# Display results
for (result in results) {
  print(result$plot)
  print(paste("High NaCl stats:", paste(result$high_salt_stats, collapse = ", ")))
  print(paste("Low NaCl stats:", paste(result$low_salt_stats, collapse = ", ")))
  print(paste("P-value:", result$p_value))
  print(paste("Pathway:", result$pathway_name))
}



################# Additional: Violinplots and Featureplots for genes and genesets of interest ##########################################

#PGK1 gene
VlnPlot(combined.object.final, features = "PGK1", group.by = "condition", pt.size = 0.1) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width = 0.4, size = 0.1)


#NFAT5 gene
VlnPlot(combined.object.final, features = "NFAT5", group.by = "condition", pt.size = 0.1) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width = 0.4, size = 0.1)

#Check if PGK1 expression is higher in NFAT5 expressing cells or low NFAT5 expressing cells

# Plot a histogram for NFAT5 expression
data_nfat5 <- FetchData(combined.object.final, vars = 'NFAT5')
hist(data_nfat5$NFAT5, col = "grey", main = "NFAT5 Expression Histogram", xlab = "NFAT5 Expression")

# Add a vertical line at the cutoff (adjust the value)
cutoff_plot <- 0.2
abline(v = cutoff_plot, col = "red", lwd = 3)  # lwd value to make the line thicker
legend("topright", legend = c("Cutoff Threshold"), fill = "red", border = "red")

# Create a binary column indicating "high" or "low" NFAT5 expression
cutoff <- 0.0
combined.object.final$NFAT5_expression <- ifelse(data_nfat5$NFAT5 > cutoff, "High NFAT5", "Low NFAT5")

# Check the levels of the NFAT5_High column
table(combined.object.final$NFAT5_expression)

# Subset Seurat object based on high NFAT5 expression
high_nfat5 <- subset(x = combined.object.final, subset = NFAT5_High == "High")

# Subset Seurat object based on low NFAT5 expression
low_nfat5 <- subset(x = combined.object.final, subset = NFAT5_High == "Low")

#PGK1 gene
VlnPlot(combined.object.final, features = "PGK1", group.by = "NFAT5_expression", pt.size = 0.1) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width = 0.4, size = 0.1)


# CD5 gene
VlnPlot(combined.object.final, features = "CD5", group.by = "condition", pt.size = 0.1) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width = 0.4, size = 0.1)


# MTOR gene
v5 <- VlnPlot(combined.object.final, features = "MTOR", group.by = "condition", pt.size = 0.1) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width = 0.4, size = 0.1)
v5
ggsave("./figures/Violinplot_MTOR_expression.pdf",v5, units = "cm", width = 15, height = 15)


#Electron transport chain ETC geneset from Gene Ontology (https://amigo.geneontology.org/amigo/term/GO:0022900)
etc_GO_genes <- c("NDUFA4", "NDUFA7", "NDUFB2", "NDUFB4", "NDUFA3", "SRD5A1", "NDUFA6", "P4HA2", "SOD2", "NDUFB6", "AKR7A3", "SLC25A22", "NDUFV2", "NDUFA1", "MYBBP1A", "BID", "QDPR", "DEGS1", "COX7B")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(etc_GO_genes),
                                        name="GO_etc")
pe <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_etc1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score Electron transport chain (ETC) [GO:0022900] ") +
    theme(plot.title = element_text(size = 10)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_etc1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score Electron transport chain (ETC) [GO:0022900] ") +
    theme(plot.title = element_text(size = 10)) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5))

pe
ggsave("./figures/Feature_and_Violinplot_ETC_high_vs_low_salt.pdf",pe, units = "cm", width = 27, height = 15)

#Check each  ETC gene by itself as well
dir.create("./figures/ETC_figs", showWarnings = FALSE)

# Loop through the genes
for(gene in etc_GO_genes) {
  # Create a Seurat object with the single gene as a feature
  single_gene_object <- combined.object.final
  single_gene_object <- AddModuleScore(single_gene_object, 
                                       features = list(gene), 
                                       name = gene)
  
  # Create a violin plot with stat_compare_means and geom_boxplot
  vln_plot <- VlnPlot(single_gene_object, 
                      features = gene, 
                      group.by = "condition",
                      pt.size = 0.1) +
    stat_compare_means(method.args = list(alternative = "two.sided")) +
    geom_boxplot(width = 0.4, size = 0.1) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
    labs(title =gene)
  
  # Save the violin plot
  ggsave(paste0("./figures/ETC_figs/ETC_genelist_", gene, "_expression.pdf"), 
         plot = vln_plot, 
         units = "cm", width = 10, height = 10)
}

# Clean up any open plots
dev.off()


# Check genes CD3G, CD3D, CD3E, TRA, TRB, and ORAI1
GeneList <- c("CD3G", "CD3D", "CD3E","ORAI1", "TRA", "TRB")

for (gene in GeneList) {
  # Check if the gene is present in your data
  if (gene %in% rownames(combined.object.final)) {
    vln_plot <- VlnPlot(combined.object.final, features = gene, group.by = "condition", pt.size = 0.1) +
      stat_compare_means(method.args = list(alternative = "two.sided")) +
      geom_boxplot(width = 0.4, size = 0.1) +
      scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
      labs(title = gene)
    
    # Save the violin plot
    ggsave(paste0("./figures/Violinplot_", gene, "_expression.pdf"), 
           plot = vln_plot, 
           units = "cm", width = 15, height = 15)
  } else {
    cat("Gene", gene, "not found in the data. Skipping...\n")
  }
}

# Subset the data for "ORAI1"
ORAI1_data <- FetchData(combined.object.final, vars = c("ORAI1", "condition"))

# Perform the Wilcoxon test
ORAI1_test_result <- wilcox.test(ORAI1 ~ condition, data = ORAI1_data, alternative = "two.sided")
ORAI1_p_value <- ORAI1_test_result$p.value
ORAI1_p_value # 7.845427e-58


#Cytotoxic geneset from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/cards/BIOCARTA_TCYTOTOXIC_PATHWAY)
cytotoxic_genes <- c("THY1","CD3G","ICAM1","ITGB2","CD3D","CD3E","CD247","ITGAL","CD8A","CD28","PTPRC","CD2")
Idents(object = combined.object.final ) <- "seurat_clusters"
combined.object.final <- AddModuleScore(combined.object.final,
                       features = list(cytotoxic_genes),
                       name="Cytotoxicity")
FeaturePlot(combined.object.final,
            features = "Cytotoxicity1", label = TRUE, repel = TRUE, pt.size = 0.6) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

DotPlot(combined.object.final, features = cytotoxic_genes) + RotatedAxis()


#Check Cytotoxicity (Positive regulation of T cell mediated cytotoxicity: GO:0001916, http://amigo.geneontology.org/amigo/term/GO:0001916)
cytotoxicity_GO_genes <- c("HLA-C","CD1D","CD1E","SLC22A13","HLA-F","MR1","P2RX7","CD1A","ULBP3","XCL1","HFE","ULBP1","ULBP2","HLA-A","HLA-DRB1","RAET1E","B2M","NECTIN2","TAP2","PTPRC","IL12B","IL12A","RAET1L","IL23A","HLA-G","CD1C","CD1B","MICA","MICB","HLA-E","HLA-DRA","HLA-H","HLA-B","PVR","STX7","FADD","RAET1G","CYRIB","IL23R","IL12RB1")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(cytotoxicity_GO_genes),
                                        name="GO_cytotoxicity")

p22 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_cytotoxicity1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score T cell mediated cytotoxicity [GO:0001916] ") +
    theme(plot.title = element_text(size = 10)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_cytotoxicity1", 
             group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score T cell mediated cytotoxicity [GO:0001916] ") +
    theme(plot.title = element_text(size = 10)) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5))

p22
ggsave("./figures/Feature_and_Violinplot_Cytotoxicity_high_vs_low_salt.pdf",p22, units = "cm", width = 27, height = 15)


p22_leiden <- VlnPlot(object = combined.object.final, features = "GO_cytotoxicity1", 
        group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p22_leiden
ggsave("./figures/Violinplot_Cytotoxicity_leiden.pdf",p22_leiden, units = "cm", width = 25, height = 15)

#Get exact p-value by manually running wilcox two-sided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("High NaCl", "Low NaCl"))
cytotoxicity_sub <- combined.object.final_sub@meta.data$GO_cytotoxicity1
cytotoxicity_test_result <- wilcox.test(cytotoxicity_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
cytotoxicity_p_value <- cytotoxicity_test_result$p.value 
cytotoxicity_p_value #1.230387e-35


#Check Activation (Positive regulation of T cell activation: GO:0050870, http://amigo.geneontology.org/amigo/term/GO:0050870)
activation_GO_genes <- c("ZP4","SMARCB1","IL6ST","CD46","RASGRP1","TESPA1","CD274","CD74","THY1","CD160","DNAJA3","NFKBID","TNFSF14","LILRB1","PRKCQ","PIK3CA","SELENOK","RARA","HLX","HLA-DQB2","JAK3","ICOS","ARID1B","RAG1","SLC7A1","CCR7","LILRB4","HSPD1","SMARCD1","BCL10","ABL1","ABL2","DPP4","KLHL25","DHPS","CD28","NCK2","PTPN11","PNP","SOCS5","DOCK8","ICOSLG","MALT1","TMIGD2","CD1D","PCK1","TNFSF13B","GLI3","GLI2","CD86","PTPN22","VNN1","ZBTB1","HLA-DRB3","IHH","HES1","ZMIZ1","IL6","VAV1","HHLA2","PYCARD","CD27","RASAL3","TRAF6","IGFBP2","KITLG","BTN2A2","CD80","CARD11","BCL6","RPS3","NKAP","CCDC88B","CD55","CCL19","SHH","SHB","KLRK1","FYN","LCK","MDK","PHF10","KAT5","CLECL1P","IL7R","SMARCD2","BAD","ZBTB16","PRKCZ","TNFRSF14","AKT1","EFNB3","IGF2","GPAM","ZAP70","ZP3","NFKBIZ","NCKAP1L","WNT10B","IL2RA","EPO","SMARCC1","DUSP10","SYK","CCR2","CD5","IL15","XCL1","EFNB2","ACTL6A","RHOA","HLA-DPB1","HLA-A","PBRM1","HLA-DQA2","HLA-DRB1","HLA-DQB1","SMARCA4","SMARCA2","B2M","MAP3K8","TNFSF9","LGALS9","CSK","HSPH1","NOD2","CD4","IL21","HLA-DRB5","STAT5B","HAVCR2","CORO1A","AGER","SIRPB1","AP3B1","CBFB","IL1B","IL1A","IFNG","FOXP3","LILRB2","CCL21","AIF1","HLA-DOA","IL2RG","RHOH","CD47","HLA-DMB","HLA-DMA","LEP","IL1RL2","RUNX3","PTPRC","IL12B","IL12A","CD70","TGFBR2","RUNX1","CD3E","CD83","PIK3R6","PPP3CA","SPTA1","NLRP3","EFNB1","ADA","RIPK2","SOX4","VSIR","LYN","YES1","VCAM1","CD209","ITPKB","PTPN6","CR1","PDCD1LG2","TYK2","EGR3","BRD7","SASH3","IL23A","XBP1","HLA-G","CD81","GATA3","FCHO1","SLAMF1","ARID1A","JAK2","IL4","HLA-DPA1","SOX13","CD6","FLOT2","IGF1","EBI3","AMBRA1","SMARCE1","HLA-DOB","HLA-DRB4","HLA-E","IL2","TNFSF11","ANXA1","ACTL6B","NCK1","SIRPG","HLA-DQA1","HLA-DRA","SART1","TNFSF4","VTCN1","IL7","IL18","CD24","SMARCC2","BMI1","IL36B","CAV1","IL4I1","CD276","TNFRSF13C","FADD","TFRC","CD40LG","HMGB1","IL4R","ZBTB7B","SRC","SOX12","ARID2","AP3D1","LGALS1","ACTB","SMARCD3","CYRIB","ADAM8","SIRPA","SOCS1","CCL5","CCL2","LEF1","FOXO3","IL23R","IL12RB1")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(activation_GO_genes),
                                        name="GO_activation")

p23 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_activation1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score T cell activation [GO:0050870] ") +
    theme(plot.title = element_text(size = 10)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_activation1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score T cell activation [GO:0050870] ") +
    theme(plot.title = element_text(size = 10)) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p23
ggsave("./figures/Feature_and_Violinplot_Activation_high_vs_low_salt.pdf",p23, units = "cm", width = 27, height = 15)

#Get exact p-value by manually running wilcox two-sided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("High NaCl", "Low NaCl"))
Activation_sub <- combined.object.final_sub@meta.data$GO_activation1
Activation_test_result <- wilcox.test(Activation_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
Activation_p_value <- Activation_test_result$p.value 
Activation_p_value #8.061671e-85

p23_leiden <- VlnPlot(object = combined.object.final, features = "GO_activation1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p23_leiden
ggsave("./figures/Violinplot_Activation_leiden.pdf",p23_leiden, units = "cm", width = 25, height = 15)


#Check Effector function
effctor_GO_genes <- c("GZMB","GZMK","IFNG","EOMES","ITGAD","ITGAX","ITGB7","CXCR3","CCR5")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(effctor_GO_genes),
                                        name="GO_effector_function")

p24 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_effector_function1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score Effector function [Adam M. Drake et al. Science] ") +
    theme(plot.title = element_text(size = 8)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_effector_function1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score Effector function [Adam M. Drake et al. Science]") +
    theme(plot.title = element_text(size = 8)) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p24
ggsave("./figures/Feature_and_Violinplot_Effector_function_high_vs_low_salt.pdf",p24, units = "cm", width = 27, height = 15)

#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("High NaCl", "Low NaCl"))
Effector_sub <- combined.object.final_sub@meta.data$GO_effector_function1
Effector_test_result <- wilcox.test(Effector_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
Effector_p_value <- Effector_test_result$p.value 
Effector_p_value # 1.600367e-61

p24_leiden <- VlnPlot(object = combined.object.final, features = "GO_effector_function1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p24_leiden
ggsave("./figures/Violinplot_Effector_function_leiden.pdf",p24_leiden, units = "cm", width = 25, height = 15)


#Check TCR signalling function (GO:0050862)
TCR_sig_GO_genes <- c("RELA","CCR7","BCL10","CD81","TESPA1","RAB29","CARD11","RPS3","TRAT1","IKBKG","LCK","CD226","PRKD2","NECTIN2","ADA","KCNN4")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(TCR_sig_GO_genes),
                                        name="GO_TCR_signalling")

p25 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_TCR_signalling1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score TCR signalling function [GO:0050862] ") +
    theme(plot.title = element_text(size = 10))  +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_TCR_signalling1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score TCR signalling function [GO:0050862]") +
    theme(plot.title = element_text(size = 10))  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p25
ggsave("./figures/Feature_and_Violinplot_TCR_signalling_high_vs_low_salt.pdf",p25, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("High NaCl", "Low NaCl"))
TCRsig_sub <- combined.object.final_sub@meta.data$GO_TCR_signalling1
TCRsig_test_result <- wilcox.test(TCRsig_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
TCRsig_p_value <- TCRsig_test_result$p.value 
TCRsig_p_value #2.438073e-31

p25_leiden <- VlnPlot(object = combined.object.final, features = "GO_TCR_signalling1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p25_leiden
ggsave("./figures/Violinplot_TCR_signalling_leiden.pdf",p25_leiden, units = "cm", width = 25, height = 15)


#Check stess activated MAPK signalling function (GO:0051403)
MAPK_GO_genes <- c("MYD88","MAPK10","ARHGEF6","LGALS9","ERCC6","SH2D3C","WNT5A","NOD2","MAP4K3","CASR","IRAK1","AGT","MAPK11","MAPK8IP2","LRRK2","MAPK1","ZFP36L1","MINK1","DUSP10","MAP4K1","ZFP36","TNF","CDC42EP5","DUSP9","PAFAH1B1","IL1B","MAP3K5","NOX1","IKBKB","MAP3K12","MAP4K2","MAP2K3","MAPK8IP1","PBK","ATF2","MAP2K4","MAPK9","MAPK8","TNFSF11","MAP2K7","TNFRSF19","TLR4","TLR3","STRADB","IRAK4","MAPK13","STK3","PTGER4","CRYAB","SH2D3A","RIPK2","NPHS1","GPS1","TLR7","SMAD3","MAP3K20","CARD9","MAP3K11","MAPK14","CCM2","MAP2K6","MAP3K7","MAPK3","DAXX","MAPKAPK2","HACD3","NFKB1","GPS2","CRKL","TRIB1","MAP3K10","MAP3K13","NOD1","TAOK2")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(MAPK_GO_genes),
                                        name="GO_Stress_activated_MAPK_cascade")

p26 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_Stress_activated_MAPK_cascade1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score MAPK signalling function [GO:0051403]") +
    theme(plot.title = element_text(size = 10)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_Stress_activated_MAPK_cascade1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score MAPK signalling function [GO:0051403]") +
    theme(plot.title = element_text(size = 10)) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p26
ggsave("./figures/Feature_and_Violinplot_MAPK_cascade_high_vs_low_salt.pdf",p26, units = "cm", width = 27, height = 15)

#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("High NaCl", "Low NaCl"))
MAPK_sub <- combined.object.final_sub@meta.data$GO_Stress_activated_MAPK_cascade1
MAPK_test_result <- wilcox.test(MAPK_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
MAPK_p_value <- MAPK_test_result$p.value 
MAPK_p_value # 1.621521e-118

p26_leiden <- VlnPlot(object = combined.object.final, features = "GO_Stress_activated_MAPK_cascade1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p26_leiden
ggsave("./figures/Violinplot_MAPK_cascade_leiden.pdf",p26_leiden, units = "cm", width = 25, height = 15)


#Check positive regulation of calcineurin-NFAT signaling cascade (GO:0070886)
NFAT_GO_genes <- c("CIB1","PPP3CC","PTBP1","SPPL3","ERBB3","STIMATE","PPP3CB","PPP3R1","TNF","IGF1","CEFIP","PPP3R2","LACRT","PPP3CA","SLC9A1","AKAP6","CAMTA1","AKAP5","CHP2","CHERP","LMCD1")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(NFAT_GO_genes),
                                        name="GO_NFAT_signalling")

p27 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_NFAT_signalling1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score NFAT signaling cascade [GO:0070886]") +
    theme(plot.title = element_text(size = 10))  +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_NFAT_signalling1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score NFAT signaling cascade [GO:0070886]") +
    theme(plot.title = element_text(size = 10))  + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p27
ggsave("./figures/MA_Feature_and_Violinplot_NFAT_signalling_high_vs_low_salt.pdf",p27, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("High NaCl", "Low NaCl"))
NFAT_sub <- combined.object.final_sub@meta.data$GO_NFAT_signalling1
NFAT_test_result <- wilcox.test(NFAT_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
NFAT_p_value <- NFAT_test_result$p.value 
NFAT_p_value #4.498713e-118

p27_leiden <- VlnPlot(object = combined.object.final, features = "GO_NFAT_signalling1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p27_leiden
ggsave("./figures/MA_Violinplot_NFAT_cascade_leiden.pdf",p27_leiden, units = "cm", width = 25, height = 15)


#Check TRM signature genes (#gratz_trm_signature_genes_CD103+CCR7-_skin_vs_CD103-CCR7-_skin_lfc>2_padj<0.05
#Klicznik MM, Morawski PA, Höllbacher B, Varkhande SR et al. Human CD4+CD103+ cutaneous resident memory T cells are found in the circulation of healthy individuals. Sci Immunol 2019 Jul 5;4(37). PMID: 31278120

TRM_GO_genes <- c("XIST","UBC","LGALS3","MT-CO2","VIM","ANKRD28","RGS1","RGCC","HSPA1B","MT-ND4","HSP90AB1","PPP1R15A")
#TRM_GO_genes <- c("IL22", "LONRF2", "NEDD4L", "RGS18", "FAM110B", "NUGGC", "ADM", "LPCAT2", "AK4", "ZNF239")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(TRM_GO_genes),
                                        name="GO_TRM")

p28 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "GO_TRM1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score TRM signature Gratz, Sci Immunol.") +
    theme(plot.title = element_text(size = 10)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "GO_TRM1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score TRM signature Gratz, Sci Immunol.") +
    theme(plot.title = element_text(size = 10)) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p28
ggsave("./figures/MA_Feature_and_Violinplot_TRM_high_vs_low_salt.pdf",p28, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("High NaCl", "Low NaCl"))
TRM_sub <- combined.object.final_sub@meta.data$GO_TRM1
TRM_test_result <- wilcox.test(TRM_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
TRM_p_value <- TRM_test_result$p.value 
TRM_p_value #5.243797e-58

p28_leiden <- VlnPlot(object = combined.object.final, features = "GO_TRM1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p28_leiden
ggsave("./figures/MA_Violinplot_TRM_leiden.pdf",p28_leiden, units = "cm", width = 25, height = 15)


#Check Exhaustion geneset (Kusnadi, Anthony et al., Science immunology.)
Exhaustion_GO_genes <- c("TOX","SNX9","ITGAE","HAVCR2","IRF4","UBE2F","LAG3","STMN1","NDFIP2","ENTPD1","PRDM1","CD63","PDCD1","CD2BP2","FKBP1A","CTLA4","RAB27A","TPI1","CD38","DUSP4","CDCA8","TIGIT","PHLDA1","NCAPG2","VCAM1","ITM2A","CDKN3","CD27","IFI35","SNAP47","ISG15","IGFLR1","STAT3","RAD51","WARS","CCNB1","SYNGR2","BUB1","GBP2","SIRPG","LYST","SEMA4A","BST2","CXCR6","PARK7","ACP5","CCL4L2","HLA-DRA","RGS2","FCRL3","NAB1","OSBPL3","ID3","ICOS","CCR5","FAM3C","GOLIM4","PTPN11","ACP5","CKS2","FUT8","GALM","HLA-DMA")
combined.object.final <- AddModuleScore(combined.object.final,
                                        features = list(Exhaustion_GO_genes),
                                        name="Exhaustion")

p29 <- plot_grid(
  FeaturePlot(object = combined.object.final, features = "Exhaustion1", 
              pt.size = 0.4, cols = c("blue", "red"), 
              combine = TRUE) +  ggtitle("Module Score for Exhaustion Kusnadi, Anthony et al., Sci. immunology.") +
    theme(plot.title = element_text(size = 10)) +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))),
  VlnPlot(object = combined.object.final, features = "Exhaustion1", 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) 
  +  ggtitle("Module Score for Exhaustion Kusnadi, Anthony et al., Sci. immunology.") +
    theme(plot.title = element_text(size = 10)) + stat_compare_means(method.args = list(alternative = "two.sided")) + geom_boxplot(width=0.8, size=0.5) )

p29
ggsave("./figures/MA_Feature_and_Violinplot_Exhaustion_high_vs_low_salt.pdf",p29, units = "cm", width = 27, height = 15)


#Get exact p-value by manually running wilcox twosided test again
combined.object.final_sub <- subset(combined.object.final, idents = c("high salt", "low salt"))
Exhaustion_sub <- combined.object.final_sub@meta.data$Exhaustion1
Exhaustion_test_result <- wilcox.test(Exhaustion_sub ~ combined.object.final_sub$orig.ident, alternative = "two.sided")
Exhaustion_p_value <- Exhaustion_test_result$p.value 
Exhaustion_p_value #7.609517e-76

p29_leiden <- VlnPlot(object = combined.object.final, features = "Exhaustion1", 
                      group.by = "seurat_clusters") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + geom_boxplot(width=0.8, size=0.5)
p29_leiden
ggsave("./figures/MA_Violinplot_Exhaustion_leiden.pdf",p28_leiden, units = "cm", width = 25, height = 15)


#Make violinplot and UMAP for specific genes of interest (ICOS, PDCD1, CTLA4, ITGAE)
features <- c("ICOS", "PDCD1", "CTLA4", "ITGAE")

p30 <- FeaturePlot(object = combined.object.final, features = features) 
p30
ggsave("./figures/MA_Featureplot_genes_of_interest_high_vs_low_salt.pdf",p30, units = "cm", width = 20, height = 15)

#ICOS
p31 <- VlnPlot(object = combined.object.final, features = c("ICOS"), 
          group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p31
ggsave("./figures/MA_Violinplot_ICOS_gene_high_vs_low_salt.pdf",p31, units = "cm", width = 13, height = 13)

#PDCD1
p32 <- VlnPlot(object = combined.object.final, features = c("PDCD1"), 
               group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p32
ggsave("./figures/MA_Violinplot_PDCD1_gene_high_vs_low_salt.pdf",p32, units = "cm", width = 13, height = 13)

#CTLA4
p33 <- VlnPlot(object = combined.object.final, features = c("CTLA4"), 
               group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p33
ggsave("./figures/MA_Violinplot_CTLA4_gene_high_vs_low_salt.pdf",p33, units = "cm", width = 13, height = 13)


#ITGAE
p34 <- VlnPlot(object = combined.object.final, features = c("ITGAE"), 
               group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p34
ggsave("./figures/MA_Violinplot_ITGAE_gene_high_vs_low_salt.pdf",p34, units = "cm", width = 13, height = 13)


#SLC9A1
p35 <- VlnPlot(object = combined.object.final, features = c("SLC9A1"), 
               group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p35
ggsave("./figures/MA_Violinplot_SLC9A1_gene_high_vs_low_salt.pdf",p35, units = "cm", width = 13, height = 13)


#CD28
p36 <- VlnPlot(object = combined.object.final, features = c("CD28"), 
               group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p36
ggsave("./figures/MA_Violinplot_CD28_gene_high_vs_low_salt.pdf",p36, units = "cm", width = 13, height = 13)


#KCNA3
p37 <- VlnPlot(object = combined.object.final, features = c("KCNA3"), 
               group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p37
ggsave("./figures/MA_Violinplot_KCNA3_gene_high_vs_low_salt.pdf",p37, units = "cm", width = 13, height = 13)


#ATP1A1
p38 <- VlnPlot(object = combined.object.final, features = c("ATP1A1"), 
               group.by = "condition") +  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + stat_compare_means(method.args = list(alternative = "two.sided"))
p38
ggsave("./figures/MA_Violinplot_ATP1A1_gene_high_vs_low_salt.pdf",p38, units = "cm", width = 13, height = 13)
