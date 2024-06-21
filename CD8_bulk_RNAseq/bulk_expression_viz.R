#===================================================================================================================================#
# Title: Vizualise top50 DEGs with heatmaps, selected genes in volcano plot
#        and pathview
# Author: Sascha Schäuble
# Figures: Figure 2b, Extended Figure 5b, Extended Figure 13a
#===================================================================================================================================#

# Description: generic
# 
# author Sascha Schäuble
# date of creation: Fri Apr  8 10:48:41 2022
# license: CC BY-SA 4.0

library("tidyverse")
library("magrittr")
library("ComplexHeatmap")
library("circlize")
library("pathview")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

DATE_STR <- format(Sys.time(), "%y%m%d")

FN1 <- "cd8_salt_diff_expression.csv"
FN2 <- "cd8_salt_cd8_highsalt_vs_cd8_lowsalt_diff_expression.csv"

# configs
SAVE_OUT <- F
SORTBYP <- T # T: p-value, F: fold change
PVAL_CUTOFF <- 0.05
FC_CUTOFF   <- 0.5
# ================================================================ #

#### functions #####################################################

#
#### Description: draw volcano plot
#
# IN:
#    df with columns: genename, logFC, adj.P.Val, Status, Labels
#
get_volc <-
  function(df,
           log_fc_cutoff = FC_CUTOFF,
           adj_pval_cutoff = PVAL_CUTOFF,
           font_size = 16,
           padding = 0.5) {
    p <- df %>%
      ggplot(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        col = Status,
        label = Labels
      )) +
      geom_point(size = 4, alpha = 0.2) +
      ggrepel::geom_label_repel(
        alpha = 0.5,
        size = 10,
        label.size = NA,
        color = "black",
        fill = "white",
        fontface = "italic",
        force = 40,
        show.legend = FALSE,
        box.padding = padding,
        max.overlaps = Inf,
        min.segment.length = 0,
        nudge_y = 0,
        nudge_x = 0,
        seed = 17
      ) +
      ggrepel::geom_label_repel(
        alpha = 1,
        size = 10,
        label.size = NA,
        color = "black",
        fill = NA,
        fontface = "italic",
        force = 40,
        show.legend = FALSE,
        box.padding = padding,
        max.overlaps = Inf,
        min.segment.length = 0,
        nudge_y = 0,
        nudge_x = 0,
        seed = 17
      ) +
      geom_vline(
        xintercept = c(-log_fc_cutoff, log_fc_cutoff),
        col = "gray50",
        lty = "dashed"
      ) +
      geom_hline(
        yintercept = -log10(adj_pval_cutoff),
        col = "gray50",
        lty = "dashed"
      ) +
      scale_color_manual(
        breaks = c("up-regulated", "down-regulated", "no significant change"),
        values = c("#ce2627", "#4d69b1", "#989898")
      ) +
      scale_x_continuous(limits = c(-(ceiling(
        max(volc_data$logFC)
      )), ceiling(max(volc_data$logFC))),
      n.breaks = 7) +
      labs(x = expression(log[2] ~ fold ~ change),
           y = expression(-log[10] ~ adjusted ~ p.value)) +
      guides(colour = guide_legend(override.aes = list(size = 6))) +
      theme_bw() +
      theme(
        text = element_text(size = font_size),
        axis.text = element_text(size = font_size),
        axis.title = element_text(size = font_size),
        legend.text = element_text(size = font_size*0.9),
        legend.title = element_blank(),
        legend.position=c(0.25, 0.9)
      ) +
      labs(x = expression(log[2](fold~change)),
           y = expression(-log[10](adjusted~p-value)
           )
      )
    return(p)
    
  }
# ================================================================ #


#
#### Description: provide heatmap
#
# IN:
#    df - expression selection from cd8 dataset
#    dat.deg - df of "cd8_salt_diff_expression.csv"
#    PVALLOG10 - transform p values to -log10 yes or no?
#    SIGNIF
#    SORT - "FC" or "P"
#
get_heatmap_cd8_goi <-
  function(df,
           PVALLOG10 = F,
           TITLE = NULL,
           TITLE_FS = 12,
           SIGNIF = T,
           SORT = "P") {
    require("ComplexHeatmap")
    
    if (SORT == "FC") {
      df %<>% arrange(LFC_cd8_highsalt_vs_cd8_lowsalt)
    } else {
      df %<>% arrange(FDR_cd8_highsalt_vs_cd8_lowsalt)
    }
    
    if (SIGNIF) {
      df %<>% filter(FDR_cd8_highsalt_vs_cd8_lowsalt <= 0.05)
    }
    
    ### prepare data
    ## convert to mat:
    df.mat <- df %>%
      dplyr::select(starts_with("cd8")) %>%
      relocate(starts_with("cd8_high"), .before = "cd8_lowsalt_1") %>% as.matrix()

    rownames(df.mat) <- df %>% pull(genename)
    ## mean expr anno
    df.base_mean <- apply(df.mat, 1, mean)
    df.base_mean %<>% log10()
    
    ## pval info
    anno.pval.vec <- df$FDR_cd8_highsalt_vs_cd8_lowsalt
    names(anno.pval.vec) <- df$genename
    ## lfc anno
    anno.lfc.vec <- df$LFC_cd8_highsalt_vs_cd8_lowsalt
    names(anno.lfc.vec) <- df$genename
    
    df.mat.scaled <- apply(df.mat, 1, scale) %>% t()
    colnames(df.mat.scaled) <- df.mat %>% colnames()
    
    ht.core <- ComplexHeatmap::Heatmap(
      df.mat.scaled,
      col = colorRamp2(breaks = c(-2, 0, 2),
                       hcl.colors(3, "Blue-Red")),
      column_title = TITLE,
      column_title_gp = gpar(fontsize = TITLE_FS),
      cluster_rows = F,
      cluster_columns = F,
      row_names_side = "left",
      heatmap_legend_param = list(title = "Expression (z-scaled)"),
      width = 5,
      row_names_gp = gpar(fontsize = 7),
    )
    
    ht.mean <- Heatmap(
      df.base_mean,
      name = "Base expr.",
      col =
        colorRamp2(breaks = seq(0, 4,
                                length.out = 100),
                   rev(hcl.colors(100, "Viridis"))),
      width = unit(5, "mm"),
      heatmap_legend_param = list(title = "Base expression\n(log10(mean))"),
      show_row_names = F
    )
    
    if (PVALLOG10) {
      # true
      ht.pval <- Heatmap(
        -log10(anno.pval.vec),
        name = "Adj. p-value",
        col =
          colorRamp2(breaks = seq(0, 100,
                                  length.out = 100),
                     rev(hcl.colors(100, "Plasma"))),
        width = unit(5, "mm"),
        heatmap_legend_param = list(title = "Adj. p-value (-log10))",
                                    at = c(0, 20, 40, 60, 80, 100)),
        show_row_names = F
      )
    } else {
      # false
      ht.pval <- Heatmap((anno.pval.vec),
                         name = "Adj. p-value",
                         col =
                           colorRamp2(breaks = seq(0.05, 0,
                                                   length.out = 100),
                                      rev(hcl.colors(100, "Plasma"))),
                         width = unit(5, "mm"),
                         heatmap_legend_param = list(title = "Adj. p-value",
                                                     at = c(0, 0.01, 0.05)),
                         show_row_names = F
      )
    }
    
    ht.lfc <-
      Heatmap(
        anno.lfc.vec,
        name = "Fold change",
        col =
          colorRamp2(breaks = seq(
            min(0,floor(min(anno.lfc.vec, na.rm = T))),
            ceiling(max(anno.lfc.vec, na.rm = T)),
                     length.out=100),
          rev(hcl.colors(100, "Rocket"))),
        width = unit(5, "mm"),
        heatmap_legend_param = list(title = "Fold change (log2)"),
        show_row_names = F
      )
    
    ht <- draw(ht.core + ht.mean + ht.pval + ht.lfc)
    return(ht)
  }


# ================================================================ #

# ================================================================ #



#### Data wrangling ################################################
dat.deg <- read_csv(FN1)
if (SORTBYP) {
  dat.deg %<>% arrange(FDR_cd8_highsalt_vs_cd8_lowsalt)
} else {
  dat.deg %<>% arrange(desc(abs(LFC_cd8_highsalt_vs_cd8_lowsalt)))
}

degs.top <- list()
degs.top$posTop50 <- dat.deg %>% filter(DE_cd8_highsalt_vs_cd8_lowsalt == "upregulated") %>% pull(geneid) %>% head(50)
degs.top$negTop50 <- dat.deg %>% filter(DE_cd8_highsalt_vs_cd8_lowsalt == "downregulated") %>% pull(geneid) %>% head(50)

dat.expr <- read_csv(FN2)
colnames(dat.expr)[1] <- "ensembl_ID"
dat.expr %<>% relocate(starts_with("cd8_high"), .after = "cd8_lowsalt_3")
dat.expr.signif <- dat.expr[match(dat.deg %>%
                              filter(FDR_cd8_highsalt_vs_cd8_lowsalt <= 0.05) %>%
                              pull(geneid), dat.expr$ensembl_ID),]

dat.expr.signif %<>% left_join(
  dat.deg %>% dplyr::select(
    geneid,
    genename,
    description,
    LFC_cd8_highsalt_vs_cd8_lowsalt,
    FDR_cd8_highsalt_vs_cd8_lowsalt
  ),
  by = c("ensembl_ID" = "geneid")
)

dat.expr.top <- list()
dat.expr.top$up50 <- dat.expr.signif %>% filter(ensembl_ID %in% degs.top$posTop50 )
dat.expr.top$dn50 <- dat.expr.signif %>% filter(ensembl_ID %in% degs.top$negTop50 )

# ================================================================ #

#### Volcano plot ##################################################
volc_data <- dat.deg %>%
  dplyr::select(genename, LFC_cd8_highsalt_vs_cd8_lowsalt, FDR_cd8_highsalt_vs_cd8_lowsalt) %>%
  dplyr::rename("logFC" = "LFC_cd8_highsalt_vs_cd8_lowsalt", "adj.P.Val" = "FDR_cd8_highsalt_vs_cd8_lowsalt") %>%
  mutate(DEG = case_when(
    abs(logFC) >= FC_CUTOFF & adj.P.Val < PVAL_CUTOFF ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(Status = case_when(
    DEG == TRUE & logFC >= FC_CUTOFF ~ "up-regulated",
    DEG == TRUE & logFC <= -FC_CUTOFF ~ "down-regulated",
    TRUE ~ "no significant change"
  ))

sel_genes <- c(
  "SLC5A3",
  "SLC35F3",
  "SLC12A8",
  "SLC29A1",
  "SGK1",
  "BATF3",
  "SLC7A5",
  "IRF4",
  "IL7R",
  "LTA",
  "HK1",
  "HK2",
  "MYC", 
  "HIF1A", 
  "CD24"
)
volc_data %<>%
  mutate(Labels = case_when(
    genename %in% sel_genes ~ genename,
    TRUE ~ NA_character_
  ))
volc_data %>% filter(!is.na(Labels))

### draw volcano plots 
p.volc <- get_volc(volc_data, padding = 0.35, font_size = 24)
if (SAVE_OUT) {
  p.volc %>%
    ggsave(filename = paste0("deg_volc", "_", DATE_STR, ".pdf"), width = 9, height = 9)
}
# ================================================================ #

#### heatmaps top50 ################################################

### up 50
p.degs.up50 <- get_heatmap_cd8_goi(
  df = dat.expr.top$up50,
  PVALLOG10 = T,
)
if (SAVE_OUT) {
  cairo_pdf(
    filename = paste0("degs_expr_topUp50_signif", "_",DATE_STR, ".pdf"),
    width = 6,
    height = 8
  )
  draw(p.degs.up50)
  dev.off()
}

### dn50
p.degs.dn50 <- get_heatmap_cd8_goi(
  df = dat.expr.top$dn50,
  PVALLOG10 = T
)
if (SAVE_OUT) {
  cairo_pdf(
    filename = paste0("degs_expr_topDown50_signif", "_",DATE_STR, ".pdf"),
    width = 6,
    height = 8
  )
  draw(p.degs.dn50)
  dev.off()
}
# ================================================================ #

#### pathview glycolysis ###########################################
dat.4pv <- dat.deg$LFC_cd8_highsalt_vs_cd8_lowsalt %>% as.matrix()
dimnames(dat.4pv) = list(dat.deg$geneid,  c("log2FC"))
pathview(gene.data =  dat.4pv
           , pathway.id = "00010" # Glycolysis / Gluconeogenesis
           , species = "hsa"
           , gene.idtype = "ENSEMBL"
           , limit = list(gene = 3)
           , out.suffix = paste0("pv_CD8_highVSlow_glycolysis_", DATE_STR)
           , multi.state = F
           , low = (gene = hcl.colors(n = 9, palette = "Blue-Red 2")[2])
           , mid = (gene = hcl.colors(n = 9, palette = "Blue-Red 2")[5])
           , high = (gene = hcl.colors(n = 9, palette = "Blue-Red 2")[8])
)
# ================================================================ #
