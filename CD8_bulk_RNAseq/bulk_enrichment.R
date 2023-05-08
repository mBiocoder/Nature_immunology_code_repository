#===================================================================================================================================#
# Title: Perform enrichment analysis with bulk RNA seq data
# Author: Sascha Schäuble
# Figures: Figure 4A,H, Supplementary Figure S12B 
#===================================================================================================================================#


# Description: 
#
# author Sascha Schäuble
# date of creation: Fri Jun 24 11:21:11 2022
# license: CC BY-SA 4.0

library("tidyverse")
library("magrittr")
library("data.table")
library("ggpubr")
library("ggsci")
library("clusterProfiler")
library("fgsea")
library("limma")
library("org.Hs.eg.db")





DATE_STR <- format(Sys.time(), "%y%m%d")

FN1 <- "cd8_salt_cd8_highsalt_vs_cd8_lowsalt_diff_expression.csv"
FN2 <- "cd8_salt_diff_expression.csv"
FN3 <- "revigo_ora_upGOBP_mitos_230418.tsv" # REVIGO reduction using 0.5 with 

# configs
SAVE_OUT <- F
ADD_EPS2DUPSEXPR <- T # some log2fc have the same value - add a small eps to resolve
REMOVE_DUPIDS <- T
PVAL_CUT <- 0.05
# ================================================================ #


#### functions #####################################################
#
#### Description: adds a small epsilon value to e.g. fold change 
#                 data in case these are identical (which is bad for
#                 GSEA barcode plots)

#
# IN:
#    df - the data frame with data from the high/low salt project
#    my_seed - a specific seed value
# OUT:
#    changed df
# PRE:
#   df holds column "LFC_cd8_highsalt_vs_cd8_lowsalt"
#
add_eps2_dups <- function(df, my_seed = 1722) {
  dummy.dupsIndex <-
    rank(df$LFC_cd8_highsalt_vs_cd8_lowsalt) %>%  duplicated()
  set.seed(my_seed)
  df$LFC_cd8_highsalt_vs_cd8_lowsalt[dummy.dupsIndex] <-
    df$LFC_cd8_highsalt_vs_cd8_lowsalt[dummy.dupsIndex] +
    sample(seq(1e-11, 1e-10, length.out = sum(dummy.dupsIndex)),
           sum(dummy.dupsIndex),
           replace = F)
  df %<>%
    arrange(desc(LFC_cd8_highsalt_vs_cd8_lowsalt))
  
  return(df)
}
# ================================================================ #


#### Description: gives back per factor the top DEGs
#
# IN:
#    df with data and comp variable for comparisons
#    criteria - which column should be arranged and the criteria?
#    direction - asc or desc?
#    topN - how many top N genes to pick (or fewer if there are less)?
#    id - which gene id? default is "geneid" which is ensembl
# OUT:
#   df with top degs, unique rows for comp and geneid
# PRE:
#   df has to have "comp" and id variable
#
pickTopGenes <- function(df, criteria, direction = "desc", topN = NULL, id = "geneid") {
  res <- tibble()
  for (c in levels(df$comp)) {
    tmp <- df %>% filter(comp == c) %>% 
      dplyr::select(comp,id) %>% 
      distinct()
    
    if (direction == "desc")
      tmp %<>% arrange(desc(criteria))
    else
      tmp %<>% arrange(criteria)
    
    if (!is.null(topN))
      tmp %<>% head(topN)
    
    res %<>% rbind(tmp)
    
  }
  return(res)
}
# ================================================================ #

#
#### Description: simple helper
#
# IN:
#    df - salt
#    id - "symbol" - SYMBOL, ensembl for any other string
# PRE:
#   genename or geneid column in df
#
getTopSaltGenes <- function(df, term = "salt_signature", n = 100,
                            id = "symbol", returnType = "list") {
  if (id == "symbol") {
    res <- tibble(signature = term,
                  id = df$genename %>% head(n))
  } else {
    res <- tibble(signature = term,
                  id = df$geneid %>% head(n))
  }
  if (returnType == "list") {
    res.list <- list()
    for (f in (res$signature %>% as_factor() %>% levels())) {
      res.list[[f]] <- res %>% filter(signature == f) %>% pull(id)
    }
    res <- res.list
  }
  return(res)
}
# ================================================================ #

#
#### Description: produces barcode plots adapted from fgsea package
#
getMyBarcode <-
  function(pathway,
           stats,
           gseaParam = 1,
           ticksSize = 0.2,
           scoreType = "std",
           minGS = 1,
           maxGS = length(stats) - 1,
           exp = 1,
           nPerm = 500,
           n1 = NA,
           n2 = NA,
           title = NULL) {
    
    # stats
    tmp.fgsea <- fgsea::fgsea(
      pathways = pathway,
      stats = stats,
      minSize = minGS,
      maxSize = maxGS,
      scoreType = scoreType,
      gseaParam = exp,
      nproc = 4,
      eps = 0
    )
    
    # plot - enabling "positive" one tailed test
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))
    pathway <-
      unname(as.vector(na.omit(match(
        pathway[[1]], names(statsAdj)
      ))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(
      statsAdj,
      selectedStats = pathway,
      returnAllExtremes = TRUE,
      scoreType = "pos"
    )
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms)) / 8
    x = y = NULL
    g <- ggplot(toPlot, aes(x = x, y = y)) +
      geom_hline(
        yintercept = max(tops),
        colour = "red",
        linetype = "dashed"
      ) + geom_hline(
        yintercept = min(bottoms),
        colour = "red",
        linetype = "dashed"
      ) +
      geom_hline(yintercept = 0, colour = "black") + geom_line(color = "green") +
      theme_pubr() +
      ggtitle(title) +
      geom_segment(
        data = data.frame(x = pathway),
        mapping = aes(
          x = x,
          y = -diff /
            2,
          xend = x,
          yend = diff / 2
        ),
        size = ticksSize
      ) +
      theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "Rank of genes", y = "Running enrichment score")  +
      geom_point(color = "green", size = 0.1) +
      annotate(
        geom = "text",
        x = 0,
        y = -0.05,
        label = paste0(
          "NES = ",
          format(tmp.fgsea$NES, digits = 3),
          ", p = ",
          format(tmp.fgsea$pval, digits = 3)
        ),
        color = "black",
        size = 4,
        hjust = "left"
      )
    if (!(is.na(n1) | is.na(n2))) {
      g <- g +
        annotate(
          geom = "text",
          x = 0,
          y = -0.1,
          label = paste0("list(n[1]==", n1, ", n[2]==", n2, ")"),
          color = "black",
          size = 4,
          hjust = "left",
          parse = T
        )
    }
    
    return(g)
  }
# ================================================================ #

# ================================================================ #



#### Data wrangling ################################################
### gene gene id map ensembl to entrez
x <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(x)
xx <- x[mapped_genes] %>% as.list()
geneID.map <-
  tibble(entrez = names(xx), ensembl = xx) %>% unnest(cols = ensembl)

### CD8 expression data
dat.expr <- read_csv(FN1)
colnames(dat.expr)[1] <- "ensemblID"
dat.expr.stats <- read_csv(FN2)

dat.expr.stats$DE_cd8_highsalt_vs_cd8_lowsalt %<>% as.factor()
dat.expr.stats %<>% left_join(geneID.map, by = c("geneid" = "ensembl")) %>% relocate(entrez, .after = geneid)
dat.expr.geneList.df <-
  dat.expr.stats %>% dplyr::select(entrez, LFC_cd8_highsalt_vs_cd8_lowsalt) %>%
  filter(!(is.na(LFC_cd8_highsalt_vs_cd8_lowsalt) |
             is.na(entrez))) %>%
  arrange(desc(LFC_cd8_highsalt_vs_cd8_lowsalt))
dat.expr.geneList.ensembl.df <-
  dat.expr.stats %>% dplyr::select(geneid, LFC_cd8_highsalt_vs_cd8_lowsalt) %>%
  filter(!(is.na(LFC_cd8_highsalt_vs_cd8_lowsalt) |
             is.na(geneid))) %>%
  arrange(desc(LFC_cd8_highsalt_vs_cd8_lowsalt))
# add a small epsilon to prevent duplicate fold changes
if (ADD_EPS2DUPSEXPR) {
  dat.expr.geneList.df <- add_eps2_dups(dat.expr.geneList.df)
  dat.expr.geneList.ensembl.df <-
    add_eps2_dups(dat.expr.geneList.ensembl.df)
}
# remove duplicate ids
if (REMOVE_DUPIDS) {
  dat.expr.geneList.df <-
    dat.expr.geneList.df[!(dat.expr.geneList.df$entrez %>% duplicated()), ]
  dat.expr.geneList.ensembl.df <-
    dat.expr.geneList.ensembl.df[!(dat.expr.geneList.ensembl.df$geneid %>%
                                     duplicated()), ]
}
dat.expr.geneList <- dat.expr.geneList.df %>% deframe()
dat.expr.geneList.ensembl <-
  dat.expr.geneList.ensembl.df %>% deframe()

dat.DEG.up <-
  dat.expr.stats %>% filter(DE_cd8_highsalt_vs_cd8_lowsalt == "upregulated") %>%
  dplyr::select(
    geneid,
    genename,
    meanExpression,
    LFC_cd8_highsalt_vs_cd8_lowsalt,
    FDR_cd8_highsalt_vs_cd8_lowsalt
  ) %>% arrange(desc(LFC_cd8_highsalt_vs_cd8_lowsalt))
dat.DEG.down <-
  dat.expr.stats %>% filter(DE_cd8_highsalt_vs_cd8_lowsalt == "downregulated") %>%
  dplyr::select(
    geneid,
    genename,
    meanExpression,
    LFC_cd8_highsalt_vs_cd8_lowsalt,
    FDR_cd8_highsalt_vs_cd8_lowsalt
  ) %>% arrange(desc(LFC_cd8_highsalt_vs_cd8_lowsalt))

# get genes for KEGG pathways
dat.kegg <- limma::getGeneKEGGLinks(species="hsa")
dat.kegg$Symbol <- mapIds(org.Hs.eg.db, dat.kegg$GeneID,
                          column="SYMBOL", keytype="ENTREZID")
dat.kegg %<>% filter(PathwayID == "hsa00010") %>% 
  dplyr::rename(comp = PathwayID)
dat.kegg$comp %<>% as_factor()
dat.kegg %<>% left_join(dat.expr.stats %>% 
                          dplyr::select(entrez,geneid, genename,
                                        FDR_cd8_highsalt_vs_cd8_lowsalt,
                                        LFC_cd8_highsalt_vs_cd8_lowsalt), 
                        by = c("Symbol" = "genename") )

# ====================================================== #


#### run ORA and GSEA ####################################
dat.cd8.4ora <- dat.expr.stats %>%
  filter(DE_cd8_highsalt_vs_cd8_lowsalt != "unchanged") %>%
  dplyr::select(entrez, LFC_cd8_highsalt_vs_cd8_lowsalt) %>%
  filter(!is.na(entrez)) %>%
  deframe() %>% sort(decreasing = T) 

gsea <- list()
gsea[["GOBP"]]  <- gseGO( geneList = dat.cd8.4ora, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gsea[["GOMF"]]  <- gseGO( geneList = dat.cd8.4ora, ont = "MF", OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gsea[["KEGG"]]  <- gseKEGG( geneList = dat.cd8.4ora, keyType = "ncbi-geneid")
gsea[["MKEGG"]] <- gseMKEGG(geneList = dat.cd8.4ora, keyType = "ncbi-geneid")

gsea.tibble <- list()
gsea.tibble[["GOBP"]] <- gsea[["GOBP"]]@result %>% as_tibble() 
gsea.tibble[["GOMF"]] <- gsea[["GOMF"]]@result %>% as_tibble() 
gsea.tibble[["KEGG"]] <- gsea[["KEGG"]]@result %>% as_tibble()
gsea.tibble[["MKEGG"]] <- gsea[["MKEGG"]]@result %>% as_tibble() 
gsea.tibble %<>% data.table::rbindlist(idcol = "category")
gsea.tibble$category %>% table()

### ORA
## all genes
ora.all <- list()
ora.all[["GOBP"]]  <- enrichGO(gene = names(dat.cd8.4ora), OrgDb="org.Hs.eg.db", ont = "BP")
ora.all[["GOMF"]]  <- enrichGO(gene = names(dat.cd8.4ora), OrgDb="org.Hs.eg.db", ont = "MF")
ora.all[["KEGG"]]  <- enrichKEGG(gene = names(dat.cd8.4ora), keyType = "ncbi-geneid")
ora.all[["MKEGG"]] <- enrichMKEGG(gene = names(dat.cd8.4ora), keyType = "ncbi-geneid")

## up-regulation
dat.expr.stats.up <- dat.expr.stats %>%
  filter(DE_cd8_highsalt_vs_cd8_lowsalt == "upregulated") %>%
  dplyr::select(entrez, LFC_cd8_highsalt_vs_cd8_lowsalt) %>% deframe() %>% sort(decreasing = T)
dat.expr.stats.down <- dat.expr.stats %>%
  filter(DE_cd8_highsalt_vs_cd8_lowsalt == "downregulated") %>%
  dplyr::select(entrez, LFC_cd8_highsalt_vs_cd8_lowsalt) %>% deframe() %>% sort(decreasing = T)

ora.all.up <- list()
ora.all.up[["GOBP"]]  <- enrichGO(gene = names(dat.expr.stats.up), OrgDb="org.Hs.eg.db", ont = "BP")
ora.all.up[["GOMF"]]  <- enrichGO(gene = names(dat.expr.stats.up), OrgDb="org.Hs.eg.db", ont = "MF")
ora.all.up[["KEGG"]]  <- enrichKEGG(gene = names(dat.expr.stats.up), keyType = "ncbi-geneid")
ora.all.up[["MKEGG"]] <- enrichMKEGG(gene = names(dat.expr.stats.up), keyType = "ncbi-geneid")

# down-regulation
ora.all.down <- list()
ora.all.down[["GOBP"]]  <- enrichGO(gene = names(dat.expr.stats.down), OrgDb="org.Hs.eg.db", ont = "BP")
ora.all.down[["GOMF"]]  <- enrichGO(gene = names(dat.expr.stats.down), OrgDb="org.Hs.eg.db", ont = "MF")
ora.all.down[["KEGG"]]  <- enrichKEGG(gene = names(dat.expr.stats.down), organism = "hsa",
                                      keyType = "ncbi-geneid")
ora.all.down[["MKEGG"]] <- enrichMKEGG(gene = names(dat.expr.stats.down), keyType = "ncbi-geneid")

ora.tibble <- list()
ora.tibble[["all_GOBP"]] <-   ora.all[["GOBP"]]@result %>% as_tibble() 
ora.tibble[["all_GOMF"]] <-   ora.all[["GOMF"]]@result %>% as_tibble() 
ora.tibble[["all_KEGG"]] <-   ora.all[["KEGG"]]@result %>% as_tibble()
ora.tibble[["all_MKEGG"]] <-  ora.all[["MKEGG"]]@result %>% as_tibble() 
ora.tibble[["up_GOBP"]] <-    ora.all.up[["GOBP"]]@result %>% as_tibble() 
ora.tibble[["up_GOMF"]] <-    ora.all.up[["GOMF"]]@result %>% as_tibble() 
ora.tibble[["up_KEGG"]] <-    ora.all.up[["KEGG"]]@result %>% as_tibble()
ora.tibble[["up_MKEGG"]] <-   ora.all.up[["MKEGG"]]@result %>% as_tibble() 
ora.tibble[["down_GOBP"]] <-  ora.all.down[["GOBP"]]@result %>% as_tibble() 
ora.tibble[["down_GOMF"]] <-  ora.all.down[["GOMF"]]@result %>% as_tibble() 
ora.tibble[["down_KEGG"]] <-  ora.all.down[["KEGG"]]@result %>% as_tibble()
ora.tibble[["down_MKEGG"]] <- ora.all.down[["MKEGG"]]@result %>% as_tibble() 
ora.tibble %<>% data.table::rbindlist(idcol = "category")
ora.tibble %<>% filter(qvalue <= 0.05)
ora.tibble$category %>% table()

if (SAVE_OUT) {
  write_tsv(gsea.tibble, paste0(RES_PATH, "gsea_all", "_", DATE_STR, ".tsv"))
  write_tsv(ora.tibble, paste0(RES_PATH, "ora_all", "_", DATE_STR, ".tsv"))
}
# ================================================================ #


#### barplot | highlighted metabolism ##############################
dat.clPrGSEA.up_KEGG.signif <-
  ora.tibble[["up_KEGG"]] %>% filter(qvalue <= 0.05)
kegg.mb <-
  c(
    "hsa00010",
    "hsa00051",
    "hsa00100",
    "hsa00240",
    "hsa00270",
    "hsa00630",
    "hsa00670",
    "hsa00900",
    "hsa01040",
    "hsa01200",
    "hsa01212",
    "hsa01230",
    "hsa01232",
    "hsa01240",
    "hsa00620"
  )
dat.clPrGSEA.up_KEGG.signif %<>% mutate(KEGG = if_else(ID %in% kegg.mb, "Metabolism", "other"))

p.keggEnrich <-
  dat.clPrGSEA.up_KEGG.signif %>% mutate(qvalue = -log10(qvalue)) %>%
  arrange(qvalue) %>%
  ggpubr::ggbarplot(
    x = "Description",
    y = "qvalue",
    fill = "KEGG",
    xlab = "",
    ylab = "-log10(q-value)"
  ) %>%
  ggpar(orientation = "horizontal")

if (SAVE_OUT) {
  ggexport(
    p.keggEnrich,
    filename = paste0("barplotKEGG4upDEGs", "_",
                      DATE_STR, ".pdf"),
    width = 10,
    height = 9
  )
}
# ================================================================ #


#### barplot | mitochondrial categories ############################
dat.enrich.revigo <- read_tsv(FN3) %>%
  filter(Representative == "null") %>% pull(Name)
dat.mito <- ora.tibble %>% 
  filter(category == "up_GOBP") %>% 
  filter(grepl("mitoch",Description)) %>% 
  mutate(qvalMinusLog10 = -log10(qvalue))
p.mito <- ggbarplot(dat.mito %>% arrange(qvalMinusLog10) %>% 
                    filter(Description %in% dat.enrich.revigo),
                    x = "Description", y = "qvalMinusLog10", 
                    orientation = "horizontal", 
                    fill = "gray",
                    xlab = "",
                    ylab = ""
)
p.mito <- p.mito + labs(y = expression(paste("-log"[10], "(q-value)")))

if (SAVE_OUT) {
  p.mito %>% ggexport(
    filename = paste0("barplot_mito_", DATE_STR, ".pdf"),
    width = 8,
    height = 4
  )
}

# ================================================================ #


#### barcode | glycolysis, GSEA ####################################

### KEGG
# using fgsea
dat.glycolysis <- pickTopGenes(df = dat.kegg %>%
                                 filter(comp =="hsa00010" &
                                          FDR_cd8_highsalt_vs_cd8_lowsalt <= PVAL_CUT) %>% 
                                 distinct(),
                               criteria = "LFC_cd8_highsalt_vs_cd8_lowsalt")
n = (dat.glycolysis %>% dim())[1]
p.glycolysis <- getMyBarcode(pathway = getTopSaltGenes(df = dat.glycolysis, n = n, id = "geneid"),
                             stats = ( dat.expr.geneList.ensembl ),
                             scoreType = "pos",
                             maxGS = 4000,
                             minGS = 10, nPerm = 1e4,
                             title = "KEGG glycolysis (ID hsa00010)"
                             
)

if (SAVE_OUT) {
  p.glycolysis %>%
    cowplot::save_plot(
      filename = paste0("KEGGpwys/",
                        "glycolysis_n",
                        n,
                        "_barcode_",
                        DATE_STR,
                        ".pdf"),
      base_height = 6
    )
}