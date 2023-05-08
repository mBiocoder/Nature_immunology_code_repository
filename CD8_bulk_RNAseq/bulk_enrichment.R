#===================================================================================================================================#
# Title: Perform enrichment analysis with bulk RNA seq data
# Author: Sascha Schäuble
# Figures: Figure 4A,H
#===================================================================================================================================#


# Description: 
#
# author Sascha Schäuble
# date of creation: Fri Jun 24 11:21:11 2022
# license: CC BY-SA 4.0

library("tidyverse")
library("magrittr")
library("ggpubr")
library("ggsci")
library("clusterProfiler")
library("org.Hs.eg.db")
library("data.table")


DATE_STR <- format(Sys.time(), "%y%m%d")

FN1 <- "cd8_salt_cd8_highsalt_vs_cd8_lowsalt_diff_expression.csv"
FN2 <- "cd8_salt_diff_expression.csv"
FN3 <- "revigo_ora_upGOBP_mitos_230418.tsv" # REVIGO reduction using 0.5 with 

# configs
SAVE_OUT <- F
ADD_EPS2DUPSEXPR <- T # some log2fc have the same value - add a small eps to resolve
REMOVE_DUPIDS <- T

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

p.mito %>% ggexport(filename = paste0("barplot_mito_", DATE_STR, ".pdf"),
                    width = 8, height = 4) 


# ================================================================ #
