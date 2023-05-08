#===================================================================================================================================#
# Title: Compare bulk RNAseq CD8+ CD45RA- T cells against public TCGA-GTEx data
# Author: Sascha Schäuble
# Figures: Figure 1D, Supplementary Figure S1
#===================================================================================================================================#

#### config ########################################################

library("tidyverse")
library("magrittr")
library("ggpubr")
library("janitor")
library("data.table")
library("fgsea")
library("cowplot")

DATE_STR <- format(Sys.time(), "%y%m%d")

SAVE_OUT <- F

FN1 <- "gencode.v23.annotation.gene.probemap"
FN2 <- "TcgaTargetGtex_rsem_gene_tpm.gz"
FN3 <- "TcgaTargetGTEX_phenotype.txt.gz"
FN4 <- "cd8_salt_diff_expression.csv"

source("src/funs.R")
# ================================================================ #


#### functions #####################################################

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



#
#### Description: simple helper function
#
# OUT:
#   sorted (named) gene vector 
#
get_scaled_deframe <- function(df = expr.all.mean.df, description) {
  
  geneList.sorted <- df %>% filter(dataDescr == description) %>%
    dplyr::select(gene, meanExprCount) %>%
    arrange(desc(meanExprCount)) %>%
    deframe() %>%
    scale(scale = T, center = T) %>% dplyr::as_tibble(rownames = "symbol") %>%
    deframe()
  
  return(geneList.sorted)
}
# ================================================================ #


#
#### Description: simple helper to select top DEGs
#
# IN:
#    df - salt expression data
#    id - "symbol" - SYMBOL, ensembl for any other string
# OUT:
#   variables
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


#### Description: simple helper function based on de-log2 data
# PRE:
#   log2fc is computed by cond1 / cond2; each gets added a small delta to prevent div/0
#
get_fc <-
  function(df = expr.all.mean.df,
           cond1,
           cond2,
           delta = 0.001,
           reduced = T,
           addMeta = T,
           map = dat.map) {
    
    dat1 <- df %>% filter(dataDescr %in% cond1)
    dat2 <- df %>% filter(dataDescr %in% cond2)
    if (length(cond2) > 1) {
      # need to get the true mean over multiple groups
      # use metaInfo for getting the sums and create new mean with correct weights
      dat2.wide <-
        dat2 %>% select(-sampletype) %>% pivot_wider(names_from = "dataDescr",
                                                     values_from = c("meanExpr", "meanExprCount"))
      nTotal = 0
      for (i in cond2) {
        dat2.wide[, paste0("meanExpr_", i)] <-
          dat2.wide[paste0("meanExpr_", i)] *
          dat.pheno.info %>% filter(descr == i) %>% pull(n)
        dat2.wide[, paste0("meanExpr_", i)] <-
          dat2.wide[paste0("meanExprCount_", i)] *
          dat.pheno.info %>% filter(descr == i) %>% pull(n)
        nTotal = nTotal + dat.pheno.info %>% filter(descr == i) %>% pull(n)
      }
      dat2.wide$meanExpr <- (dat2.wide %>%
                               deplyr::select(paste0("meanExpr_", cond2)) %>%
                               rowSums()) / nTotal
      dat2.wide$meanExprCount <- (dat2.wide %>%
                                    deplyr::select(paste0("meanExprCount_", cond2)) %>%
                                    rowSums()) / nTotal
    }
    dat.all <-
      left_join(dat1,
                dat2,
                by = "sample",
                suffix = c(".dat1", ".dat2"))
    dat.all %<>% mutate(log2fc = log2((dat.all$meanExprCount.dat1 + delta) /
                                        (dat.all$meanExprCount.dat2 + delta)
    )) %>%
      arrange(desc(log2fc))
    
    if (reduced) {
      dat.all %<>% dplyr::select(sample, sampletype.dat1, sampletype.dat2, log2fc)
    }
    if (addMeta) {
      dat.all %<>%
        left_join(map %>% dplyr::select(id, gene), by = c("sample" = "id")) %>%
        relocate(gene, .after = "sample")
    }
    return(dat.all)
  }
# ================================================================ #


#
#### Description: simple helper function to compute/set mean expression
#
get_meanExpr <-
  function(df = dat.expr,
           pheno = dat.pheno,
           site,
           sampletype) {
    
    df.mean <- df %>% select(all_of(c(
      "sample",
      pheno %>%
        filter(primary_site == site &
                 sample_type %in% sampletype) %>%
        pull(sample)
    )))
    df.mean <- tibble(sample = df.mean$sample,
                      meanExpr = df.mean %>%
                        select(-sample) %>% rowMeans())
    
    return(df.mean)
  }
# ================================================================ #

# ================================================================ #


#### Data wrangling ################################################
### data download
## retrieve tcgatargetgtex pheno data
if (!file.exists("TcgaTargetGTEX_phenotype.txt.gz")) {
  system("wget https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz")
} 
dat.pheno <- read_tsv(FN3, col_types = c("cffffff"))
dat.pheno %<>% clean_names()

### filter data 
## reduce to primary tumor, solid tissue normal and normal tissue (GTex) and ignore TARGET
dat.pheno %<>% filter(
  study %in% c("TCGA", "GTEX") &
    sample_type %in% c("Primary Tumor", "Normal Tissue", "Solid Tissue Normal") &
    !is.na(primary_site)
)
dat.pheno$sample_type %<>% fct_drop()
dat.pheno$study %<>% fct_drop()

# check data
dat.pheno %>% summary()
dat.pheno$primary_site %>% levels() %>% sort()
dat.pheno$primary_site %<>% fct_recode("Adrenal Gland" = "Adrenal gland") %>% fct_drop()
dat.pheno$primary_site %<>% fct_recode("Cervix" = "Cervix Uteri") %>% fct_drop()
dat.pheno$primary_site %<>% fct_recode("Thyroid Gland" = "Thyroid") %>% fct_drop()
dat.pheno$primary_site %>% levels() %>% sort()
# filter tumor sites
dat.pheno %<>% filter(!(
  primary_site %in% c(
    "Fallopian Tube",
    "Blood",
    "Blood Vessel",
    "Pituitary",
    "Heart",
    "Nerve",
    "Small Intestine",
    "Spleen",
    "Salivary Gland",
    "Vagina",
    "Adipose Tissue",
    "Muscle"
  )
))



### retrieve tcgatargetgtex expression data
if (!file.exists("TcgaTargetGTEX_phenotype.txt.gz")) {
  system("wget https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz")
}
dat.map <- read_tsv(FN1)
dat.expr <- data.table::fread(FN2)
stopifnot(dplyr::setequal(dat.map$id, dat.expr$sample))
# filter dat.expr for valid tcga ids
dat.expr %<>% select(all_of(c("sample", dat.pheno$sample)) )

### salt data
dat.salt.upDEGs <- read_csv(FN4) %>%
  filter(DE_cd8_highsalt_vs_cd8_lowsalt == "upregulated") %>%
  dplyr::select(
    geneid,
    genename,
    meanExpression,
    LFC_cd8_highsalt_vs_cd8_lowsalt,
    FDR_cd8_highsalt_vs_cd8_lowsalt
  ) %>%
  arrange(desc(LFC_cd8_highsalt_vs_cd8_lowsalt))

n_salt_up <- dim(dat.salt.upDEGs)[1]
# ================================================================ #


#### get means of expression #######################################
expr.all.mean <- list()
for (i in 1:dim(dat.pheno.info)[1]) {
  dummy <-
    paste(dat.pheno.info$primary_site[i], 
          dat.pheno.info$sample_type[i], sep = "_") %>%
    str_replace_all(" ", "_")
  
  expr.all.mean[[dummy]] <-
    get_meanExpr(site = dat.pheno.info$primary_site[i],
                 sampletype = dat.pheno.info$sample_type[i])
}
# add "healthy" - combination of "solid normal" and "normal" 
for (i in 1:length(dat.pheno.info$primary_site %>% levels())) {
  site.tmp <- (dat.pheno.info$primary_site %>% levels())[i]
  dummy <-
    paste(site.tmp, "normalAll", sep = "_")
  dummyStype <- dat.pheno.info %>% 
    filter(primary_site == site.tmp &
             sample_type != "Primary Tumor") %>% pull(sample_type) %>% as.character()
  if (length(dummyStype) > 0) {
    expr.all.mean[[dummy]] <-
      get_meanExpr(site = site.tmp,
                   sampletype = dummyStype)
  }
}

# as we do not use it anymore, we can get rid of dat.expr
rm(dat.expr)

# put together the mean expression table
expr.all.mean.df <- rbindlist(expr.all.mean, idcol = "dataDescr")
expr.all.mean.df$dataDescr %<>% str_replace_all(" ", "_") 
expr.all.mean.df$dataDescr %<>% as_factor %>% fct_drop() 

expr.all.mean.df %<>% mutate(site = dataDescr %>% 
                               str_remove("_{1}Primary_Tumor|Solid_Tissue_Normal|Normal_Tissue|normalAll") %>% 
                               as_factor(),
                             sampletype = dataDescr %>% 
                               str_extract("Primary_Tumor|Solid_Tissue_Normal|Normal_Tissue|normalAll") %>% 
                               as_factor())
expr.all.mean.df %<>% mutate(meanExprCount = 2^meanExpr-0.001)
expr.all.mean.df$dataDescr %<>% as_factor()
expr.all.mean.df$site %<>% str_remove("_$") %>% fct_drop()
expr.all.mean.df$site %<>% fct_recode(Adrenal_Gland = "Adrenal_gland") %>% 
  fct_drop()
expr.all.mean.df %>% summary()
expr.all.mean.df %<>% left_join(dat.map %>% select(id, gene), by = c("sample" = "id")) 

# ================================================================ #

#### filter tissues and adapt phenotype info #######################
# filter for tissues to consider
tissues2consider <-
  expr.all.mean.df %>% filter(sampletype == "normalAll") %>%
  pull(site) %>% unique() %>% fct_drop() %>% levels()

dat.pheno %<>% filter(primary_site %in% tissues2consider)
dat.pheno$primary_site %<>% fct_drop()

# add information about patient id
dat.pheno %<>% mutate(patientID = sample) %>% relocate(patientID, .after = sample)
dat.pheno %<>% mutate(
  patientID = if_else(
    dat.pheno$sample_type != "Normal Tissue",
    dat.pheno$sample %>% str_remove("-[0-9]{2}$"),
    dat.pheno$sample %>% str_extract("^GTEX-[^-]+")
  )
)
if (SAVE_OUT) {
  dat.pheno %>%
    write_tsv(file = paste0(RES_PATH, "xena_infoSamples_", DATE_STR, ".tsv"))
}


# n over site and types
dat.pheno.info <-
  dat.pheno %>% group_by(primary_site, sample_type) %>%
  dplyr::summarise(n = dplyr::n())
dat.pheno$primary_site %<>% str_replace_all(" ", "_") %>% as_factor() %>% fct_drop()
dat.pheno.info %<>% mutate(descr = paste(primary_site, sample_type, sep = "_") %>%
                             str_replace_all(" ", "_"))
dat.pheno.info.patientsBased <- dat.pheno %>%
  dplyr::select(patientID, primary_site, sample_type) %>% distinct() %>%
  group_by(primary_site, sample_type) %>%
  dplyr::summarise(n_Patients = dplyr::n())
dat.pheno.info.patientsBased %<>%
  mutate(descr = paste(primary_site, sample_type, sep = "_") %>%
           str_replace_all(" ", "_"))
dat.pheno.info %<>% left_join(dat.pheno.info.patientsBased %>% ungroup() %>%
                                dplyr::select(n_Patients, descr),
                              by = "descr")
dat.pheno.info %>% print(n = 64)

if (SAVE_OUT) {
  dat.pheno.info %>%
    write_tsv(file = paste0(RES_PATH, "xena_infoSampleNo.tsv"))
}
# ================================================================ #


#### check sample and patient numbers ##############################
# total sample number
dat.pheno %>% 
  dplyr::summarise(n = dplyr::n())
# total patient number
dat.pheno %>% dplyr::select(patientID) %>% unique() %>% 
  dplyr::summarise(n = dplyr::n())
# sample number over types
dat.pheno %>% group_by(sample_type) %>% 
  dplyr::summarise(n = dplyr::n())
# patient number over types
dat.pheno %>% dplyr::select(patientID, sample_type) %>% distinct() %>% 
  group_by(sample_type) %>% 
  dplyr::summarise(n = dplyr::n())
# how many unique TCGA Solid Tissue Normal IDs compared to tumor?
setdiff(dat.pheno %>% filter(sample_type == "Solid Tissue Normal") %>% 
          pull(patientID), dat.pheno %>% filter(sample_type == "Primary Tumor") %>% 
          pull(patientID)) %>% length()
# sample n over TCGA patient IDs
dat.pheno %>% filter(sample_type != "Normal Tissue") %>% 
  pull(patientID) %>% unique() %>% length()
dat.pheno %>% filter(sample_type == "Normal Tissue") %>% 
  pull(patientID) %>% 
  unique() %>% length() # %>% table()
# sample n over site
dat.pheno %>% group_by(primary_site) %>% 
  dplyr::summarise(n = dplyr::n()) %>% print(n=25) %>% summary()
# patient n over site
dat.pheno %>% dplyr::select(patientID, primary_site) %>% distinct() %>% 
  group_by(primary_site) %>% 
  dplyr::summarise(n = dplyr::n()) %>% print(n=25) %>% summary()
# ================================================================ #


#### breast cancer #################################################
dat.pheno %>% filter(primary_site == "Breast") %>% summary()

expr.log2fc.mean.breast.tumorVSnormalAll <- get_fc(cond1 = "Breast_Primary_Tumor",
                                                     cond2 = c("Breast_normalAll"))
expr.log2fc.mean.breast.tumorVSnormalAll$sample %<>% str_remove(".[0-9]+$")

n = n_salt_up # use alle upregulated DEGs

p.breast <- getMyBarcode(pathway = getTopSaltGenes(df = dat.salt.upDEGs, 
                                                   n = n, id = "ensembl"),
                        stats = ( expr.log2fc.mean.breast.tumorVSnormalAll %>%
                                    dplyr::select(sample, log2fc) %>%
                                    deframe() ),
                        n1 = ( dat.pheno.info %>%
                                 filter(
                                   primary_site == "Breast" &
                                     sample_type == "Primary Tumor"
                                 ) %>%
                                 pull(n_Patients) ),
                        n2 = ( dat.pheno.info %>%
                                 filter(
                                   primary_site == "Breast" &
                                     sample_type != "Primary Tumor"
                                 ) %>%
                                 pull(n_Patients) %>% sum() ),
                        scoreType = "neg",
                        maxGS = 4000,
                        minGS = 10,
                        nPerm = 1e4,
                        title = "Breast"
                        
)

if (SAVE_OUT) {
  p.breast %>%
    cowplot::save_plot(
      filename = paste0(
        "breast/log2fc_n",
        n,
        "_tumorVSnormalAll_",
        "barcode_",
        DATE_STR,
        ".pdf"
      ),
      base_height = 6
    )
}
# ================================================================ #


#### all other tissues #############################################
n = n_salt_up
tissues2consider <-
  expr.all.mean.df %>% filter(sampletype == "normalAll" &
                                site != "Breast") %>%
  pull(site) %>% unique() %>% fct_drop() %>% levels()
plots.tissues <- list()

for (i in 1:length(tissues2consider)) {
  expr.log2fc.dummy <-
    get_fc(
      cond1 = paste0(tissues2consider[i], "_Primary_Tumor"),
      cond2 = paste0(tissues2consider[i], "_normalAll")
    )
  expr.log2fc.dummy$sample %<>% str_remove(".[0-9]+$")
  
  plots.tissues[[tissues2consider[i]]] <- getMyBarcode(
    pathway = getTopSaltGenes(df = dat.salt.upDEGs,
                              n = n, id = "ensembl"),
    stats = (
      expr.log2fc.dummy %>%
        dplyr::select(sample, log2fc) %>%
        deframe()
    ),
    n1 = (
      dat.pheno.info %>%
        filter(
          primary_site == tissues2consider[i] &
            sample_type == "Primary Tumor"
        ) %>%
        pull(n_Patients)
    ),
    n2 = (
      dat.pheno.info %>%
        filter(
          primary_site == tissues2consider[i] &
            sample_type != "Primary Tumor"
        ) %>%
        pull(n_Patients) %>% sum()
    ),
    scoreType = "pos",
    maxGS = 4000,
    minGS = 10,
    nPerm = 1e4,
    title = tissues2consider[i]
    
  )
}

p.cow <- cowplot::plot_grid(plotlist = plots.tissues, align = "hv", 
                            ncol = 4
                            )
save_plot(plot = p.cow, filename = 
            paste0(RES_PATH, "non_breast/facetted_exp1_", DATE_STR, ".pdf"),
          base_asp = 1/(sqrt(2)),
          base_height = 30
          )
# ================================================================ #
