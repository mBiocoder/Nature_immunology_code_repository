#===================================================================================================================================#
# Title: mTOR interaction network analysis with bulk RNA seq data and BIOGRID PPI information
# Author: Sascha Schäuble
# Figures: Extended Data 17b 
#===================================================================================================================================#

library("tidyverse")
library("magrittr")
library("here")
library("janitor")
library("ggpubr")
library("ComplexHeatmap")
library("circlize")
library("igraph")

DAT_PATH <- "./PPI_data/"

DATE_STR <- format(Sys.time(), "%y%m%d")

FN1 <- "cd8_salt_cd8_highsalt_vs_cd8_lowsalt_diff_expression.csv"
FN2 <- "cd8_salt_diff_expression.csv"
FN.ppi.1 <- "gProfiler_gconvert_all_interactors.csv"
FN.ppi.2 <- "mtor_BIOGRID-GENE-108757-4.4.229.tab3.txt"
FN.ppi.3 <- "sgk1_BIOGRID-GENE-112344-4.4.229.tab3.txt"

# ================================================================ #

#### Functions #####################################################
#'
#' @description: removing duplicates and resolving preferred gene names
#'
#' @param df - PPI BioGRID information for which respective column names are assumed
#' @param MAP - named char vector with names being the unpreferred and entries the preferred ID
#'
#' @return df
cleanBioGRID <- function(df, MAP = id.map) {
  # all to uppercase
  df %<>% 
    mutate(official_symbol_interactor_a = toupper(official_symbol_interactor_a)) %>%
    mutate(official_symbol_interactor_b = toupper(official_symbol_interactor_b))
  
  # correct some of the names
  df$official_symbol_interactor_a[which(df$official_symbol_interactor_a %in% names(MAP))] <-
    MAP[match(df$official_symbol_interactor_a[which(df$official_symbol_interactor_a %in% names(MAP))], names(MAP))]
  
  df$official_symbol_interactor_b[which(df$official_symbol_interactor_b %in% names(MAP))] <-
    MAP[match(df$official_symbol_interactor_b[which(df$official_symbol_interactor_b %in% names(MAP))], names(MAP))]
  
  # sort and reduce for duplicates
  df %<>% rowwise() %>%
    mutate(sorted_row = list(sort(
      c(
        official_symbol_interactor_a,
        official_symbol_interactor_b
      )
    ))) %>%
    select(sorted_row) %>% unnest_wider(col = sorted_row, names_sep = "fool") %>%
    rename(
      "official_symbol_interactor_a" = "sorted_rowfool1",
      "official_symbol_interactor_b" = "sorted_rowfool2"
    ) %>%
    filter(official_symbol_interactor_a != official_symbol_interactor_b) %>%
    distinct()
  
  return(df)

}
# ================================================================ #

#'
#' @description: gather correct meta data information
#'
#' @param DFEDGES
#'
#' @return df.meta
get_metaInfo <- function(DFEDGES, DFSTATS = dat.expr.ppi.stats) {
  # add some metadata
  df.meta <- tibble(symbol = union(
    DFEDGES$official_symbol_interactor_a,
    DFEDGES$official_symbol_interactor_b)
  )
  
  df.meta %<>%
    dplyr::left_join(
      y = DFSTATS %>%
        dplyr::select(
          genename,
          DE_cd8_highsalt_vs_cd8_lowsalt,
          LFC_cd8_highsalt_vs_cd8_lowsalt
        ) %>%
        rename(
          "symbol" = "genename",
          "salt_DEG_status" = "DE_cd8_highsalt_vs_cd8_lowsalt",
          "log2fc" = "LFC_cd8_highsalt_vs_cd8_lowsalt"
        ),
      by = "symbol"
    ) %>% distinct()
  return(df.meta)
}
# ================================================================ #


#'
#' @description: compute correlations based on a search SUBSTRING
#'
#' @param df
#' @param SUBSTRING - filter for dplyr::contains ; any of "high", "cd8" or "HvsL"
#'
#' @return df with augmented corr and corrColor information
get_cor <- function(df, SUBSTRING, EXPR = dat.expr) {
  df %<>%
    rowwise() %>%
    mutate(corr = tryCatch({
      cor(
        x = (
          EXPR %>% filter(genename == official_symbol_interactor_a) %>%
            dplyr::select(contains(SUBSTRING)) %>% unlist()
        ),
        y = (
          EXPR %>% filter(genename == official_symbol_interactor_b) %>%
            dplyr::select(contains(SUBSTRING)) %>% unlist()
        )
      )
    }, error = function(e) {
      # Handle the error, you can print a message or take other actions
      warning("Error occurred while computing correlation: ", conditionMessage(e))
      return(NA)  # Default value in case of an error
    })
    )
  color_fun.corr <- colorRamp2(breaks = c(-1,0,1), 
                               hcl_palette = "Blue-Red")
  # add color information to df
  df$corrColor <- color_fun.corr( df$corr )
  df$corrColor[is.na(df$corrColor)] <- "grey80"
  return(df)
}
# ================================================================ #

#'
#' @description: compute igraph object and necessary attributes for plotting
#'
#' @param df data frame with edge information in first two columns and optional meta info in further columns
#' @param QUANTS quantiles for node color cutoff (log2fc for high salt)
#' @param FILTER char if union computes "union" of node ids (OR for (MIN_EDGEDEGR, DEGs), otherwise intersect
#'
#' @return list with graphic objects for displaying
get_igraph <- function(df,
                       META,
                       MIN_EDGEDEGR = 1,
                       DEGS = F,
                       QUANTS = 0.95,
                       MIN_NODESIZE = 2,
                       SCALENODESIZE = T,
                       SCALEMODE = 1,
                       FILTER = "union") {
  g <- igraph::graph_from_data_frame(d = df,
                                     directed = F,
                                     vertices = META)
  
  ## filter information
  # edge degree
  edge_degrees <- igraph::degree(g, mode = "all")
  # DEGs
  if (DEGS) {
    ids.degs <- match(
      META %>% filter(
        salt_DEG_status == "upregulated" |
          salt_DEG_status == "downregulated"
      ) %>%
        filter(symbol %in% names(edge_degrees)) %>%
        pull(symbol),
      names(edge_degrees)
    )
    if (FILTER == "union") {
      node_ids <- union(which(edge_degrees >= MIN_EDGEDEGR), ids.degs)
    } else {
      node_ids <- intersect(which(edge_degrees >= MIN_EDGEDEGR), ids.degs)
    }
  } else {
    node_ids <- which(edge_degrees >= MIN_EDGEDEGR)
  }
  
  gsub <- induced_subgraph(graph = g, vids = node_ids)
  nodenames <- V(gsub)$name
  
  # Compute centrality measures
  degree_centrality <- igraph::degree(gsub)
  betweenness_centrality <- betweenness(gsub)
  closeness_centrality <- closeness(gsub)
  
  # 2. Scale node sizes relative to network and add minimal node size
  if (SCALENODESIZE) {
    if (SCALEMODE == 1) {
      # degree: btw 0-100 + minsize - default
      scaled_node_sizes <-
        (100 * (degree_centrality / sum(degree_centrality)) + MIN_NODESIZE)
    } else
      # degree - scaled by max degree and minsize
      if (SCALEMODE == 2) {
        scaled_node_sizes <-
          ((10 * degree_centrality / max(degree_centrality)) + MIN_NODESIZE)
      } else
        # betweenness: btw 0-100 + minsize
        if (SCALEMODE == 3) {
          scaled_node_sizes <-
            ((
              100 * betweenness_centrality / sum(betweenness_centrality)
            ) + MIN_NODESIZE)
        } else
          # betweenness - scaled by max degree and minsize
          if (SCALEMODE == 4) {
            scaled_node_sizes <-
              ((
                10 * betweenness_centrality / max(betweenness_centrality)
              ) + MIN_NODESIZE)
          } else
            # closeness: btw 0-100 + minsize
            if (SCALEMODE == 5) {
              scaled_node_sizes <-
                ((100 * closeness_centrality / sum(closeness_centrality)) + MIN_NODESIZE)
            } else
              # closeness - scaled by max degree and minsize
              if (SCALEMODE == 6) {
                scaled_node_sizes <-
                  ((10 * closeness_centrality / max(closeness_centrality)) + MIN_NODESIZE)
              }
    #
  } else {
    scaled_node_sizes <- MIN_NODESIZE
  }
  # create color ramp for remaining nodes
  tmp_dat <- META %>% filter(symbol %in% V(gsub)$name)
  tmpCut <-
    max(abs(quantile(
      tmp_dat$log2fc,
      na.rm = T,
      probs = c(1 - QUANTS, QUANTS)
    )))
  color_fun2 <- colorRamp2(breaks = c(-tmpCut, 0, tmpCut),
                           hcl_palette = "Green_Orange")
  # node color
  node_colors <-
    setNames(color_fun2(tmp_dat$log2fc), tmp_dat$symbol)
  node_label_color <-
    setNames(if_else(
      condition = (tmp_dat$symbol %in% c("MTOR", "SGK1")),
      true = "black",
      false = NA,
      missing = NA
    ),
    nm = tmp_dat$symbol)
  node_label <-
    setNames(if_else(
      condition = (tmp_dat$symbol %in% c("MTOR", "SGK1")),
      true = tmp_dat$symbol,
      false = NA,
      missing = NA
    ),
    nm = tmp_dat$symbol)
  
  vertex_colors <- node_colors[V(gsub)$name]
  
  return(
    list(
      igraph = gsub,
      NODENAMES = nodenames,
      VERTEXSIZE = scaled_node_sizes,
      LABEL = node_label,
      LABELCOL = node_label_color,
      VERTEXCOL = vertex_colors,
      EDGECOL = E(gsub)$corrColor,
      VERTEXCOLFUN = color_fun2
    )
  )
}
# ================================================================ #

#'
#' @description: provides standard output for stats of interest; currently with Wilc test
#'
#' @param variables
#'
#' @return list
get_stats <- function(DF.EXPR, DF.CORR) {
  
  stats <- list()
  # correlation
  stats[["corr"]] <- DF.CORR %>% ungroup() %>% 
    summarise(p = wilcox.test(corr, mu = 0)$p.value, 
              p.format = format(x = wilcox.test(corr, mu = 0)$p.value, digits = 2),
              yPos = max(corr, na.rm = T)+0.1)
  # log2fc
  stats[["expr"]] <- DF.EXPR %>% ungroup() %>% 
    summarise(p = wilcox.test(log2fc, mu = 0)$p.value, 
              p.format = format(x = wilcox.test(log2fc, mu = 0)$p.value, digits = 2),
              yPos = max(log2fc, na.rm = T)+0.1)
  # visually
  stats[["p.corr"]] <- DF.CORR %>% mutate(xlab = "bulk RNAseq") %>% 
    drop_na %>% ggviolin( x = "xlab",
                          y = "corr", 
                          add = c("boxplot"),
                          # add = c("jitter", "boxplot"),
                          font.tickslab = 14,
                          xlab = "", 
                          ylab = "Correlations of PPIs" ) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = paste("Wilc.-test, p =", p.format), x = 1, y = yPos+1), 
              inherit.aes = F, size = 6, data = stats[["corr"]])
  
  stats[["p.expr"]] <- DF.EXPR %>% mutate(xlab = "bulk RNAseq") %>% 
    drop_na %>% ggviolin( x = "xlab",
                          y = "log2fc", 
                          # add = c("jitter", "boxplot"), 
                          add = c("boxplot"), 
                          font.tickslab = 14,
                          xlab = "", 
                          ylab = "log2(fold-change)\nHigh vs. Low salt" ) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = paste("Wilc.-test, p =", p.format), x = 1, y = yPos+1), 
              inherit.aes = F, size = 6, data = stats[["expr"]])
  
  return(stats)
}
# ================================================================ #

#'
#' @description: convenience function to combine some plot elements
#'  currently adds the legends to the main plot
#'
#' @param variables
#'
#' @return variables
save_iGraphObj <- function(iGrObj,BP.EXPR,BP.CORR, PATH = "./", BASEFILENAME, SEED = 221705) {
  cairo_pdf(filename = paste0(PATH, BASEFILENAME, "_", DATE_STR, ".pdf"))
  set.seed(SEED)
  igraph::plot.igraph(x = iGrObj[["igraph"]],
                      # layout = layout_with_fr(ppi.mtor_sgk1$igraph),
                      vertex.size = iGrObj[["VERTEXSIZE"]],
                      vertex.label = iGrObj[["LABEL"]],
                      vertex.label.dist = 1,
                      vertex.label.font = 2,
                      vertex.label.cex = 0.8,
                      vertex.label.color = iGrObj[["LABELCOL"]],
                      vertex.color = iGrObj[["VERTEXCOL"]],
                      edge.width = 1,
                      edge.color = iGrObj[["EDGECOL"]]
  )
  draw(Legend(col_fun = iGrObj$VERTEXCOLFUN, title = "log2(FC)"),
       x = unit(0.2, "npc"), y = unit(0.7, "npc"), just = "left")
  draw(Legend(col_fun = colorRamp2(breaks = c(-1,0,1), 
                                   hcl_palette = "Blue-Red"), 
              title = "Correlation"),
       x = unit(0.2, "npc"), y = unit(0.5, "npc"), just = "left")
  dev.off()
  BP.EXPR %>% ggexport(
    filename = paste0(PATH, BASEFILENAME, "_", DATE_STR, "_expr.pdf"),
    width = 3.5, height = 4)
  BP.CORR %>% ggexport(
    filename = paste0(PATH, BASEFILENAME,"_", DATE_STR, "_corr.pdf"),
    width = 3.5, height = 4)
}
# ================================================================ #


# ================================================================ #


#### Data wrangling ################################################
# bulk expression data
dat.expr <- read_csv(FN1) %>% clean_names() %>% rename("ensembl" = "x1")
dat.stats <- read_csv(FN2)
dat.expr %<>% left_join(dat.stats %>% dplyr::select(c(geneid,genename)),
                        by = c("ensembl" = "geneid") ) %>% 
  relocate(genename, .after = "ensembl")
offset <- 1
dat.expr %<>% mutate(log2HvsL_rep1 = log2((cd8_highsalt_1+offset)/(cd8_lowsalt_1+offset)),
                     log2HvsL_rep2 = log2((cd8_highsalt_2+offset)/(cd8_lowsalt_2+offset)), 
                     log2HvsL_rep3 = log2((cd8_highsalt_3+offset)/(cd8_lowsalt_3+offset))) 

# pr-pr interaction network information
ppi.genes <- read_csv(paste0(DAT_PATH, FN.ppi.1), col_names = T)
biogrid.mtor <- read_tsv(paste0(DAT_PATH, FN.ppi.2)) %>% clean_names()
biogrid.sgk1 <- read_tsv(paste0(DAT_PATH, FN.ppi.3)) %>% clean_names()

id.map <- c(
  "ARSE" = "ARSL",
  "ATP5A1" = "ATP5F1A",
  "DDX58" = "RIGI",
  "Hsp22" = "HSPB8",
  "HYPM" = "H2AP",
  "IARS" = "IARS1",
  "ITFG3" = "FAM234A",
  "KIAA0226" = "RUBCN",
  "KIAA1429" = "VIRMA",
  "KIAA1715" = "LNPK",
  "Kif19a" = "KIF19",
  "LARS" = "LARS1",
  "MFSD4" = "MFSD4A",
  "SEPT2" = "SEPTIN2",
  "SMEK1" = "PPP4R3A",
  "TMEM206" = "PACC1",
  "TMEM261" = "DMAC1"
)
dat.expr %<>% mutate(genename = if_else(genename == "ARSE", "ARSL", genename))

# extract relevant expr info
dat.expr.ppi <- dat.expr %>% filter(ensembl %in% ppi.genes$converted_alias)
dat.expr.ppi.stats <- dat.stats %>% filter(geneid %in% ppi.genes$converted_alias)

## igraph data structure
ppi.edges <- list()
ppi.edges[["mtor"]] <- biogrid.mtor %>%
  dplyr::select(official_symbol_interactor_a,
                official_symbol_interactor_b) %>% distinct()
ppi.edges[["sgk1"]] <- biogrid.sgk1 %>%
  dplyr::select(official_symbol_interactor_a,
                official_symbol_interactor_b) %>% distinct()
ppi.edges[["mtor_sgk1"]] <- rbind(ppi.edges$mtor,
                                  ppi.edges$sgk1)
ppi.edges[["mtor_sgk1"]] <- cleanBioGRID(df = ppi.edges[["mtor_sgk1"]])

## BioGRID PPI | correlations
ppi.edges$mtor_sgk1 <- get_cor(df = ppi.edges$mtor_sgk1, SUBSTRING = "HvsL")
# ================================================================ #

#### igraph PPI ####################################################

## mtor-sgk1 network
ppi.mtor_sgk1.meta <- get_metaInfo(DFEDGES = ppi.edges$mtor_sgk1)
ppi.mtor_sgk1 <-
  get_igraph(
    df = ppi.edges$mtor_sgk1,
    META = ppi.mtor_sgk1.meta,
    MIN_EDGEDEGR = 1,
    DEGS = F,
    MIN_NODESIZE = 4, 
    QUANTS = c(1),
    SCALENODESIZE = T,
    SCALEMODE = 2,
    FILTER = "union"
  )
# stats 
ppi.mtor_sgk1.stats <- get_stats(DF.EXPR = ppi.mtor_sgk1.meta %>% 
                                   filter(symbol %in% V(ppi.mtor_sgk1$igraph)$name), 
                                 DF.CORR = ppi.edges$mtor_sgk1 %>% ungroup() %>%  
                                   filter(official_symbol_interactor_a %in% ppi.mtor_sgk1$NODENAMES |
                                            official_symbol_interactor_b %in% ppi.mtor_sgk1$NODENAMES)
                                 )

plotObj <- ppi.mtor_sgk1

save_iGraphObj(
  iGrObj = ppi.mtor_sgk1,
  BP.EXPR = ppi.mtor_sgk1.stats$p.expr,
  BP.CORR = ppi.mtor_sgk1.stats$p.corr,
  BASEFILENAME = "ppi_network_mtor_sgk1_bpWoJitter"
)

# ================================================================ #

