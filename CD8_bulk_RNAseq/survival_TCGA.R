#===================================================================================================================================#
# Title: survival analysis by Kaplan-Meier analysis using TCGA data
# Author: Sascha Schäuble
# Figures: Figure 6h, Extended Data 6i, Extended Data 23d
#===================================================================================================================================#

library("tidyverse")
library("magrittr")
library("data.table")
library("writexl")
library("rstatix")
library("ggpubr")
library("ggsci")
library("janitor")
library("survival")
library("survminer")
library("cowplot")

DATE_STR <- format(Sys.time(), "%y%m%d")

SAVE_OUT <- F

FN1 <- "gencode.v23.annotation.gene.probemap"

# ================================================================ #


#### Functions #####################################################

#'
#' @description: convenience function to make sure, we always use the same configs - use after loading data
#'
#' @param variables
#'
#' @return variables
refresh_base_config <- function() {
  PROJ_PATH <<-
    "/home/schaeuble/projects/funginet_inf/C7_Zielinski/EX0009/"
  DAT_PATH <<- paste0(PROJ_PATH, "dat/")
  RES_PATH <<- paste0(PROJ_PATH, "res/goi/ATP1A1/")
  
  DATE_STR <<- format(Sys.time(), "%y%m%d")
}
# ================================================================ #

#'
#' @description: convenience function to run KM analysis for both NFAT5 and SGK1
#'
#' @param sfitObj sfit object coming from surv
#' @param TITLE char
#' @param SUBTITLE char
#' @param CAPTION char
#' @param LEGENDTITLE char
#' @param LEGENDLABSdf char vec with length number of categories to differ (normally 2)
#' @param XLAB char
#' @param DATA df used for creating the fit
#' @param P char p value of computed survival or T to use default
#' @param PMETH char method to compute survival or T to use default
#'
#'
#' @return ggsurvplot obj
#'
get_km_plot <- function(sfitObj,
                        TITLE,
                        SUBTITLE = NULL,
                        CAPTION = NULL,
                        LEGENDTITLE = NULL,
                        LEGENDLABS = NULL,
                        XLAB,
                        DATA,
                        P = T,
                        PMETH = T,
                        CENSOR = F) {
  p <-
    ggsurvplot(
      fit = sfitObj,
      title = TITLE,
      legend.title = LEGENDTITLE,
      legend.labs = LEGENDLABS,
      pval.method = T,
      pval = P,
      risk.table = "abs_pct",
      risk.table.y.text = T,
      risk.table.height = 0.2,
      palette = c("black", "black"),
      linetype = c("solid", "dashed"),
      xlab = XLAB,
      data = DATA,
      censor = CENSOR,
      xlim = c(0, (sfitObj$time %>% max())),
      ggtheme = theme_pubr()
    )
  
  return(p)
}
# ================================================================ #

#'
#' @description: provide all necessary data as list; of note, different time formats
#' and censored data are compared correctly as of 240201
#'
#' @param df data frame: base data frame with "os" column for survival (1 dead, 0 alive)
#' @param COLS char vec: for which columns should KM data be generated?
#'
#' @return list of dfs with name of col from COLS per df
get_KM_data <- function(df, COLS, MINPROP = 0.1, FACTOR) {
  dat.km <- list()
  
  statusTimeMap <- c(
    "os_time" = "os",
    "new_tumor_event_dx_days_to" = "os",
    "dss_time" = "dss",
    "dfi_time" = "dfi",
    "pfi_time" = "pfi"
  )
  
  for (c in COLS) {
    dat.km[[c]] <- df %>%
      # drop_na(any_of(c(c, "os"))) %>%
      tidyr::drop_na(all_of(c(c, statusTimeMap[c] %>% unname))) %>%
      mutate(vital_status  = vital_status %>% as_factor() %>% fct_drop()) %>%
      dplyr::select(all_of(c(
        # c, "tumor_status", "os", "nfat5", "sgk1"
        c, "tumor_status", statusTimeMap[c] %>% unname, FACTOR
      ))) %>%
      rename("os" = (statusTimeMap[c] %>% unname))
    dat.km[[c]] <- dat.km[[c]] %>% cbind(
      surv_categorize(
        x = surv_cutpoint(
          data = dat.km[[c]],
          time = c,
          event = "os",
          variables = c(FACTOR),
          minprop = MINPROP
        ),
        variables = c(FACTOR),
        labels = c("low", "high")
      ) %>%
        as_tibble() %>%
        dplyr::select(c(FACTOR)) %>% dplyr::rename("optCut" = FACTOR)
    )
  }
  
  return(dat.km)
}

# ================================================================ #

#'
#' @description: convenience function to run one full KM analysis given a filtered dataset
#'
#' @param LIST.DF list of data frame: filtered datasets to run analysis with; name of list items will influence results labels
#' @param SAVE
#' @param RES_PATH
#' @param PREFIX
#' @param SURVCUT
#'
#'
#' @return list of sfit and p.surv fit list objects
#'
get_km_analysis <- function(LIST.DF,
                            SAVE = F,
                            SAVEPNG = F,
                            RES_PATH = "./",
                            PREFIX = "",
                            METHOD = "survdiff",
                            WTumSTATUS = F,
                            SURVCUT = c(500, 1000, 2000, 10000),
                            MODCENSOR = F,
                            WIDTH = 11,
                            HEIGHT = 10,
                            FACTOR) {
  sfit <- list()
  sfit.stats <- list()
  p.surv <- list()
  
  stepi <- 0
  pb <-
    txtProgressBar(
      min = 0,
      max = length(names(LIST.DF)) * length(SURVCUT) * if_else(WTumSTATUS, 2, 1),
      initial = stepi,
      style = 3
    )
  
  
  for (c in names(LIST.DF)) {
    for (i in SURVCUT) {
      if (MODCENSOR) {
        dummy <- LIST.DF[[c]] %>%
          mutate(os = if_else(.[[c]] > i, 0, os)) %>%
          mutate(!!c := if_else(.[[c]] > i, i, .[[c]]))
      } else {
        dummy <- LIST.DF[[c]] %>%
          filter(.[[c]] <= i)
      }
      
      dummy %<>% dplyr::mutate(dummyTime = dummy[[c]])
      
      ### Optimal cutoff as provided by surv_cutpoint
      # only expression status
      sfit[[paste(c, FACTOR, "_optCut", i, sep = "_")]] <-
        survfit(Surv(
          time = dummyTime,
          event = os,
          type = "right"
        ) ~ optCut,
        data = dummy)
      
      survdiff(Surv(dummyTime, os) ~ optCut, data = dummy)
      
      sfit.stats[[paste(c, FACTOR, "_optCut", i, sep = "_")]] <-
        survminer::surv_pvalue(fit = sfit[[paste(c, FACTOR, "_optCut", i, sep = "_")]],
                               method = METHOD,
                               data = dummy)
      
      p.surv[[paste(c, FACTOR, "_optCut", i, sep = "_")]] <-
        get_km_plot(
          sfitObj = sfit[[paste(c, FACTOR, "_optCut", i, sep = "_")]],
          TITLE = paste0(FACTOR, " separated by optimal cutoff"),
          SUBTITLE = c,
          LEGENDTITLE = "",
          CAPTION = paste0("Survival <", i, " days"),
          XLAB = c,
          DATA = dummy,
          P = paste0(sfit.stats[[paste(c, FACTOR, "_optCut", i, sep = "_")]]$method, "\n",
                     sfit.stats[[paste(c, FACTOR, "_optCut", i, sep = "_")]]$pval.txt)
        )
      setTxtProgressBar(pb = pb, value = (stepi = stepi + 1))
      
      
      ## with tumor status
      if (WTumSTATUS) {
        dummy %<>% drop_na(any_of(c("tumor_status")))
        sfit[[paste(c, FACTOR, "_optCut_wTStat", i, sep = "_")]] <-
          survfit(Surv(
            time = dummyTime,
            event = os,
            type = "right"
          ) ~ optCut + tumor_status,
          data = dummy)
        sfit.stats[[paste(c, FACTOR, "_optCut_wTStat", i, sep = "_")]] <-
          survminer::surv_pvalue(fit = sfit[[paste(c, FACTOR, "_optCut_wTStat", i, sep = "_")]],
                                 method = METHOD,
                                 data = dummy)
        p.surv[[paste(c, FACTOR, "_optCut_wTStat", i, sep = "_")]] <-
          get_km_plot(
            sfitObj = sfit[[paste(c, FACTOR, "_optCut_wTStat", i, sep = "_")]],
            TITLE = paste0(FACTOR, " separated by optimal cutoff"),
            SUBTITLE = c,
            CAPTION = paste0("Survival <", i, " days"),
            XLAB = c,
            DATA = dummy,
            P = paste0(sfit[[paste(c, FACTOR, "_optCut_wTStat", i, sep = "_")]]$method, "\n",
                       sfit[[paste(c, FACTOR, "_optCut_wTStat", i, sep = "_")]]$pval.txt)
          )
        setTxtProgressBar(pb = pb, value = (stepi = stepi + 1))
        
      }
    }
  }
  close(pb)
  
  res <- list(sfit, sfit.stats, p.surv)
  
  if (SAVE) {
    print("Saving result stats...")
    write_xlsx(
      x = rbindlist(sfit.stats, idcol = "mode"),
      path = paste0(
        RES_PATH,
        PREFIX,
        "km_",
        FACTOR,
        "_",
        DATE_STR,
        ".xlsx"
      )
    )
    print("Saving result plots...")
    ggexport(
      plotlist = p.surv,
      filename = paste0(RES_PATH, PREFIX, "km_", FACTOR, "_", DATE_STR, ".pdf"),
      width = WIDTH,
      height = HEIGHT
    )
    if (SAVEPNG) {
      ggexport(
        plotlist = p.surv,
        filename = paste0(
          RES_PATH,
          PREFIX,
          "km_",
          FACTOR,
          "_",
          DATE_STR,
          ".png"
        ),
        width = 1200,
        height = 800
      )
    }
    
  }
  
  return(res)
}

#'
#' @description: text
#'
#' @param variables
#'
#' @return variables
merge_meta <-
  function(DF = dat.expr,
           PHENO = dat.pheno,
           SURV = dat.survival) {
    dat.expr.long <-
      DF %>% dplyr::rename("ensembl" = "sample") %>%
      pivot_longer(-c(ensembl, gene), names_to = "sample") %>%
      dplyr::select(-ensembl) %>%
      pivot_wider(names_from = gene, values_from = value)
    # add metadata
    dat.expr.long.ext <- dat.expr.long %>% left_join(
      PHENO %>%
        dplyr::select(sample, patientID,
                      study, primary_site, sample_type, gender),
      by = "sample"
    ) %>%
      # relocate(value, .after = gene) %>%
      left_join(
        SURV %>%
          dplyr::select(
            sample,
            `cancer type abbreviation`,
            age_at_initial_pathologic_diagnosis,
            ajcc_pathologic_tumor_stage,
            vital_status,
            OS,
            DSS,
            DFI,
            DFI,
            PFI,
            tumor_status,
            new_tumor_event_dx_days_to,
            treatment_outcome_first_course,
            OS.time,
            DSS.time,
            DFI.time,
            PFI.time
          ),
        by = "sample"
      ) %>% clean_names()
    dat.expr.long.ext$ajcc_pathologic_tumor_stage %<>% as_factor() %>% fct_relevel(sort)
    
    return(dat.expr.long.ext)
  }
# ================================================================ #

# ================================================================ #


#### data wrangling ################################################
refresh_base_config()
## load meta data from TCGA
load("./survival/tcga_meta.RData")

## TCGA data from xena server will be loaded from RData files
# these are based on data available at
# https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=http%3A%2F%2F127.0.0.1%3A7222
# including
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_survival_data
# https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz

# set to FALSE if loading and processing original data
USERDATA <- TRUE

if (!USERDATA) {
  system("wget https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz")
  system("wget https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz")
  
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
  
  # add information about patient id
  dat.pheno %<>% mutate(patientID = sample) %>% relocate(patientID, .after = sample)
  dat.pheno %<>% mutate(
    patientID = if_else(
      dat.pheno$sample_type != "Normal Tissue",
      dat.pheno$sample %>% str_remove("-[0-9]{2}$"),
      dat.pheno$sample %>% str_extract("^GTEX-[^-]+")
    )
  )
  
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
  dat.pheno.info %<>% left_join(
    dat.pheno.info.patientsBased %>% ungroup() %>%
      dplyr::select(n_Patients, descr),
    by = "descr"
  )
  
  dat.map <- read_tsv(FN1)
  dat.expr <- data.table::fread(FN2)
  stopifnot(dplyr::setequal(dat.map$id, dat.expr$sample))
  
  dat.expr %<>% select(all_of(c("sample", dat.pheno$sample)))
}
# ================================================================ #

# ================================================================ #

#### pancreas KM ###################################################
# load data for convenience
if (USERDATA)
  load("./survival/tcga_xena_panc.RData")

dat.expr.all <- dat.expr

## NFAT5
dat.expr <-
  dat.expr.all %>% filter(sample == "ENSG00000102908") %>%
  mutate(gene = if (sample == "ENSG00000102908")
    "NFAT5")

dat.expr.panc <-
  merge_meta(DF = dat.expr) %>% dplyr::filter(primary_site == "Pancreas") %>%
  drop_na(nfat5) %>%
  mutate(vital_status = vital_status %>% as_factor())

dat.expr.panc %>%
  filter(sample_type == "Primary Tumor") %>%
  filter(tumor_status == "WITH TUMOR") %>%
  group_by(sample_type, vital_status, os, gender) %>% summarise(group_size = n())

dat.tmp <- get_KM_data(
  dat.expr.panc %>%
    filter(sample_type == "Primary Tumor") %>%
    filter(tumor_status == "WITH TUMOR"),
  COLS = c("new_tumor_event_dx_days_to"),
  FACTOR = "nfat5"
)

res.surv.paad.tumor.wTum.sel <- get_km_analysis(
  LIST.DF = dat.tmp,
  METHOD = "S1",
  SAVE = T,
  SURVCUT = c(10000),
  FACTOR = "nfat5",
  PREFIX = "panc_"
)


## ATP1A1
dat.expr <-
  dat.expr.all %>% filter(sample == "ENSG00000163399") %>%
  mutate(gene = if (sample == "ENSG00000163399")
    "ATP1A1")
dat.expr$gene

dat.expr.panc <-
  merge_meta(DF = dat.expr) %>% dplyr::filter(primary_site == "Pancreas") %>%
  drop_na(atp1a1) %>%
  mutate(vital_status = vital_status %>% as_factor())

dat.expr.panc %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(sample_type, vital_status, os, gender) %>% summarise(group_size = n())

dat.tmp <- get_KM_data(
  dat.expr.panc %>%
    filter(sample_type == "Primary Tumor"),
  #%>%
  COLS = c("new_tumor_event_dx_days_to"),
  FACTOR = "atp1a1"
)

res.surv.panc.tumorDead.wTum.sel <- get_km_analysis(
  LIST.DF = dat.tmp,
  METHOD = "S1",
  FACTOR = "atp1a1",
  SAVE = T,
  SURVCUT = c(10000),
  PREFIX = "panc_"
)

# ================================================================ #


#### breast cancer #################################################
if (USERDATA)
  load("./survival/tcga_xena_breast.RData")

dat.expr <- dat.expr %>% filter(sample == "ENSG00000102908") %>%
  mutate(gene = if (sample == "ENSG00000102908")
    "NFAT5")
dat.expr$gene

dat.expr.breast <-
  merge_meta(DF = dat.expr) %>% dplyr::filter(primary_site == "Breast") %>%
  drop_na(nfat5) %>%
  mutate(vital_status = vital_status %>% as_factor())

dat.expr.breast %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(sample_type, vital_status, os, gender) %>% summarise(group_size = n())

dat.tmp <- get_KM_data(
  dat.expr.breast %>%
    filter(sample_type == "Primary Tumor"),
  #%>%
  COLS = c("new_tumor_event_dx_days_to"),
  FACTOR = "nfat5"
)

res.surv.panc.tumorDead.wTum.sel <- get_km_analysis(
  LIST.DF = dat.tmp,
  METHOD = "S1",
  PREFIX = "breast_",
  FACTOR = "nfat5",
  SAVE = T,
  SURVCUT = c(10000)
)

# ================================================================ #
