#Pipeline:

#- Data Clean
#- Imputation (KNN-TN)
#- Univariate Analysis (Wilcoxon, PCA (fold-change), Loading, Volcano (fold change + classes) - see double bounds
#- Multivariate Analysis (OPLS-DA, WGCNA)
#- Linear Regression Models (volcano, forest plot corr)
#- Enrichment Analysis

# 
#======================================
#OPLS-DA = 
#- Score plot 
#- S Plot
#- Permutation Plot

#WGCNA = 
#- lipid modules Plot
#- Network modules-plot (circular)
#- Modules Eigen correlation with outcomes (Heatmap)
#- Network module principal (p.ex: pink, ou um para cada tipo de outcome)
#- Membership plot (scatter plot, um para cada module)
#####################################################################################

#Clear all
rm(list = ls())
library(maplet)
library(tidyverse)

#purrr::zap()
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/')
#file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "maplet")
file_data <- 'C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/data_Plasma_RAW_v2.xlsx'

# STARTING ----------------------------------------------------

# Load Data ----
D <-
  # load data - this function loads the assay data only
  mt_load_xls(file=file_data, sheet="Sheet3", samples_in_row=T, id_col="RID") %>%
  # load metabolite (rowData) annotations
  #mt_anno_xls(file=file_data, sheet="lipid",anno_type="features", anno_id_col="name", data_id_col = "Lipids") %>%
  # load clinical (colData) annotations
  mt_anno_xls(file=file_data, sheet="pheno", anno_type="samples", anno_id_col ="RID", data_id_col ="RID") %>%
  # record raw data info
  mt_reporting_data()%>%
  {.}

# Normalization ----
D %<>%
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  # normalize abundances using probabilistic quotient
  mt_pre_norm_quot() %>%
  # post-normalization sample boxplot
  {.}


# Data Transfomration ----
#D %<>%
#  mt_reporting_heading(heading = "Data Transformation", lvl=2) %>%
  # Log2 transformation
#  mt_pre_trans_log() %>%
  # scale and center data
#  mt_pre_trans_scale() %>%
#  {.}

# Imputation ----
D %<>%
  mt_reporting_heading(heading = "Imputation", lvl = 2) %>%
  #mt_plots_sample_boxplot(title = "Before imputation") %>%
  #mt_pre_impute_knn() %>%
  {.}


D <- D %>%
  mt_write_se_xls(file = "PreprocessedData_done2.xlsx") %>%
  {.}


