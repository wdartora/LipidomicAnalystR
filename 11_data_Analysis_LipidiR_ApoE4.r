# Limpar o ambiente
rm(list = ls())

# Carregar bibliotecas necessárias
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(limma)
library(lipidr)
library(dplyr)
library(reshape2)

# Definir o diretório de trabalho
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/')

data_clin <- read.csv("data_clin.csv")

#####################################################################################################
# 1) Convert to ipidomicsExperiment
d <- as_lipidomics_experiment(read.csv("data_matrix.csv"))
d <- add_sample_annotation(d, data_clin)
d

# list non_parsed molecules
non_parsed_molecules(d)

# We can have a look at the clinical data, which was conveniently extracted from Metabolomics Workbench by lipidr.
colData(d)


#Performing PCA ----------------------->>>>
mvaresults = mva(d, measure="Area", method="PCA")

# Salvar o heatmap em um arquivo PNG
pca_plot <- plot_mva(mvaresults, color_by="ApoE4", components = c(1,2))
pca_plot <- pca_plot + theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("17-PCA_Plasma_ApoE4.png", plot = pca_plot, width = 6, height = 8)


#Performing Univariate Analysis ----------------------->>>>
# Definir os contrastes corretamente
two_group <- de_analysis(d, MCI-CN, AD-CN)
p_volcano<- plot_results_volcano(two_group)
p_volcano<- p_volcano+theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("18-Volcano_LipidiR_Plasma_ApoE4.png", plot = p_volcano, width = 6, height = 8)


#####################################################################################################
##################################################################################################

data_clin <- read.csv("data_clin_ApoE4.csv")


# 1) Convert to ipidomicsExperiment
d <- as_lipidomics_experiment(read.csv("data_matrix.csv"))
d <- add_sample_annotation(d, data_clin)
d

# list non_parsed molecules
non_parsed_molecules(d)

# We can have a look at the clinical data, which was conveniently extracted from Metabolomics Workbench by lipidr.
colData(d)


#Performing PCA ----------------------->>>>
mvaresults = mva(d, measure="Area", method="PCA")

# Salvar o heatmap em um arquivo PNG
pca_plot <- plot_mva(mvaresults, color_by="ApoE4", components = c(1,2))
pca_plot <- pca_plot + theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("17-PCA_Plasma_ApoE4.png", plot = pca_plot, width = 6, height = 8)


#Performing Univariate Analysis ----------------------->>>>
# Definir os contrastes corretamente
two_group <- de_analysis(d, Carrier-NonCarrier)
p_volcano<- plot_results_volcano(two_group)
p_volcano<- p_volcano+theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("18-Volcano_LipidiR_Plasma_ApoE4.png", plot = p_volcano, width = 6, height = 8)
