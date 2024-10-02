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

data_clin <- read.csv("data_clin_v2.csv")

#####################################################################################################
# 1) Convert to LipidomicsExperiment
d <- as_lipidomics_experiment(read.csv("data_matrix_v2.csv"))
d <- add_sample_annotation(d, data_clin)
d

# list non_parsed molecules
non_parsed_molecules(d)


# All of them are Ceramides, written with full chemical name
# We can replace the first part with "Cer" using RegEx
#non_parsed <- non_parsed_molecules(d)
#new_names <- sub("CAR 14:0;O", "AC(12:0)", non_parsed)
#d <- update_molecule_names(d, old = non_parsed, new = new_names)

# We can check once again to make sure all molecules were parsed correctly
#non_parsed_molecules(d)


# We can have a look at the clinical data, which was conveniently extracted from Metabolomics Workbench by lipidr.
colData(d)


#Performing PCA ----------------------->>>>
mvaresults = mva(d, measure="Area", method="PCA")

# Salvar o heatmap em um arquivo PNG
pca_plot <- plot_mva(mvaresults, color_by="DX_text", components = c(1,3))
pca_plot <- pca_plot + theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("17-PCA_Plasma_DX.png", plot = pca_plot, width = 6, height = 8)


#Performing Univariate Analysis ----------------------->>>>
# Definir os contrastes corretamente
two_group <- de_analysis(d, MCI-CN, AD-CN)
p_volcano<- plot_results_volcano(two_group)
p_volcano<- p_volcano+theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("18-Volcano_LipidiR_Plasma_DX.png", plot = p_volcano, width = 6, height = 8)

# Multi-group comparison
multi_group <- de_design(d, ~ ApoE4)
significant_molecules(multi_group)

# Factorial analysis
factorial_de <- de_design(d, ~ ApoE4 + DX, coef = "DX")
significant_molecules(factorial_de)
# $DX
# [1] "CA"         "deDE(20:4)"

p<- plot_results_volcano(factorial_de)
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("19-Volcano_Sig_DX.png", plot = p, width = 6, height = 8)


# Multivariate analysis
# Orthogonal multivariate analysis
# OPLS and OPLS-DA can be performed to determine which lipids are associated with a group (y-variable) of interest. 
mvaresults = mva(d, method = "OPLS-DA", group_col = "DX_text", groups=c("CN", "MCI"))
w<-plot_mva(mvaresults, color_by="DX_text")
w<-w+theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("20-PCA_OPLS_DX.png", plot = w, width = 6, height = 8)


# Loading Plot -------------------------------------------->>>>>
load<- plot_mva_loadings(mvaresults, color_by="Class", top.n=10)
load<-load+theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("21-Loading_DX.png", plot = load, width = 6, height = 8)


# Top 10 lipids
pa<- top_lipids(mvaresults, top.n=10)
pa

#    filename                                  Molecule
#1  dataframe                      TG(O-54:2) [NL-18:1]
#2  dataframe                      TG(O-52:1) [NL-18:1]
#3  dataframe                      TG(O-52:1) [NL-16:0]
#4  dataframe                           LPC(15:0) [sn1]
#5  dataframe                               LPC(O-24:2)
#6  dataframe                             LPC(19:1) (c)
#7  dataframe       LPC(19:0) [sn1] (a) & LPC(19:0) [sn2] (b)
#8  dataframe                           LPC(20:1) [sn1]
#9  dataframe       LPC(17:1) [sn1] (a) & LPC(17:1) [sn2] (b)
#10 dataframe                           LPC(15:0) [sn2]

#Supervised multivariate analysis with continuous response variable
stage <- d$ApoE4_cat
stage <- as.numeric(stage)
stage

# We can see stage contains missing values. We should filter them out first.
d_filtered <- d[, !is.na(stage)]
stage <- stage[!is.na(stage)]

mvaresults = mva(d_filtered, method = "OPLS", group_col = stage )

#plot
p<-plot_mva(mvaresults)
p<-p+theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("22-Supervised_plot_DX.png", plot = p, width = 6, height = 8)


# Enrichment analysis ----------------->>>
enrich_results = lsea(two_group, rank.by = "logFC")
significant_lipidsets(enrich_results)
# > significant_lipidsets(enrich_results)
#$`AD - CN`
#[1] "total_cs_6"  "Class_TGO"   "Class_PE"    "Class_LPC"   "total_cl_18"
#[6] "total_cs_5"  "total_cs_7"  "total_cs_1"  "Class_LPCO"

#$`MCI - CN`
# [1] "Class_LPC"   "Class_PC"    "Class_PE"    "Class_TGO"   "total_cl_20"
# [6] "total_cl_54" "total_cl_52" "total_cl_17" "total_cl_18" "total_cs_0"
#[11] "Class_AC"    "total_cl_22" "Class_DE"    "Class_SM"    "total_cs_7"
#[16] "Class_PA"    "total_cl_15" "Class_CE"    "Class_LPE"   "total_cl_16"

# Visualization of enrichment analysis results. The enriched lipid classes are highlighted.
p<-plot_enrichment(two_group, significant_lipidsets(enrich_results), annotation="class")
p<-p+theme_light()
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("23-enrichment_plot_DX.png", plot = p, width = 32, height = 8)


# Lipid chain analysis --------------------------->>>>
p<-plot_trend(two_group)
p<-p+theme_light()+
    scale_color_manual(values = c("blue", "red"))
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("24-Chain_plot_DX.png", plot = p, width = 12, height = 8)


# Carregar a lista de classes de acyls
acyl_classes <- read.xlsx("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/ACYL_class_list.xlsx")

# Adicionar acyls ao two_group
two_group <- merge(two_group, acyl_classes, by = "Class", all.x = TRUE)

# Dividir two_group em quatro sub-conjuntos com base na coluna acyls
two_group_1 <- subset(two_group, acyls == 1)
two_group_2 <- subset(two_group, acyls == 2)
two_group_3 <- subset(two_group, acyls == 3)


# Função para gerar e salvar gráficos
plot_and_save <- function(data, subset_name) {
  p <- plot_trend(data) +
       theme_light() +
       scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish)
  
  # Salvar o gráfico
  setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
  ggsave(paste0("24-dClass_Length_chain_ADvsCN_", subset_name, ".png"), plot = p, width = 8, height = 8)
}

# Gerar e salvar gráficos para cada sub-conjunto
plot_and_save(two_group_1, "DX1")
plot_and_save(two_group_2, "DX2")
plot_and_save(two_group_3, "DX3")



# Obter lista única de classes
unique_classes <- unique(two_group$Class)

# Função para gerar e salvar gráficos
plot_and_save <- function(data, class_name) {
  p <- plot_trend(data) +
       theme_light() +
       scale_color_manual(values = c("AD - CN" = "red", "MCI - CN" = "blue"))
  
  # Salvar o gráfico
  setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/class_acyl_chain/')
  ggsave(paste0("Class_chain_length_chain_", class_name, ".png"), plot = p, width = 8, height = 8)
}

# Gerar e salvar gráficos para cada classe
for (class_name in unique_classes) {
  class_data <- subset(two_group, Class == class_name)
  plot_and_save(class_data, class_name)
}



# Definir limiar de significância
significance_threshold <- 0.05

# Obter lista única de classes
unique_classes <- unique(two_group$Class)

# Função para gerar e salvar gráficos
plot_and_save <- function(data, class_name) {
  # Adicionar coluna de significância
  data$Significant <- data$adj.P.Val < significance_threshold

  # Gerar o gráfico original usando plot_trend
  p <- plot_trend(data) +
       theme_light() +
       scale_color_manual(name = "Comparison", values = c("AD - CN" = "red", "MCI - CN" = "blue")) +
       geom_point(data = data, aes(color = Significant), size = 3) +
       scale_color_manual(values = c("TRUE" = "#3cff00", "FALSE" = "grey"), guide = "none") +  # Mantém a legenda dos grupos, mas não dos pontos
       geom_text_repel(data = subset(data, Significant), aes(label = Molecule), size = 3)

  # Adicionar novamente a escala de cores para as comparações
  p <- p + scale_color_manual(name = "Comparison", values = c("AD - CN" = "red", "MCI - CN" = "blue"))

  # Salvar o gráfico
  setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/class_acyl_chain/sig/')
  ggsave(paste0("Class_chain_sig_", class_name, ".png"), plot = p, width = 8, height = 8)
}

# Gerar e salvar gráficos para cada classe
for (class_name in unique_classes) {
  class_data <- subset(two_group, Class == class_name)
  plot_and_save(class_data, class_name)
}


#Performing Univariate Analysis ----------------------->>>>
# Definir os contrastes corretamente
A_group <- de_analysis(d,  AD-CN)
B_group <- de_analysis(d,  MCI-CN)

#AD vs CN
p<-plot_chain_distribution(A_group, contrast = NULL, measure = "logFC")
p<-p+theme_light() +
     scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish)

# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("25-Class_chain_inf_ADvsCN_DX.png", plot = p, width = 18, height = 24)


# Extrair informações de todos os lipídios
lipid_info <- data.frame(
  Molecule = A_group$Molecule,
  Class = A_group$Class,
  logFC = A_group$logFC,
  p.value = A_group$P.Value,
  adj.p.value = A_group$adj.P.Val,
  total_chain_length = A_group$total_cl,
  total_chain_unsaturation = A_group$total_cs
)

# Exportar para um arquivo Excel
write.xlsx(lipid_info, "C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/results/25_all_lipids_ADvsCN.xlsx")

# Exibir as primeiras linhas do dataframe exportado
head(lipid_info)




#MCI vs CN
p<-plot_chain_distribution(B_group, contrast = NULL, measure = "logFC")
p<-p+theme_light() +
     scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0, limits = c(-0.5, 0.5), oob = scales::squish)
# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("25-Class_chain_inf_MCIvsCN_DX.png", plot = p, width = 18, height = 24)

# Extrair informações de todos os lipídios
lipid_info <- data.frame(
  Molecule = B_group$Molecule,
  Class = B_group$Class,
  logFC = B_group$logFC,
  p.value = B_group$P.Value,
  adj.p.value = B_group$adj.P.Val,
  total_chain_length = B_group$total_cl,
  total_chain_unsaturation = B_group$total_cs
)

# Exportar para um arquivo Excel
write.xlsx(lipid_info, "C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/results/25_all_lipids_MCIvsCN.xlsx")

# Exibir as primeiras linhas do dataframe exportado
head(lipid_info)
