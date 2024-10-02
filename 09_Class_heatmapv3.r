# Clear all
rm(list = ls())
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/')
library(openxlsx)
library(Hmisc)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)

# Defina uma paleta de cores personalizada
my_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(50)

# Ler dados brutos
assay_data <- read.xlsx("Filtered_Data.xlsx", sheet = "wgcna.Plasma")
classes <- read.xlsx("class_plasma.xlsx")

# Combinar os dados de abundância com as classes
data_long <- assay_data %>%
  pivot_longer(-Lipids, names_to = "Sample", values_to = "Abundance") %>%
  inner_join(classes, by = "Lipids")

# Calcular a média da abundância para cada classe
data_avg <- data_long %>%
  group_by(Class, Sample) %>%
  summarise(Average_Abundance = mean(Abundance, na.rm = TRUE)) %>%
  pivot_wider(names_from = Class, values_from = Average_Abundance)

# Transformar em matriz
data_avg_matrix <- as.matrix(data_avg[,-1])
rownames(data_avg_matrix) <- data_avg$Sample

# Calcular a matriz de correlação de Spearman entre as classes
cor_matrix <- rcorr(data_avg_matrix, type = "spearman")$r

# Tratar valores NA/NaN/Inf na matriz de correlação
cor_matrix[is.na(cor_matrix) | is.nan(cor_matrix) | is.infinite(cor_matrix)] <- 0

# Dendrograma
mat_cluster_cols <- hclust(dist(cor_matrix))

# Salvar o heatmap em um arquivo PNG
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
png(paste0("16_Plasma_heatmap_Class.png"), width = 8, height = 8, units = "in", res = 500)

# Heatmap
p <- pheatmap(
  mat = cor_matrix,
  color = inferno(10),
  border_color = NA,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize = 10,
  main = "Heatmap - Class Averages"
)

p
dev.off()
