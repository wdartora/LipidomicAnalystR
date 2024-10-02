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

# 1) Convert to LipidomicsExperiment
d <- as_lipidomics_experiment(read.csv("data_matrix_v2.csv"))
d <- add_sample_annotation(d, data_clin)

# Listar moléculas não analisadas
non_parsed_molecules(d)

# Visualizar dados clínicos
colData(d)

# Verificar se as colunas necessárias estão presentes
if (!all(c("ApoE4_cat", "DX") %in% colnames(colData(d)))) {
  stop("As colunas 'ApoE4_cat' e 'DX' devem estar presentes nos dados clínicos.")
}

# Carregar o arquivo class_plasma.xlsx
class_plasma <- read.xlsx("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/class_plasma.xlsx")

# Dividir os dados em ApoE4 carriers e non-carriers
carriers <- d[, colData(d)$ApoE4_cat == 1]
non_carriers <- d[, colData(d)$ApoE4_cat == 0]

# Função para realizar as análises específicas e gerar gráficos
analyze_subset <- function(d, subset_name) {
  # Carregar a lista de classes de acyls
  acyl_classes <- read.xlsx("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/ACYL_class_list.xlsx")
  
  # Perform univariate analysis
  two_group <- tryCatch(de_analysis(d, MCI-CN, AD-CN), error = function(e) NULL)
  
  if (is.null(two_group)) {
    warning(paste("Análise univariada falhou para o grupo", subset_name))
    return(NULL)
  }
  
  # Reclassificar a coluna Class do two_group com base no code_var do class_plasma
  two_group <- merge(two_group, class_plasma[, c("code_var", "Class")], by.x = "Molecule", by.y = "code_var", all.x = TRUE)
  
  # Atualizar a coluna Class com os novos valores
  two_group$Class <- two_group$Class.y
  two_group <- two_group[, !names(two_group) %in% c("Class.x", "Class.y")]
  
  # Adicionar acyls ao two_group
  two_group <- merge(two_group, acyl_classes, by = "Class", all.x = TRUE)
  
  # Definir limiar de significância
  significance_threshold <- 0.05
  
  # Obter lista única de classes
  unique_classes <- unique(two_group$Class)
  
  # Função para gerar e salvar gráficos com significância
  plot_and_save <- function(data, class_name) {
    data$Significant <- data$adj.P.Val < significance_threshold

    p <- plot_trend(data) +
      theme_light() +
      scale_color_manual(name = "Comparison", values = c("AD - CN" = "red", "MCI - CN" = "blue")) +
      geom_point(data = data, aes(color = Significant), size = 3) +
      scale_color_manual(values = c("TRUE" = "#3cff00", "FALSE" = "grey"), guide = "none") +
      geom_text_repel(data = subset(data, Significant), aes(label = Molecule), size = 3)

    p <- p + scale_color_manual(name = "Comparison", values = c("AD - CN" = "red", "MCI - CN" = "blue"))

    ggsave(paste0("Class_chain_sig_", subset_name, "_", class_name, ".png"), plot = p, width = 8, height = 8, path = "\\\\Wat-em004jy3\\D\\william\\AAIC_2024\\lipidiR\\sig")
  }

  # Gerar e salvar gráficos para cada classe
  for (class_name in unique_classes) {
    class_data <- subset(two_group, Class == class_name)
    plot_and_save(class_data, class_name)
  }
}

# Realizar análise para cada subset
analyze_subset(carriers, "Carriers")
analyze_subset(non_carriers, "NonCarriers")

# Criar uma tabela de contagem de participantes por grupo ApoE4
apoE4_table <- as.data.frame(table(data_clin$ApoE4_cat))
colnames(apoE4_table) <- c("ApoE4_cat", "Count")

# Salvar a tabela em um arquivo Excel
write.xlsx(apoE4_table, "\\\\Wat-em004jy3\\D\\william\\AAIC_2024\\lipidiR\\sig\\ApoE4_counts.xlsx")
