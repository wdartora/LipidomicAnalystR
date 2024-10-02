# Limpar o ambiente
rm(list = ls())

# Carregar bibliotecas necessárias
library(WGCNA)
library(pheatmap)
library(openxlsx)
########################################################
# 1) Plasma Module vs Serum Particle

# Configurações básicas
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Definir diretórios e carregar dados de expressão
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data")

# Carregar dados de soro a partir do arquivo Excel
serum_df <- read.xlsx("Filtered_Data.xlsx", sheet = "wgcna.Serum")

# Extrair a coluna de lipídios e os dados de expressão
lipids <- serum_df$Lipids
datExprSerum <- as.data.frame(t(serum_df[, -1]))
colnames(datExprSerum) <- lipids

# Carregar dados de plasma a partir de arquivo RData
load("03_output_Input_plasma.RData")  # Carregar dados de plasma
datExprPlasma = datExpr

# Carregar dados da rede e identificar módulos
load("04_output_networkConstruction_Plasma.RData")  # Carregar rede de plasma
moduleColorsPlasma = moduleColors
MEsPlasma = MEs

# Ordenar os eigengenes
MEsPlasma = orderMEs(MEsPlasma)

# Calcular a correlação entre os módulos de plasma e as features de soro
crossCorr = cor(datExprSerum, MEsPlasma, use = "p")

# Preparar a matriz de anotação
textMatrix = paste(signif(crossCorr, 2), "\n(", signif(corPvalueStudent(crossCorr, nrow(datExprSerum)), 1), ")", sep = "")
dim(textMatrix) = dim(crossCorr)

# Definir o diretório de saída
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots")

# Plotar o heatmap
png(file="12-cross_correlation_heatmap_particle.png", width=18, height=48, units="in", res=300)
pheatmap(crossCorr,
         display_numbers = textMatrix,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Cross-Correlation between Serum Features and Plasma Modules",
         fontsize_number = 5)
dev.off()




########################################################
# 2) Plasma Module vs Serum subClass
# Limpar o ambiente
rm(list = ls())

# Carregar bibliotecas necessárias
library(WGCNA)
library(pheatmap)
library(openxlsx)
library(dplyr)

# Configurações básicas
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Definir diretórios e carregar dados de expressão
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data")

# Carregar dados de soro a partir do arquivo Excel
serum_df <- read.xlsx("Filtered_Data.xlsx", sheet = "wgcna.Serum")

# Carregar informações sobre os grupos de lipídios
lipid_info <- read.xlsx("Class_Serum_list.xlsx")

# Combinar as informações dos grupos com a matriz de expressão de soro
serum_df <- serum_df %>%
  left_join(lipid_info, by = "Lipids")

# Calcular a média dos lipídios por grupo (classe)
serum_avg <- serum_df %>%
  select(-Biomarker_name, -Unit) %>%
  group_by(Group) %>%
  summarise(across(starts_with("s_"), ~mean(.x, na.rm = TRUE)))

# Transpor a matriz para que as amostras fiquem nas linhas
datExprSerum <- as.data.frame(t(serum_avg[, -1]))
colnames(datExprSerum) <- serum_avg$Group

# Carregar dados de plasma a partir de arquivo RData
load("03_output_Input_plasma.RData")  # Carregar dados de plasma
datExprPlasma = datExpr

# Carregar dados da rede e identificar módulos
load("04_output_networkConstruction_Plasma.RData")  # Carregar rede de plasma
moduleColorsPlasma = moduleColors
MEsPlasma = MEs

# Ordenar os eigengenes
MEsPlasma = orderMEs(MEsPlasma)

# Calcular a correlação entre os módulos de plasma e as médias dos grupos de lipídios de soro
crossCorr = cor(MEsPlasma, datExprSerum, use = "p")

# Preparar a matriz de anotação
textMatrix = paste(signif(crossCorr, 2), "\n(", signif(corPvalueStudent(crossCorr, nrow(datExprSerum)), 1), ")", sep = "")
dim(textMatrix) = dim(crossCorr)

# Definir o diretório de saída
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots")

# Plotar o heatmap
png(file="13-cross_correlation_heatmap_class.png", width=15, height=12, units="in", res=300)
pheatmap(crossCorr,
         display_numbers = textMatrix,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Cross-Correlation between Serum Lipid Subclass and Plasma Modules",
         fontsize_number = 10)
dev.off()



########################################################
# 3) ---------------------------------------->>>>
# Plasma Modules vs Serum Cluster 3 from UMAP
# Limpar o ambiente
rm(list = ls())

# Carregar bibliotecas necessárias
library(WGCNA)
library(pheatmap)
library(openxlsx)
library(dplyr)

# Configurações básicas
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Definir diretórios e carregar dados de expressão
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots")

# Carregar dados de abundância de lipídios a partir do arquivo Excel
lipid_abundance_df <- read.xlsx("Lipid_Abundance_in_Cluster_3_Grouped_Serum.xlsx")

setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data")

# Carregar dados de plasma a partir de arquivo RData
load("03_output_Input_plasma.RData")  # Carregar dados de plasma
datExprPlasma <- datExpr

# Carregar dados da rede e identificar módulos
load("04_output_networkConstruction_Plasma.RData")  # Carregar rede de plasma
moduleColorsPlasma <- moduleColors
MEsPlasma <- MEs

# Ordenar os eigengenes
MEsPlasma <- orderMEs(MEsPlasma)

# Ajustar os RIDs nos dados de abundância de lipídios
lipid_abundance_df$RID <- paste0("s_", lipid_abundance_df$RID)

# Reorganizar os dados de abundância para que cada RID tenha uma linha com médias de abundância por classe
serum_avg <- lipid_abundance_df %>%
  group_by(RID, Class) %>%
  summarise(Mean_Abundance = mean(Total_Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Class, values_from = Mean_Abundance)

# Verificar os RIDs presentes em ambos os conjuntos de dados
common_RIDs <- intersect(rownames(MEsPlasma), serum_avg$RID)
print(paste("Number of common RIDs:", length(common_RIDs)))

# Se houver RIDs comuns, continuar com a análise
if (length(common_RIDs) > 0) {
  # Filtrar os dados para manter apenas os RIDs comuns
  MEsPlasma <- MEsPlasma[common_RIDs, ]
  serum_avg <- serum_avg %>% filter(RID %in% common_RIDs)
  
  # Remover a coluna RID para manter a consistência com as dimensões dos módulos de plasma
  serum_avg_matrix <- serum_avg %>%
    select(-RID)
  
  # Verificar as dimensões
  print(dim(MEsPlasma))
  print(dim(serum_avg_matrix))
  
  # Calcular a correlação entre os módulos de plasma e as classes de lipídios do soro
  crossCorr <- cor(MEsPlasma, serum_avg_matrix, use = "p")
  
  # Preparar a matriz de anotação
  textMatrix <- paste(signif(crossCorr, 2), "\n(", signif(corPvalueStudent(crossCorr, nrow(MEsPlasma)), 1), ")", sep = "")
  dim(textMatrix) <- dim(crossCorr)
  
  # Definir o diretório de saída
  setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots")
  
  # Plotar o heatmap
  png(file = "14-cross_correlation_heatmap.png", width = 15, height = 12, units = "in", res = 300)
  pheatmap(crossCorr,
           display_numbers = textMatrix,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           cluster_rows = TRUE, cluster_cols = TRUE,
           main = "Cross-Correlation between Serum Lipid Classes (Cluster 3) and Plasma Modules",
           fontsize_number = 10)
  dev.off()
} else {
  print("Não há RIDs comuns entre os conjuntos de dados.")
}
