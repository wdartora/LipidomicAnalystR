# Limpar o ambiente
rm(list = ls())
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/")
library(openxlsx)
library(dplyr)


datExprPlasma <-  read.xlsx("data1_Plasma.xlsx", sheet = "plasma.wgcna")
datExprSerum <- read.xlsx("data2_Serum.xlsx", sheet = "serum.wgcna")



# Função para normalizar dados usando z-score
normalize_zscore <- function(data) {
  lipids <- data[, 1]  # Preservar a coluna de lipídios
  numeric_cols <- data[, -1]  # Selecionar apenas colunas numéricas
  normalized_data <- as.data.frame(scale(numeric_cols, center = TRUE, scale = TRUE))
  normalized_data <- cbind(lipids, normalized_data)  # Re-inserir a coluna de lipídios
  colnames(normalized_data)[1] <- colnames(data)[1]  # Manter o nome original da coluna de lipídios
  return(normalized_data)
}

# Normalizar os dados
normalizedPlasma <- normalize_zscore(datExprPlasma)
normalizedSerum <- normalize_zscore(datExprSerum)

head(normalizedPlasma)
head(normalizedSerum)
# Salvar os dados normalizados em novas planilhas
write.xlsx(normalizedPlasma, "normalized_Plasma.xlsx")
write.xlsx(normalizedSerum, "normalized_Serum.xlsx")
