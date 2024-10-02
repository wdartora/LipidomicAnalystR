# Carregar pacotes necessários
library(openxlsx)
library(dplyr)
library(tidyr)

# Definir o diretório de trabalho
setwd("your_location/data")

# 1) --------------------------------------------------------------------------------->>>>
# Carregar os dados
df_AA <- read.xlsx("normalized_Plasma.xlsx")
df_BB <- read.xlsx("normalized_Serum.xlsx")
df_CC <- read.xlsx("pheno.xlsx", sheet = "pheno.wgcna")

# Identificar os IDs na base CC
ids_CC <- df_CC$Lipids

# Transpor as bases AA e BB para que os IDs se tornem linhas
df_AA_transposed <- df_AA %>%
  pivot_longer(cols = -Lipids, names_to = "ID", values_to = "Value") %>%
  pivot_wider(names_from = Lipids, values_from = Value)

df_BB_transposed <- df_BB %>%
  pivot_longer(cols = -Lipids, names_to = "ID", values_to = "Value") %>%
  pivot_wider(names_from = Lipids, values_from = Value)

# Manter apenas os IDs que estão presentes em todas as três bases
common_ids <- intersect(ids_CC, intersect(df_AA_transposed$ID, df_BB_transposed$ID))

# Filtrar as linhas das bases transpostas com base nos IDs comuns
filtered_AA <- df_AA_transposed %>%
  filter(ID %in% common_ids)

filtered_BB <- df_BB_transposed %>%
  filter(ID %in% common_ids)

filtered_CC <- df_CC %>%
  filter(Lipids %in% common_ids)

# Transpor de volta as bases filtradas para o formato original
filtered_AA <- filtered_AA %>%
  pivot_longer(cols = -ID, names_to = "Lipids", values_to = "Value") %>%
  pivot_wider(names_from = ID, values_from = Value)

filtered_BB <- filtered_BB %>%
  pivot_longer(cols = -ID, names_to = "Lipids", values_to = "Value") %>%
  pivot_wider(names_from = ID, values_from = Value)

# 2) --------------------------------------------------------------------------------->>>>
# Transformar a base de dados para o formato desejado
transform_data <- function(df) {
  df_transformed <- df %>%
    pivot_longer(cols = -Lipids, names_to = "RID", values_to = "Value") %>%
    mutate(RID = as.numeric(gsub("s_", "", RID))) %>%
    pivot_wider(names_from = Lipids, values_from = Value)
  return(df_transformed)
}

# Aplicar a transformação
df_AA_transformed <- transform_data(filtered_AA)
df_BB_transformed <- transform_data(filtered_BB)

# Criar um novo workbook
wb <- createWorkbook()

# Adicionar abas com dados filtrados e transformados
addWorksheet(wb, "wgcna.Plasma")
writeData(wb, "wgcna.Plasma", filtered_AA)

addWorksheet(wb, "Filtered_Plasma")
writeData(wb, "Filtered_Plasma", df_AA_transformed)

addWorksheet(wb, "wgcna.Serum")
writeData(wb, "wgcna.Serum", filtered_BB)

addWorksheet(wb, "Filtered_Serum")
writeData(wb, "Filtered_Serum", df_BB_transformed)

addWorksheet(wb, "Filtered_Pheno")
writeData(wb, "Filtered_Pheno", filtered_CC)

# Salvar o workbook
saveWorkbook(wb, "Filtered_Data.xlsx", overwrite = TRUE)
