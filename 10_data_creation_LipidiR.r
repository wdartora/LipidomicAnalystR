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
setwd('your_location/data/')

# Carregar os dados dos arquivos Excel
plasma_aa <- read.xlsx("Filtered_Data.xlsx", sheet = "wgcna.Plasma")
metadata <- read.xlsx("Filtered_Data.xlsx", sheet = "Filtered_Pheno")
classes <- read.xlsx("class_plasma.xlsx")

# Substituir os nomes dos lipídios pelo code_var
plasma_aa$Lipids <- classes$code_var[match(plasma_aa$Lipids, classes$Lipids)]

# Adicionar um espaço entre o nome do lipídio e o parêntese
#plasma_aa$Lipids <- gsub("([A-Za-z])\\(", "\\1 (", plasma_aa$Lipids)

# Replace  (COL) 
colnames(plasma_aa)[which(names(plasma_aa)=="Lipids")]<-"lipids"
colnames(metadata)[which(names(metadata)=="Lipids")]<-"Sample"

#Re-order
metadata<-metadata %>% 
  relocate(any_of(c("Sample")))

# Sub-set metadata
metadata<- (metadata[, c('Sample', 'ApoE4_cat', 'DX','SEX', 'Sex_text')])

#Create apoe4 char and DX char
metadata$ApoE4 <- ifelse(metadata$ApoE4 == 1, 'carrier', 'non-carrier')
metadata$DX_text <- ifelse(metadata$DX == 0, 'CN',
                ifelse(metadata$DX ==  1, 'MCI','AD'))

#Re-order
metadata<-metadata %>% 
  relocate(any_of(c("Sample", "DX_text","ApoE4","ApoE4_cat",  "DX",  "SEX", "Sex_text")))

# Save in csv
write.csv(plasma_aa, "data_matrix_v2.csv",row.names = FALSE)
write.csv(metadata, "data_clin_v2.csv",row.names = FALSE)

#data_clin_ApoE4
metadata$ApoE4 <- ifelse(metadata$ApoE4 == 'non-carrier', 'NonCarrier', 'Carrier')
#Re-order
metadata<-metadata %>% 
  relocate(any_of(c("Sample", "ApoE4","ApoE4","ApoE4_cat",  "DX",  "SEX", "Sex_text")))
write.csv(metadata, "data_clin_ApoE4_v2.csv",row.names = FALSE)
