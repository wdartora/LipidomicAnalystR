#################################################################################
#################################################################################
#################################################################################
# Carrier

# Limpa o ambiente
rm(list = ls())
# Carrega as bibliotecas necessárias
library(lme4)
library(lmerTest)  # Fornece valores-p para modelos mistos
library(openxlsx)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(sandwich)
library(msm)
library(lmtest)

# Define o diretório de trabalho e carrega os dados
setwd('your_location/data/')
lipid_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Plasma')
phenotype_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Pheno')
class_data <- read.xlsx('class_plasma.xlsx')
class_data$Variable <- class_data$Lipids


# Combinando as bases de dados
lipid_data$RID <- as.double(lipid_data$RID)
phenotype_data$RID <- as.double(phenotype_data$RID)
combined_data <- phenotype_data %>% left_join(lipid_data, by = "RID")

combined_data2 <- subset(combined_data, ApoE4_cat == 1)
combined_data2 <- subset(combined_data2, DX == 0)

ls(phenotype_data)

# Lista de desfechos a serem analisados-------------------------------------- |||||||
outcomes <- c("ADAS.Cog13", "MMSE", "CDRSB")

#-------------------------------------- |||||||-------------------------------------- |||||||
for (outcome in outcomes) {
  results_df <- data.frame(Variable = character(0), 
                           Estimate = numeric(0), 
                           SE = numeric(0), 
                           P.value = numeric(0), 
                           P.value.adjusted = numeric(0))
  
  variables <- setdiff(colnames(lipid_data), "RID")
  
  for (var in variables) {
    formula_str <- as.formula(paste(var, "~", outcome, "+ Education + Age + SEX + BMI +  statin"))
    model <- glm(formula_str, family = "gaussian"(link = "identity"),  data = combined_data2)
    
    model_summary <- summary(model)
    se <- sqrt(diag(vcovHC(model, type = "HC0")))
    est <- coef(model)
    
    temp_df <- data.frame(Variable = var,
                          Estimate = est[outcome], 
                          SE = se[outcome])
    
    
    temp_df$P.value <- 2 * pnorm(-abs(est[outcome]/se[outcome]))
    
    results_df <- rbind(results_df, temp_df)
  }
  
  results_df$P.value.adjusted <- p.adjust(results_df$P.value, method = "BH")
  results_df <- results_df[, c("Variable", "Estimate", "SE", "P.value", "P.value.adjusted")]
  # Combina os resultados com a base de dados class_data baseada na coluna 'Variable'
  results_df <- merge(results_df, class_data, by = "Variable", all.x = TRUE)
  
  # Salva os resultados em um arquivo Excel, com o nome baseado no desfecho analisado
  filename <- paste("GLM_ApoE4_carrier_DX_cat_CN_", outcome, "_adj.xlsx", sep = "")
  write.xlsx(results_df, file.path('your_location/results/DX/', filename))
}


#################################################################################

# Limpa o ambiente
rm(list = ls())
# Carrega as bibliotecas necessárias
library(lme4)
library(lmerTest)  # Fornece valores-p para modelos mistos
library(openxlsx)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(sandwich)
library(msm)
library(lmtest)

# Define o diretório de trabalho e carrega os dados
setwd('your_location/data/')
lipid_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Plasma')
phenotype_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Pheno')
class_data <- read.xlsx('class_plasma.xlsx')
class_data$Variable <- class_data$Lipids


# Combinando as bases de dados
lipid_data$RID <- as.double(lipid_data$RID)
phenotype_data$RID <- as.double(phenotype_data$RID)
combined_data <- phenotype_data %>% left_join(lipid_data, by = "RID")

combined_data2 <- subset(combined_data, ApoE4_cat == 1)
combined_data2 <- subset(combined_data2, DX == 1)

ls(phenotype_data)

# Lista de desfechos a serem analisados-------------------------------------- |||||||
outcomes <- c("ADAS.Cog13", "MMSE", "CDRSB")

#-------------------------------------- |||||||-------------------------------------- |||||||
for (outcome in outcomes) {
  results_df <- data.frame(Variable = character(0), 
                           Estimate = numeric(0), 
                           SE = numeric(0), 
                           P.value = numeric(0), 
                           P.value.adjusted = numeric(0))
  
  variables <- setdiff(colnames(lipid_data), "RID")
  
  for (var in variables) {
    formula_str <- as.formula(paste(var, "~", outcome, "+ Education + Age + SEX + BMI +  statin"))
    model <- glm(formula_str, family = "gaussian"(link = "identity"),  data = combined_data2)
    
    model_summary <- summary(model)
    se <- sqrt(diag(vcovHC(model, type = "HC0")))
    est <- coef(model)
    
    temp_df <- data.frame(Variable = var,
                          Estimate = est[outcome], 
                          SE = se[outcome])
    
    
    temp_df$P.value <- 2 * pnorm(-abs(est[outcome]/se[outcome]))
    
    results_df <- rbind(results_df, temp_df)
  }
  
  results_df$P.value.adjusted <- p.adjust(results_df$P.value, method = "BH")
  results_df <- results_df[, c("Variable", "Estimate", "SE", "P.value", "P.value.adjusted")]
  # Combina os resultados com a base de dados class_data baseada na coluna 'Variable'
  results_df <- merge(results_df, class_data, by = "Variable", all.x = TRUE)
  
  # Salva os resultados em um arquivo Excel, com o nome baseado no desfecho analisado
  filename <- paste("GLM_ApoE4_carrier_DX_cat_MCI_", outcome, "_adj.xlsx", sep = "")
  write.xlsx(results_df, file.path('your_location/results/DX/', filename))
}



#################################################################################

# Limpa o ambiente
rm(list = ls())
# Carrega as bibliotecas necessárias
library(lme4)
library(lmerTest)  # Fornece valores-p para modelos mistos
library(openxlsx)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(sandwich)
library(msm)
library(lmtest)

# Define o diretório de trabalho e carrega os dados
setwd('your_location/data/')
lipid_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Plasma')
phenotype_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Pheno')
class_data <- read.xlsx('class_plasma.xlsx')
class_data$Variable <- class_data$Lipids


# Combinando as bases de dados
lipid_data$RID <- as.double(lipid_data$RID)
phenotype_data$RID <- as.double(phenotype_data$RID)
combined_data <- phenotype_data %>% left_join(lipid_data, by = "RID")

combined_data2 <- subset(combined_data, ApoE4_cat == 1)
combined_data2 <- subset(combined_data2, DX == 2)

ls(phenotype_data)

# Lista de desfechos a serem analisados-------------------------------------- |||||||
outcomes <- c("ADAS.Cog13", "MMSE", "CDRSB")

#-------------------------------------- |||||||-------------------------------------- |||||||
for (outcome in outcomes) {
  results_df <- data.frame(Variable = character(0), 
                           Estimate = numeric(0), 
                           SE = numeric(0), 
                           P.value = numeric(0), 
                           P.value.adjusted = numeric(0))
  
  variables <- setdiff(colnames(lipid_data), "RID")
  
  for (var in variables) {
    formula_str <- as.formula(paste(var, "~", outcome, "+ Education + Age + SEX + BMI +  statin"))
    model <- glm(formula_str, family = "gaussian"(link = "identity"),  data = combined_data2)
    
    model_summary <- summary(model)
    se <- sqrt(diag(vcovHC(model, type = "HC0")))
    est <- coef(model)
    
    temp_df <- data.frame(Variable = var,
                          Estimate = est[outcome], 
                          SE = se[outcome])
    
    
    temp_df$P.value <- 2 * pnorm(-abs(est[outcome]/se[outcome]))
    
    results_df <- rbind(results_df, temp_df)
  }
  
  results_df$P.value.adjusted <- p.adjust(results_df$P.value, method = "BH")
  results_df <- results_df[, c("Variable", "Estimate", "SE", "P.value", "P.value.adjusted")]
  # Combina os resultados com a base de dados class_data baseada na coluna 'Variable'
  results_df <- merge(results_df, class_data, by = "Variable", all.x = TRUE)
  
  # Salva os resultados em um arquivo Excel, com o nome baseado no desfecho analisado
  filename <- paste("GLM_ApoE4_carrier_DX_cat_AD_", outcome, "_adj.xlsx", sep = "")
  write.xlsx(results_df, file.path('your_location/results/DX/', filename))
}


#################################################################################
#################################################################################
#################################################################################
# Non-Carrier

#################################################################################

# Limpa o ambiente
rm(list = ls())
# Carrega as bibliotecas necessárias
library(lme4)
library(lmerTest)  # Fornece valores-p para modelos mistos
library(openxlsx)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(sandwich)
library(msm)
library(lmtest)

# Define o diretório de trabalho e carrega os dados
setwd('your_location/Proj12/data/')
lipid_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Plasma')
phenotype_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Pheno')
class_data <- read.xlsx('class_plasma.xlsx')
class_data$Variable <- class_data$Lipids


# Combinando as bases de dados
lipid_data$RID <- as.double(lipid_data$RID)
phenotype_data$RID <- as.double(phenotype_data$RID)
combined_data <- phenotype_data %>% left_join(lipid_data, by = "RID")

combined_data2 <- subset(combined_data, ApoE4_cat == 0)
combined_data2 <- subset(combined_data2, DX == 0)

ls(phenotype_data)

# Lista de desfechos a serem analisados-------------------------------------- |||||||
outcomes <- c("ADAS.Cog13", "MMSE", "CDRSB")

#-------------------------------------- |||||||-------------------------------------- |||||||
for (outcome in outcomes) {
  results_df <- data.frame(Variable = character(0), 
                           Estimate = numeric(0), 
                           SE = numeric(0), 
                           P.value = numeric(0), 
                           P.value.adjusted = numeric(0))
  
  variables <- setdiff(colnames(lipid_data), "RID")
  
  for (var in variables) {
    formula_str <- as.formula(paste(var, "~", outcome, "+ Education + Age + SEX + BMI +  statin"))
    model <- glm(formula_str, family = "gaussian"(link = "identity"),  data = combined_data2)
    
    model_summary <- summary(model)
    se <- sqrt(diag(vcovHC(model, type = "HC0")))
    est <- coef(model)
    
    temp_df <- data.frame(Variable = var,
                          Estimate = est[outcome], 
                          SE = se[outcome])
    
    
    temp_df$P.value <- 2 * pnorm(-abs(est[outcome]/se[outcome]))
    
    results_df <- rbind(results_df, temp_df)
  }
  
  results_df$P.value.adjusted <- p.adjust(results_df$P.value, method = "BH")
  results_df <- results_df[, c("Variable", "Estimate", "SE", "P.value", "P.value.adjusted")]
  # Combina os resultados com a base de dados class_data baseada na coluna 'Variable'
  results_df <- merge(results_df, class_data, by = "Variable", all.x = TRUE)
  
  # Salva os resultados em um arquivo Excel, com o nome baseado no desfecho analisado
  filename <- paste("GLM_ApoE4_non_carrier_DX_cat_CN_", outcome, "_adj.xlsx", sep = "")
  write.xlsx(results_df, file.path('your_location/results/DX/', filename))
}


#################################################################################

# Limpa o ambiente
rm(list = ls())
# Carrega as bibliotecas necessárias
library(lme4)
library(lmerTest)  # Fornece valores-p para modelos mistos
library(openxlsx)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(sandwich)
library(msm)
library(lmtest)

# Define o diretório de trabalho e carrega os dados
setwd('your_location/data/')
lipid_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Plasma')
phenotype_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Pheno')
class_data <- read.xlsx('class_plasma.xlsx')
class_data$Variable <- class_data$Lipids


# Combinando as bases de dados
lipid_data$RID <- as.double(lipid_data$RID)
phenotype_data$RID <- as.double(phenotype_data$RID)
combined_data <- phenotype_data %>% left_join(lipid_data, by = "RID")

combined_data2 <- subset(combined_data, ApoE4_cat == 0)
combined_data2 <- subset(combined_data2, DX == 1)

ls(phenotype_data)

# Lista de desfechos a serem analisados-------------------------------------- |||||||
outcomes <- c("ADAS.Cog13", "MMSE", "CDRSB")

#-------------------------------------- |||||||-------------------------------------- |||||||
for (outcome in outcomes) {
  results_df <- data.frame(Variable = character(0), 
                           Estimate = numeric(0), 
                           SE = numeric(0), 
                           P.value = numeric(0), 
                           P.value.adjusted = numeric(0))
  
  variables <- setdiff(colnames(lipid_data), "RID")
  
  for (var in variables) {
    formula_str <- as.formula(paste(var, "~", outcome, "+ Education + Age + SEX + BMI +  statin"))
    model <- glm(formula_str, family = "gaussian"(link = "identity"),  data = combined_data2)
    
    model_summary <- summary(model)
    se <- sqrt(diag(vcovHC(model, type = "HC0")))
    est <- coef(model)
    
    temp_df <- data.frame(Variable = var,
                          Estimate = est[outcome], 
                          SE = se[outcome])
    
    
    temp_df$P.value <- 2 * pnorm(-abs(est[outcome]/se[outcome]))
    
    results_df <- rbind(results_df, temp_df)
  }
  
  results_df$P.value.adjusted <- p.adjust(results_df$P.value, method = "BH")
  results_df <- results_df[, c("Variable", "Estimate", "SE", "P.value", "P.value.adjusted")]
  # Combina os resultados com a base de dados class_data baseada na coluna 'Variable'
  results_df <- merge(results_df, class_data, by = "Variable", all.x = TRUE)
  
  # Salva os resultados em um arquivo Excel, com o nome baseado no desfecho analisado
  filename <- paste("GLM_ApoE4_non_carrier_DX_cat_MCI_", outcome, "_adj.xlsx", sep = "")
  write.xlsx(results_df, file.path('your_location/results/DX/', filename))
}



#################################################################################

# Limpa o ambiente
rm(list = ls())
# Carrega as bibliotecas necessárias
library(lme4)
library(lmerTest)  # Fornece valores-p para modelos mistos
library(openxlsx)
library(broom.mixed)
library(dplyr)
library(ggplot2)
library(sandwich)
library(msm)
library(lmtest)

# Define o diretório de trabalho e carrega os dados
setwd('your_location/data/')
lipid_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Plasma')
phenotype_data <- read.xlsx('Filtered_Data.xlsx', sheet = 'Filtered_Pheno')
class_data <- read.xlsx('class_plasma.xlsx')
class_data$Variable <- class_data$Lipids


# Combinando as bases de dados
lipid_data$RID <- as.double(lipid_data$RID)
phenotype_data$RID <- as.double(phenotype_data$RID)
combined_data <- phenotype_data %>% left_join(lipid_data, by = "RID")

combined_data2 <- subset(combined_data, ApoE4_cat == 0)
combined_data2 <- subset(combined_data2, DX == 2)

ls(phenotype_data)

# Lista de desfechos a serem analisados-------------------------------------- |||||||
outcomes <- c("ADAS.Cog13", "MMSE", "CDRSB")

#-------------------------------------- |||||||-------------------------------------- |||||||
for (outcome in outcomes) {
  results_df <- data.frame(Variable = character(0), 
                           Estimate = numeric(0), 
                           SE = numeric(0), 
                           P.value = numeric(0), 
                           P.value.adjusted = numeric(0))
  
  variables <- setdiff(colnames(lipid_data), "RID")
  
  for (var in variables) {
    formula_str <- as.formula(paste(var, "~", outcome, "+ Education + Age + SEX + BMI +  statin"))
    model <- glm(formula_str, family = "gaussian"(link = "identity"),  data = combined_data2)
    
    model_summary <- summary(model)
    se <- sqrt(diag(vcovHC(model, type = "HC0")))
    est <- coef(model)
    
    temp_df <- data.frame(Variable = var,
                          Estimate = est[outcome], 
                          SE = se[outcome])
    
    
    temp_df$P.value <- 2 * pnorm(-abs(est[outcome]/se[outcome]))
    
    results_df <- rbind(results_df, temp_df)
  }
  
  results_df$P.value.adjusted <- p.adjust(results_df$P.value, method = "BH")
  results_df <- results_df[, c("Variable", "Estimate", "SE", "P.value", "P.value.adjusted")]
  # Combina os resultados com a base de dados class_data baseada na coluna 'Variable'
  results_df <- merge(results_df, class_data, by = "Variable", all.x = TRUE)
  
  # Salva os resultados em um arquivo Excel, com o nome baseado no desfecho analisado
  filename <- paste("GLM_ApoE4_non_carrier_DX_cat_AD_", outcome, "_adj.xlsx", sep = "")
  write.xlsx(results_df, file.path('your_location/results/DX/', filename))
}
