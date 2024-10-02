################################################################################################
# 1) Plasma - Volcano plot liso
# Clear all
rm(list = ls())
library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggrepel)

# Caminho dos arquivos
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/')

# Carregar os dados dos arquivos Excel
plasma_aa <- read.xlsx("Filtered_Data.xlsx", sheet = "wgcna.Plasma")
metadata <- read.xlsx("Filtered_Data.xlsx", sheet = "Filtered_Pheno")
metadata$Sample <- metadata$Lipids
metadata <- metadata[, c('Sample', 'ApoE4_cat')]

# Transformar dados de Plasma_AA para formato longo
plasma_aa_long <- plasma_aa %>%
  pivot_longer(-Lipids, names_to = "Sample", values_to = "Concentration")

# Unir os dados com Metadata
merged_data <- merge(plasma_aa_long, metadata, by.x = "Sample", by.y = "Sample")

# Calcular log2fold change e p-values
results <- merged_data %>%
  group_by(Lipids) %>%
  summarise(log2FC = log2(mean(Concentration[ApoE4_cat == 1]) / mean(Concentration[ApoE4_cat == 0])),
            p_value = t.test(Concentration ~ ApoE4_cat)$p.value)

# Definir cores para pontos significativos e não significativos
results <- results %>%
  mutate(fill = case_when(
    p_value < 0.05 & abs(log2FC) > 1 ~ "red",
    p_value < 0.05 & abs(log2FC) <= 1 ~ "blue",
    TRUE ~ "gray"
  ))

# Filtrar os lipídios para rotulagem
label_data <- results %>%
  filter(p_value < 0.05 & (log2FC > 3 | log2FC < -3))

# Plotar Volcano Plot
volcano_plot <- ggplot(results, aes(x = log2FC, y = -log10(p_value), fill = fill)) +
  geom_point(shape = 21, color = "black", alpha = 0.85, stroke = 1, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  scale_fill_identity() +
  geom_label_repel(data = label_data, aes(label = Lipids), size = 3, max.overlaps = 10, box.padding = 0.5, fill = "white") +
  annotate("label", x = -3.5, y = Inf, label = "Non-carrier", vjust = 2, size = 4, color = "black", fill = "#e9e9e9", label.padding = unit(0.5, "lines"), label.size = 0.25) +
  annotate("label", x = 3.5, y = Inf, label = "Carrier", vjust = 2, size = 4, color = "black", fill = "#e9e9e9", label.padding = unit(0.5, "lines"), label.size = 0.25) +
  labs(title = "ApoE4 - Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14))

# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("15-Plasma_volcano_plot1.png", plot = volcano_plot, width = 6, height = 8)


###########################################################################################
################################################################################################
# 1) Plasma - Volcano plot liso
# Clear all
rm(list = ls())
library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggrepel)

# Caminho dos arquivos
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/')

# Carregar os dados dos arquivos Excel
plasma_aa <- read.xlsx("Filtered_Data.xlsx", sheet = "wgcna.Plasma")
metadata <- read.xlsx("Filtered_Data.xlsx", sheet = "Filtered_Pheno")
metadata$Sample <- metadata$Lipids

#Create apoe4 char and DX char
#metadata$ApoE4 <- ifelse(metadata$ApoE4 == 1, 'carrier', 'non-carrier')
metadata$DX_text <- ifelse(metadata$DX == 0, 'CN',
                ifelse(metadata$DX ==  1, 'MCI','AD'))

# Keep var                
metadata <- metadata[, c('Sample', 'DX', 'DX_text')]

#sub-set (Comparision between AD vs CN)
metadata <- metadata %>%
    filter(DX_text != "MCI")

# Transformar dados de Plasma_AA para formato longo
plasma_aa_long <- plasma_aa %>%
  pivot_longer(-Lipids, names_to = "Sample", values_to = "Concentration")

# Unir os dados com Metadata
merged_data <- merge(plasma_aa_long, metadata, by.x = "Sample", by.y = "Sample")

# Calcular log2fold change e p-values
results <- merged_data %>%
  group_by(Lipids) %>%
  summarise(log2FC = log2(mean(Concentration[DX == 2]) / mean(Concentration[DX == 0])),
            p_value = t.test(Concentration ~ DX)$p.value)

# Definir cores para pontos significativos e não significativos
results <- results %>%
  mutate(fill = case_when(
    p_value < 0.05 & abs(log2FC) > 1 ~ "red",
    p_value < 0.05 & abs(log2FC) <= 1 ~ "blue",
    TRUE ~ "gray"
  ))

# Filtrar os lipídios para rotulagem
label_data <- results %>%
  filter(p_value < 0.05 & (log2FC > 3 | log2FC < -3))

# Plotar Volcano Plot
volcano_plot <- ggplot(results, aes(x = log2FC, y = -log10(p_value), fill = fill)) +
  geom_point(shape = 21, color = "black", alpha = 0.85, stroke = 1, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  scale_fill_identity() +
  geom_label_repel(data = label_data, aes(label = Lipids), size = 3, max.overlaps = 10, box.padding = 0.5, fill = "white") +
  annotate("label", x = -3.5, y = Inf, label = "CN", vjust = 2, size = 4, color = "black", fill = "#e9e9e9", label.padding = unit(0.5, "lines"), label.size = 0.25) +
  annotate("label", x = 3.5, y = Inf, label = "AD", vjust = 2, size = 4, color = "black", fill = "#e9e9e9", label.padding = unit(0.5, "lines"), label.size = 0.25) +
  labs(title = "Diagnosis - Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14))

# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("15-Plasma_volcano_DX_ADvsCN.png", plot = volcano_plot, width = 6, height = 8)



#--------------------------------------->>>>
# Comparision betweem MCI vs AD
# Caminho dos arquivos
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data/')

# Carregar os dados dos arquivos Excel
metadata <- read.xlsx("Filtered_Data.xlsx", sheet = "Filtered_Pheno")
metadata$Sample <- metadata$Lipids

#Create apoe4 char and DX char
#metadata$ApoE4 <- ifelse(metadata$ApoE4 == 1, 'carrier', 'non-carrier')
metadata$DX_text <- ifelse(metadata$DX == 0, 'CN',
                ifelse(metadata$DX ==  1, 'MCI','AD'))

# Keep var                
metadata <- metadata[, c('Sample', 'DX', 'DX_text')]

#sub-set (Comparision between AD vs CN)
metadata <- metadata %>%
    filter(DX_text != "AD")

# Transformar dados de Plasma_AA para formato longo
plasma_aa_long <- plasma_aa %>%
  pivot_longer(-Lipids, names_to = "Sample", values_to = "Concentration")

# Unir os dados com Metadata
merged_data <- merge(plasma_aa_long, metadata, by.x = "Sample", by.y = "Sample")

# Calcular log2fold change e p-values
results <- merged_data %>%
  group_by(Lipids) %>%
  summarise(log2FC = log2(mean(Concentration[DX == 1]) / mean(Concentration[DX == 0])),
            p_value = t.test(Concentration ~ DX)$p.value)

# Definir cores para pontos significativos e não significativos
results <- results %>%
  mutate(fill = case_when(
    p_value < 0.05 & abs(log2FC) > 1 ~ "red",
    p_value < 0.05 & abs(log2FC) <= 1 ~ "blue",
    TRUE ~ "gray"
  ))

# Filtrar os lipídios para rotulagem
label_data <- results %>%
  filter(p_value < 0.05 & (log2FC > 3 | log2FC < -3))

# Plotar Volcano Plot
volcano_plot <- ggplot(results, aes(x = log2FC, y = -log10(p_value), fill = fill)) +
  geom_point(shape = 21, color = "black", alpha = 0.85, stroke = 1, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  scale_fill_identity() +
  geom_label_repel(data = label_data, aes(label = Lipids), size = 3, max.overlaps = 10, box.padding = 0.5, fill = "white") +
  annotate("label", x = -3.5, y = Inf, label = "CN", vjust = 2, size = 4, color = "black", fill = "#e9e9e9", label.padding = unit(0.5, "lines"), label.size = 0.25) +
  annotate("label", x = 3.5, y = Inf, label = "MCI", vjust = 2, size = 4, color = "black", fill = "#e9e9e9", label.padding = unit(0.5, "lines"), label.size = 0.25) +
  labs(title = "Diagnosis - Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14))

# Salvar o gráfico
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')
ggsave("15-Plasma_volcano_DX_MCIvsCN.png", plot = volcano_plot, width = 6, height = 8)


