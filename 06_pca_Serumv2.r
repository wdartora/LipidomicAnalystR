# Limpar o ambiente
rm(list = ls())

# Carregar pacotes necessários
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(umap)
library(factoextra)
library(Seurat)

# Definir o diretório
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data")

# Carregar os dados
df_BB <- read.xlsx("Filtered_Data.xlsx", sheet = 'Filtered_Serum')
df_CC <- read.xlsx("Filtered_Data.xlsx", sheet = 'Filtered_Pheno')
lipid_classes <- read.xlsx('Class_Serum_list.xlsx')

# Garantir que as colunas RID sejam caracteres para junção
df_BB$RID <- as.character(df_BB$RID)
df_CC$RID <- as.character(df_CC$RID)

# Adicionar a coluna DX a partir de Diagnosis em df_CC, combinando todas as categorias MCI em uma única categoria "MCI"
df_CC$DX <- ifelse(df_CC$Diagnosis == "AD", "AD", 
                   ifelse(df_CC$Diagnosis == "CN", "CN", "MCI"))

# Converter ApoE4_cat em rótulos "non carrier" e "carrier"
df_CC$ApoE4_cat <- ifelse(df_CC$ApoE4_cat == 1, "carrier", "non carrier")

# Verificar os valores únicos na coluna DX e ApoE4_cat para garantir que não há NA indesejados
print(unique(df_CC$DX))
print(unique(df_CC$ApoE4_cat))

# Separar os dados e o fenótipo
data_BB <- df_BB %>%
  select(-RID) %>%
  mutate(across(everything(), as.numeric))

phenotype <- df_CC %>%
  select(RID, DX, ApoE4_cat)

# Verificar se há NA introduzidos durante a conversão e removê-los
data_BB <- na.omit(data_BB)

# Realizar UMAP
set.seed(123) # Para reprodutibilidade
umap_result <- umap(data_BB)

# Adicionar resultados do UMAP ao fenótipo
umap_df <- as.data.frame(umap_result$layout)
umap_df$RID <- df_BB$RID
umap_df <- umap_df %>%
  left_join(phenotype, by = "RID")

# Verificar se DX tem NA após o merge
print(unique(umap_df$DX))

# Definir o diretório de saída
setwd('C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/')

# Escolher o número ideal de clusters usando o método do cotovelo com factoextra
elbow_plot <- fviz_nbclust(data_BB, kmeans, method = "wss") +
  labs(title = "Elbow Method for Choosing Number of Clusters", x = "Number of Clusters (k)", y = "Total Within Sum of Squares") +
  theme_light()

# Definir o número ideal de clusters
optimal_clusters <- 6 # Defina isso com base no gráfico gerado

# Adicionar uma linha vertical indicando o número ideal de clusters
elbow_plot <- elbow_plot + geom_vline(xintercept = optimal_clusters, linetype = "dashed", color = "red")
#ggsave("10_A_Serum_elbow.png", elbow_plot, height = 6, width = 10, dpi = 300)

# Realizar clustering usando K-means
set.seed(123) # Para reprodutibilidade
kmeans_result <- kmeans(data_BB, centers = optimal_clusters)

# Adicionar os clusters ao DataFrame do UMAP
umap_df$Cluster <- as.factor(kmeans_result$cluster)

# Plot UMAP colorido pelos clusters
up <- ggplot(umap_df, aes(V1, V2, color = Cluster)) +
  geom_point(size = 3, shape = 19, alpha = 0.85, stroke = 1) +
  labs(title = "UMAP of Serum Data (K-means Clusters)", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  theme_light()
up
#ggsave("10_UMAP_Serum.png", up, height = 6, width = 10, dpi = 300)

# Identificar o cluster de interesse (por exemplo, Cluster 3)
cluster_of_interest <- 3

# Filtrar os dados do cluster de interesse
cluster_data <- umap_df %>% filter(Cluster == cluster_of_interest)

# Análise descritiva das características do cluster de interesse
summary(cluster_data)

# Ordenar DX para garantir que CN, MCI e AD estejam na ordem correta
cluster_data$DX <- factor(cluster_data$DX, levels = c("CN", "MCI", "AD"))

# Visualizar características fenotípicas do cluster de interesse
p <- ggplot(cluster_data, aes(x = DX, fill = as.factor(ApoE4_cat))) +
  geom_bar(position = "dodge") +
  labs(title = paste("Distribution of Diagnosis and ApoE4 Status in Cluster", cluster_of_interest), 
       x = "Diagnosis", 
       y = "Count", 
       fill = "ApoE4 Status") +
  scale_fill_brewer(palette = "Set1") +
  theme_light()
p
#ggsave("10_UMAP_Serum_cluster_1.png", p, height = 6, width = 10, dpi = 300)

# Identificar os RIDs presentes no cluster de interesse
rids_in_cluster <- cluster_data$RID

# Filtrar os dados fenotípicos correspondentes a esses RIDs
phenotype_in_cluster <- phenotype %>%
  filter(RID %in% rids_in_cluster)

# Calcular a abundância total dos lipídios por grupo, diagnóstico e status ApoE4
abundance_data <- df_BB %>%
  filter(RID %in% rids_in_cluster) %>%
  gather(key = "Lipids", value = "Abundance", -RID) %>%
  left_join(lipid_classes, by = "Lipids") %>%
  left_join(phenotype_in_cluster, by = "RID") %>%
  group_by(RID, Group, DX, ApoE4_cat) %>%
  summarise(Total_Abundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

# Replace  (COL) 
colnames(abundance_data)[which(names(abundance_data)=="Group")]<-"Class"

#RID numeric
abundance_data$RID <- as.numeric(abundance_data$RID)
# Salvar os dados de abundância total em um arquivo Excel
write.xlsx(abundance_data, "Lipid_Abundance_in_Cluster_3_Grouped_Serum.xlsx")
head(abundance_data)
# Exibir as categorias de DX e ApoE4_cat presentes no cluster de interesse
#print(phenotype_in_cluster)

# Salvar as anotações dos fenotipos em um arquivo
#write.xlsx(phenotype_in_cluster, "Phenotype_in_Serum_Cluster_3.xlsx")
