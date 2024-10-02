################################################################################

# S e r u m  ------------------------------------------------------------>>>>>>

################################################################################


rm(list = ls())

#set the folder
setwd("your_location/data")

#Loading the packages
library(WGCNA)
library(tidyverse)
library(ggplot2)
library(openxlsx)


# loading the female liver data set
datExpr = read.xlsx("Filtered_Data.xlsx", sheet = 'wgcna.Serum')
#col 1 =  (lipids)
#col 2... = comecam os participantes

#abaixo eu transformo a coluna 1 em nomes de linhas
rownames (datExpr) = datExpr$Lipids

#Removo as primeiras colunas, e mantenho apenas os samples
datExpr = datExpr[, -1]


#Transpose o data, pra que sejam os lipidios nas colunas, e o sample no nome das linhas
datExpr = t(datExpr)

#=====================================================================================
# Checking data for missing values

#Flagging lipids and samples with too many missing values..
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
#[1] TRUE
#=====================================================================================

#Funcao para remover os missings absolutos
if (!gsg$allOK)
{ # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing lipids:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

#=====================================================================================
# cluster the samples to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr), method = "average")

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
#setando alguns parametros 
par(cex = 0.6);
par(mar = c(0,4,2,0))

png(file="your_location/plots/01-cluster_outlier_Serum.png",
   width=10, height=10, units="in", res=300)
#Hierarchical Plot para observar outliers 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()

#apenas printando os nomes das linhas (sao os samples os nomes)
rownames(datExpr)

#=====================================================================================

#Phenotype data -----//
datTraits = read.xlsx("Filtered_Data.xlsx", sheet = 'Filtered_Pheno')
#datTraits = read.csv("02_input_ClinicalTraits.csv")


#Criando os estratos de analise:
# Definindo primeiramente os dados de analise para ApoE4 carrier. 

#ApoE4 carrier
#datTraits<- datTraits %>% 
#  filter(ApoE4 == 1 )

print(colnames(datTraits))
# remove columns that hold information we do not need.
datTraits = datTraits[, c(3:6, 9:11,  21:31)]
# Remove a coluna "composite"
#datTraits <- datTraits[, !names(datTraits) %in% "composite"]

#Re-order
ls(datTraits)
datTraits<-datTraits %>% 
  relocate(any_of(c("Lipids",
                    "ADAS.Cog13","CDRSB","MMSE",
                    "LeftHippocampus","RightHippocampus","Hippvol",
                    "cingulate","frontal","parietal","temporal",
                    "TAU","PTAU","TAU_Abeta42","Abeta40_csf","Abeta42_csf","Abeta42_A40_csf")))

#renomear linhas com os samples
datTraits <- column_to_rownames(datTraits, var = 'Lipids')

# Formar um quadro de dados semelhante aos dados de expressão que conterá as características clínicas.
rownames(datExpr)
rownames(datTraits)
#Observamos que existem alguns Samples presentes no Phenotype data, que nao possuem dados de genes
all(rownames(datExpr) == rownames (datTraits))
#[1] FALSE

#Aqui eu adiciono esse trecho para nao ter erros na frente. Preciso manter na matrix,
#os participantes que se encontram no Phenotype data.
#> ============
#> #abaixo eu transformo a coluna 1 em nomes de linhas
#rownames (datTraits) = datTraits$Lipid
#Create the list with ApoE4 carriers
nomes_sample <- rownames(datTraits)
#Create a subset of lipids matrix with Phenotype data samples
datExpr <-subset(datExpr, rownames(datExpr) %in% nomes_sample)
#> ============


# =============================================================>>>>>>>>>>>
#Aqui iremos observar quem sao os samples que nao possuem dados em comum entre as bases de dados.
all(rownames(datExpr) %in% rownames(datTraits))
#TRUE = Significa que todos os samples da base Lipids, estao no phenotype data.
all(rownames(datTraits) %in% rownames(datExpr))
#FALSE = Significa que Nao sao todos os samples da base phenotype data, estao no Genes (Lipids).
#Manter na base Phenotype apenas os samples da base de Genes (Lipids)
datTraits = datTraits[rownames(datExpr), ];
#1407 observacoes
# =============================================================>>>>>>>>>>>


#Novamente conferimos
all(rownames(datTraits) %in% rownames(datExpr))
#TRUE = Significa que todos os samples da base phenotype data, estao no Genes (Lipids).
all(rownames(datTraits) == rownames (datExpr))
#TRUE = Confirmado.

#=====================================================================================

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);


# Plot the sample dendrogram and the colors underneath.
png(file="your_location/plots/02-dendrogram_and_trait_Serum.png",
    width=10, height=10, units="in", res=300)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
#Podemos observar que temos alguns missings (cinza) para alguns phenotype

#=====================================================================================
#Save o database
save(datExpr, datTraits, file = "03_output_Input_serum.RData")




