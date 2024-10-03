#Clear all
rm(list = ls())
# Load the WGCNA package
library(WGCNA)
setwd("your_location/data")
# Load the expression and trait data saved in the first part
lnames = load(file = "03_output_Input_serum.RData")
#Matriz de genes (lipids)[datExpr] + Phenotype data [datTraits]

# Load network data saved in the second part.
lnames = load(file = "04_output_networkConstruction_Serum.RData")
# Hierarchical gene (Lipids) tree[geneTree] + 19 Modulos [MEs]

#================================================================================


### 1. Quantifying module trait associations================================================
# Recalculate MEs with color labels [ entre os modulos e o Phenotype]
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# Correlations values

# Quantify module-trait correlations [ entre os modulos e o Phenotype]
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr))
# p-value from correlations

# Display the correlation values within a heatmap plot
#convertendo o data em text matrix, para criar o Heatmap e poder identificar melhor a correlacao
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(", 
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)


#setando parametros do grafico
sizeGrWindow(20,16)
png(file="your_location/plots/08-Heatmap_moduleColor_Serum.png",
    width=18, height=10, units="in", res=300)
par(mar = c(10, 8, 5, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               xLabelsAdj = 0.9,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               zlim = c(-0.6,0.6),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               main = paste("Module-trait relationships - Serum "))
# Aqui no plot (heatmap), sera mostrado o valor de correlacao e entre parenteses o p-value, entre
# os modulos e os Phenotype

dev.off()
#Correlaction between each model with each traits



#########################################################################################
#================================================================================

#Abeta40_csf
#### 2. Gene relationship to trait and important modules==================================
# Define variable Abeta40_csf containing the Abeta40_csf column of datTrait
#Precisamos retirar as informacoes do tratamento que consideramos mais significativo, observando
#o heatmap
TAU_Abeta42 = as.data.frame(datTraits$TAU_Abeta42)
names(TAU_Abeta42) = "Abeta40_csf"
#Apenas renomeando a coluna

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, TAU_Abeta42, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

names(geneTraitSignificance) = paste("GS.", names(TAU_Abeta42), sep="")
#Calcular o pvalue para os Genes Significants (Lipids Significannts)
names(GSPvalue) = paste("p.GS.", names(TAU_Abeta42), sep="")


# Todaessa etapa serve para concatenar os dados entre p-value, genes e phenotype data



######________________________________________________________________________________________>>>
# 3. Intramodular analysis: identifying genes with high GS and MM =========================
module = "brown"
column = match(module, modNames)
moduleGenes = moduleColors==module

#> Iremos plotar o genes (lipids) significance VS o modulo membership
sizeGrWindow(7, 7)
png(file="your_location/plots/09-scatter_brown_TAU_Abeta42_Serum.png",
    width=18, height=10, units="in", res=300)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in brown module"),
                   ylab = "Lipid significance for TAU_Abeta42",
                   main = paste("Module membership vs. Lipid significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# 4. Save gene Information===================================================


######________________________________________________________________________________________>>>
# 4. Save ALL Modules ===================================================

nomes_df<- data.frame(datExpr)
names(nomes_df)


#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="turquoise"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/turquoise_Serum.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="green"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/green_Serum.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="blue"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/blue_Serum.xlsx")

#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="black"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/black_Serum.xlsx")

#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="grey"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/grey_Serum.xlsx")


