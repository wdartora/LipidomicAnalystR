#Clear all
rm(list = ls())
# Load the WGCNA package
library(WGCNA)
setwd("your_location/data")
# Load the expression and trait data saved in the first part
lnames = load(file = "03_output_Input_plasma.RData")
#Matriz de genes (lipids)[datExpr] + Phenotype data [datTraits]

# Load network data saved in the second part.
lnames = load(file = "04_output_networkConstruction_Plasma.RData")
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
png(file="your_location/plots/08-Heatmap_moduleColor_Plasma.png",
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
               main = paste("Module-trait relationships - Plasma "))
# Aqui no plot (heatmap), sera mostrado o valor de correlacao e entre parenteses o p-value, entre
# os modulos e os Phenotype

dev.off()
#Correlaction between each model with each traits



#########################################################################################
#================================================================================

#########################################################################################
#================================================================================

#Abeta42_A40_csf
#### 2. Gene relationship to trait and important modules==================================
# Define variable Abeta42_A40_csf containing the Abeta42_A40_csf column of datTrait
#Precisamos retirar as informacoes do tratamento que consideramos mais significativo, observando
#o heatmap
Abeta42_A40_csf = as.data.frame(datTraits$Abeta42_A40_csf)
names(Abeta42_A40_csf) = "Abeta42_A40_csf"
#Apenas renomeando a coluna

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, Abeta42_A40_csf, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

names(geneTraitSignificance) = paste("GS.", names(Abeta42_A40_csf), sep="")
#Calcular o pvalue para os Genes Significants (Lipids Significannts)
names(GSPvalue) = paste("p.GS.", names(Abeta42_A40_csf), sep="")


######________________________________________________________________________________________>>>
# 3. Intramodular analysis: identifying genes with high GS and MM =========================
sizeGrWindow(7, 7)
png(file="your_location/plots/09-scatter_grey_Ab42_40_csf_Plasma.png",
    width=18, height=10, units="in", res=300)

par(mfrow = c(1,1))

# Plot inicial sem borda
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in brown module"),
                   ylab = "Lipid significance for ADAS.Cog13 (csf)",
                   main = paste("Module membership vs. Lipid significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)

# Adicionar borda preta aos pontos sem preencher o interior
points(abs(geneModuleMembership[moduleGenes, column]),
       abs(geneTraitSignificance[moduleGenes, 1]),
       col = "darkgrey", pch = 16, cex = 1.5, lwd = 1.5)

dev.off()






#################################################################################################

#Abeta40_csf
#### 2. Gene relationship to trait and important modules==================================
# Define variable Abeta40_csf containing the Abeta40_csf column of datTrait
#Precisamos retirar as informacoes do tratamento que consideramos mais significativo, observando
#o heatmap
Abeta40_csf = as.data.frame(datTraits$Abeta40_csf)
names(Abeta40_csf) = "Abeta40_csf"
#Apenas renomeando a coluna

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, Abeta40_csf, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

names(geneTraitSignificance) = paste("GS.", names(Abeta40_csf), sep="")
#Calcular o pvalue para os Genes Significants (Lipids Significannts)
names(GSPvalue) = paste("p.GS.", names(Abeta40_csf), sep="")


# Todaessa etapa serve para concatenar os dados entre p-value, genes e phenotype data




# 4. Save gene Information===================================================

nomes_df<- data.frame(datExpr)
names(nomes_df)
#=====================================================================================
blue_ALL<-names(nomes_df)[moduleColors=="blue"]
blue_ALL<-data.frame(blue_ALL)
#=====================================================================================
write.xlsx(blue_ALL, 
           "your_location/results/blue_ALL.xlsx")




# Todaessa etapa serve para concatenar os dados entre p-value, genes e phenotype data



######________________________________________________________________________________________>>>
# 4. Save ALL Modules ===================================================

nomes_df<- data.frame(datExpr)
names(nomes_df)
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="greenyellow"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/greenyellow_Plasma.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="black"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/black_Plasma.xlsx")

#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="magenta"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/magenta_Plasma.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="purple"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/purple_Plasma.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="turquoise"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/turquoise_Plasma.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="green"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/green_Plasma.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="blue"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/blue_Plasma.xlsx")
#=====================================================================================

#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="red"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/red_Plasma.xlsx")

#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="brown"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/brown_Plasma.xlsx")
#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="pink"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/pink_Plasma.xlsx")

#=====================================================================================
module_ALL<-names(nomes_df)[moduleColors=="grey"]
module_ALL<-data.frame(module_ALL)
#=====================================================================================
write.xlsx(module_ALL, 
           "your_location/results/grey_Plasma.xlsx")


