#Run:
#script 0, 1 and 2 (without ) - rm(list = ls())

#Exporting the network to a cytoscape format
#Recalculating topological overlap, if necessary
#TOM = TOMsimilarityFromExpr(expression0, power = 10);
#Select the modules
#modules = c("brown", "red"); #chose modules that u want to export
#Select the gene modules
genes = colnames(datExpr)

#if you want export specific colors, substitute the second modulecolors by above modules
inModule = is.finite(match(moduleColors, moduleColors))
modGenes = genes[inModule]

#Select the corresponding topologic overlap 
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)
modTOMSignificantes = which(modTOM>0.4)
#####warnings()

#Organize the genes by importance inside the module
genes = colnames(datExpr)
#sum(is.na(genes))
#It must return 0.


#Create the dataframe since the beginning
geneInfo0 = data.frame(ESTs = genes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

#Order the modules by the significance by a character Ex: 
modOrder = order(-abs(cor(MEs, Hippvol, use = "p")))

#Add information of the members of the module in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

#Order the genes of geneinfo variable first by the color of the module, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Hippvol))
geneInfo = geneInfo0[geneOrder, ]

setwd('your_location/results/')
#write the file with the ordered values
write.csv(geneInfo, file = "geneInfo_04.csv")

#Export the network in list files os n edges that cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeEdgeFilev2.txt",
                               nodeFile = "CytoscapeNodeFilev2.txt",
                               weighted = TRUE,
                               threshold = 0.05,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
