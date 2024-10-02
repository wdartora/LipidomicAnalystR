
#Clear all
rm(list = ls())
# Load the WGCNA package
library(WGCNA)
library(ggplot2)
setwd("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/data")
# Load the data saved in the first part
lnames = load(file = "03_output_Input_plasma.RData")

#>
#=====================================================================================
# set a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 5, to=25, by=2))
#threshold power = 1 a 10, 
#TAOS increased by = 2 a 20

# Calculate a soft-thresholding power used for downstream analysis
#No banco da Matriz gene (Lipid)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft[["powerEstimate"]]
#13 =  poder estimado (poderemos verificar se faz sentido, olhando os graficos abaixo)


# Plot Scale-free topology fit index as a function of the soft-thresholding power
# Topology plot (X=Power vs Y=Scale free topology)
p<- ggplot(sft$fitIndices, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

p
ggsave("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/03-Topology_plot_Plasma.png",
       p,  width = 8, height = 8, dpi = 300)
#(X=Power vs Y=Mean connectivity
# Plot Mean connectivity as a function of the soft-thresholding power
d<- ggplot(sft$fitIndices, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
d

ggsave("C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/04-Mean_connectivity_plot_Plasma.png",
       d,  width = 8, height = 8, dpi = 300)

#=====================================================================================
# use soft thresholding power 7 and calculate the adjacencies
adjacency = adjacency(datExpr, power = 6)
# considerei o plot = 04-Mean_connectivity_plot_Plasma
#=====================================================================================
# To minimize effects of noise and spurious associations, we transform the adjacency 
# into Topological Overlap Matrix (TOM)
TOM = TOMsimilarity(adjacency)
#..matrix multiplication (system BLAS)..
#..normalization..
# then calculate the corresponding dissimilarity
dissTOM = 1-TOM

#=====================================================================================
# use the hierarchical clustering function to produce a hierarchical clustering 
# gene tree (dendrogram)
geneTree = hclust(as.dist(dissTOM), method = "average")

# 
#=====================================================================================
# set the minimum lipid in a module (module size)
minModuleSize = 30
#minModuleSize = 40
#Pequeno numero, cria mais Models,
#Grande numero, cria menos Models
# Module identification using dynamic tree cut:

#Identificacao dos models
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, 
                            pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

#cutHeight=0.99,
#deepSplit = TRUE,
#Tabela para ver quantos models temos e o numero de genes (Lipids) em cada model
table(dynamicMods)
#9 models no total 
#dynamicMods
#   0   1   2   3   4   5   6   7   8   9  10  11
#   173  96  67  66  65  46  46  42  40  39  36  33
#=====================================================================================
# plot the dynamic modules under the gene tree
dynamicMods
#modulos NUMERICOS


# first Convert numeric labels (Muda de numero para cores os models) into colors
dynamicColors = labels2colors(dynamicMods)
dynamicColors
table(dynamicColors)
#9 models (cores)
#dynamicColors
#      black        blue       brown       green greenyellow        grey
#         42          67          66          46          33         173
#    magenta        pink      purple         red   turquoise      yellow
#         39          40          36          46          96          65

# Plot the dendrogram and colors underneath

sizeGrWindow(8,6)

png(file="C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/06-dendrogram_moduleColor_Plasma.png",
    width=10, height=10, units="in", res=300)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Lipid dendrogram and module colors - Plasma")
dev.off()

#No grafico e possivel identificar que:
# blue color, brown color, e lightcyan, representam os modulos com mais genes (lipids)

#=====================================================================================
#Podemos identificar os modulos com maiores similiaridades na co-expression genes (lipid co-expression)
# Calculate Module Eigengenes (ME) and cluster them on their correlation to quantify 
# co-expression similarity of entire modules,
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result from Eigen Genes (eigen lipids)
sizeGrWindow(7, 6)

png(file="C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/07-eigen_cluster_moduleColor2_Plasma.png",
    width=10, height=10, units="in", res=300)
plot(METree, main = "Clustering of module eigensLipids",
     xlab = "", sub = "")
MEDissThres = 0.4
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

#Deveria ser 0.25

# Call an automatic merging function [ para considerar agregar os modulos de baixo e de cima do 
# ponto de corte]
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# Assim fico com menos modulos que acima, fazendo com que os modulos altamente correlacionados,
# estejam juntos (mais que 0.75).

# The merged module colors
mergedColors = merge$colors
table (mergedColors)


# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#=====================================================================================
sizeGrWindow(12, 9)
png(file="C:/Users/wjd4002/Documents/William/Project/ADNI/ADNI_wjd/Proj12/plots/06-dendrogram_moduleColor_Serumv2.png",
    width=10, height=10, units="in", res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Agora podemos observar o Dynamic Tree e o Merged dynamic, onde podemos ver que alguns modulo
# se combinaram (merged).
dev.off()
#=====================================================================================
# Save module colors and labels for use in subsequent parts
# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file = "04_output_networkConstruction_Plasma.RData")


































