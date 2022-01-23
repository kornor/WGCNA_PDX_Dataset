
## SCript for analysing PDX data to look for expression modules as in the TCGA / Metabric set

#set for the PDX folder
setwd("~/Bioinformatics Work/Meth & RNA/PDX data");
# Load the package
library(WGCNA);
library(flashClust)
# String settings
options(stringsAsFactors = FALSE);

# load the metabric / tcga stuff

load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))

#datExpr1 = TCGA
#datExpre2 = Metabric

## Load in the PDX data
pdx_e <- read.table("ExpressionSamples_PDX.txt", 
                    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE) #PDX
pdx_e<- as.data.frame(pdx_e)

### preprocess the data using goodgenes (WCGNA)

gsg = goodSamplesGenes(pdx_e,verbose = 5);
gsg$allOK

## if return is true, all good to go

########## check that gene names match 

commongenes <- intersect(rownames(datExpr1), rownames(pdx_e))
pdx1 <- pdx_e[commongenes,]

#transpose again?
datExpr1_t <- as.data.frame(t(datExpr1))
pdx_t <- as.data.frame(t(pdx1))

######### Calculation of module preservation

setLabels = c("TCGA","PDX");
multiExpr = list(TCGA= list(data = datExpr1_t), PDX = list(data = pdx_t));
multiColor = list(TCGA = modules1);

### check for preservation

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
  referenceNetworks = 1,
  networkType="signed",
  nPermutations = 30,
  maxGoldModuleSize=100,
  maxModuleSize = 400,
  randomSeed = 1,
  quickCor = 0,
  verbose = 3)
} );
# Save the results
save(mp, file = "modulePreservation.RData");
load(file =  "modulePreservation.RData")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

##graphing
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold", "turquoise", "greenyellow","pink", "purple"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "All samples");
# Start the plot
sizeGrWindow(10, 5);
pdf(file ="Short_PDX-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust plotting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 3,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.5, cex.axis = 1.5, cex.main =1.8)
  #labelPoints(moduleSizes[plotMods], plotData[plotMods,p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2, lwd = 3)
    abline(h=10, col = "darkgreen", lty = 2, lwd = 3)
  }
}
# If plotting into a file, close it
dev.off();


#### now to calculate the module eigengene for each module as in TCGA:

## first need to trim so the genes in pdx_t are the same as in module colours (otherwise it will fail)

mod_gene <- as.data.frame(modules1)
row.names(mod_gene) <- row.names(datExpr1)

mod_pdx <- mod_gene[commongenes,]


###now the calcululation:

PCs_pdx <- moduleEigengenes(pdx_t, 
                 mod_pdx, 
                 impute = TRUE, 
                 nPC = 1, 
                 align = "along average", 
                 excludeGrey = FALSE, 
                 grey = if (is.numeric(mod_pdx)) 0 else "grey",
                 subHubs = TRUE,
                 softPower = 6,
                 scale = TRUE,
                 verbose = 0, indent = 0)

####works!!

ME_pdx    = PCs_pdx$eigengenes
distPC_pdx = 1-abs(cor(ME_pdx,use="p"))
distPC_pdx = ifelse(is.na(distPC_pdx), 0, distPC_pdx)
pcTree_pdx = hclust(as.dist(distPC_pdx),method="a") 
MDS_pdx  = cmdscale(as.dist(distPC_pdx),2)
colorsPDX = names(table(mod_pdx))



### add in the rownames

rownames(ME_pdx) <- colnames(pdx_e)

write.table(ME_pdx, file = "ME_PDX.txt", sep = "\t")
save(PCs_pdx, file = 'PDX_moduledata.Rdata')



##### MDS of the module expression in the PDX group:

par()
plot(MDS_pdx, col= colorsPDX,  main="MDS plot", cex=2, pch=19,
     xlab = "Principal component 1", ylab = "Principal component 2")


######### look at the samples clustering?

## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(pdx_t), method = "average");


# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)



#########################
load(file = "modulePreservation.RData")
## plot preservation etc with just the ones we like




#######################

#Module preservation without the initial tumour samples (in case they are overpowering)
library(tidyverse)
columns_to_remove <- grep("-T", colnames(pdx_e))
my_df <- pdx_e[,-columns_to_remove]
pdx_t <- as.data.frame(t(my_df))
######### Calculation of module preservation

setLabels = c("TCGA","PDX");
multiExpr = list(TCGA= list(data = datExpr1_t), PDX = list(data = pdx_t));
multiColor = list(TCGA = modules1);

### check for preservation

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          networkType="signed",
                          nPermutations = 30,
                          maxGoldModuleSize=100,
                          maxModuleSize = 400,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
# Save the results
save(mp, file = "modulePreservation_Xonly.RData");

load(file = "modulePreservation_Xonly.RData");

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

##graphing
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold", "turquoise", "greenyellow","pink", "purple"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Xenograft samples only");
# Start the plot
sizeGrWindow(10, 5);
pdf(file ="Short_trimmed_PDX-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust plotting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 3,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.5, cex.axis = 1.5, cex.main =1.8)
  #labelPoints(moduleSizes[plotMods], plotData[plotMods,p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2, lwd = 3)
    abline(h=10, col = "darkgreen", lty = 2, lwd = 3)
  }
}
# If plotting into a file, close it
dev.off();

####################################
#alternative plotting


modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize

modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold","turquoise", "greenyellow","pink", "purple"))
# Text labels for points
point.label = modColors[selectModules]

#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres

medianRank = mp$preservation$observed[[ref]][[test]][, 2]
Zsummary =  mp$preservation$Z[[ref]][[test]][, 2]

par(mfrow=c(1,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSizes[selectModules],medianRank[selectModules],col=1,
     bg=modColors[selectModules],pch = 21,main="medianRank Preservation",
     cex = 4, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSizes[selectModules],medianRank[selectModules],point.label,cex=1.5,offs=0.09)

# plot Zsummary versus module size
plot(moduleSizes[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,main="Zsummary Preservation",
     cex=4,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSizes[selectModules],Zsummary[selectModules],point.label,cex=1.5, offs = 0.1)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2); abline(h=10, col = "red", lty = 3)
