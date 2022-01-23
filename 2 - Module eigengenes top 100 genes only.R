

##### Module eigengenes of top 100 genes only
# for TCGA, METABRIC + PDX

### Need to gather
## list of top 100 genes in "both"
## create exp file for the 3 datasets just using the top 100
# create colours file for just the top 100


#set for the PDX folder
setwd("~/Bioinformatics Work/Meth & RNA/PDX data");
# Load the package
library(WGCNA);
library(flashClust)
library(dplyr)
library(plyr)
library(tidyr)
# String settings
options(stringsAsFactors = FALSE);


load(file = "MetaAnalysis_trimmed_input.RData")
load(file = "Modules_DS0.RData")
colorsA1 = names(table(modules1))
#datExpr1 = TCGA
#datExpre2 = Metabric


## Load in the PDX data
pdx_e <- read.table("ExpressionSamples.txt", 
                    sep = "\t", header = TRUE, row.names = 1, check.names = FALSE) #PDX
pdx_e<- as.data.frame(pdx_e)


##load in the top 100 genes:
# Note that the way the top 100 are calculated means that some genes may appear in multiple lists 
# because the modules are small and we have chosen the top 100, for eg in module tan, which as 33 genes, 
# genes from other modules are coming up in the top 100 genes in that module!
top100 <- read.table("Top100_all_modules_both_networks_BYG.txt", sep = "\t", header = TRUE, row.names = 1)

##Trim the data expression files by the genes in the top100

commongenes <- intersect(rownames(datExpr1), rownames(top100))
tcga100 <- datExpr1[commongenes,]
metabric100 <- datExpr2[commongenes,]
pdx100 <- pdx_e[commongenes,]
mod100 <- top100[commongenes,]

##need the genes as columns for this calculation
tcga100 <- as.data.frame(t(tcga100))  


###
tcgaME <- moduleEigengenes(tcga100, 
                            mod100, 
                            impute = TRUE, 
                            nPC = 1, 
                            align = "along average", 
                            subHubs = TRUE,
                            softPower = 6,
                            scale = TRUE,
                            verbose = 0, indent = 0)

##works


ME_tcga   = tcgaME$eigengenes
distPC_tcga = 1-abs(cor(ME_tcga,use="p"))
distPC_tcga = ifelse(is.na(distPC_tcga), 0, distPC_tcga)
pcTree_tcga = hclust(as.dist(distPC_tcga),method="a") 
MDS_tcga  = cmdscale(as.dist(distPC_tcga),2)
colorstcga100 = names(table(mod_pdx))



### add in the rownames

rownames(ME_tcga) <- rownames(tcga100)

write.table(ME_tcga, file = "ME_tcag_100.txt", sep = "\t")
save(tcgaME, file = 'TCGA100_moduledata.Rdata')




###pdx
pdx100 <- as.data.frame(t(pdx100))  


pdxME <- moduleEigengenes(pdx100, 
                           mod100, 
                           impute = TRUE, 
                           nPC = 1, 
                           align = "along average", 
                           subHubs = TRUE,
                           softPower = 6,
                           scale = TRUE,
                           verbose = 0, indent = 0)

##works


ME_pdx100   = pdxME$eigengenes
distPC_tcga = 1-abs(cor(ME_tcga,use="p"))
distPC_tcga = ifelse(is.na(distPC_tcga), 0, distPC_tcga)
pcTree_tcga = hclust(as.dist(distPC_tcga),method="a") 
MDS_tcga  = cmdscale(as.dist(distPC_tcga),2)

rownames(ME_pdx100) <- rownames(pdx100)

write.table(ME_pdx100, file = "ME_pdx_100.txt", sep = "\t")


### metabric
metabric100 <- as.data.frame(t(metabric100))  


metabricME <- moduleEigengenes(metabric100, 
                          mod100, 
                          impute = TRUE, 
                          nPC = 1, 
                          align = "along average", 
                          subHubs = TRUE,
                          softPower = 6,
                          scale = TRUE,
                          verbose = 0, indent = 0)

##works


ME_metabric100   = metabricME$eigengenes
distPC_tcga = 1-abs(cor(ME_tcga,use="p"))
distPC_tcga = ifelse(is.na(distPC_tcga), 0, distPC_tcga)
pcTree_tcga = hclust(as.dist(distPC_tcga),method="a") 
MDS_tcga  = cmdscale(as.dist(distPC_tcga),2)

rownames(ME_pdx100) <- rownames(pdx100)

write.table(ME_metabric100, file = "ME_metabric_100.txt", sep = "\t")

rownames(ME_1A) <- colnames(datExpr1)


save(adjacency1, adjacency2, dissTOM1, dissTOM2, file = "Adjacency_and_distance_table.RData")
write.table(ME_1A, file = "ME_TCGA_all.txt", sep = "\t")
