#libraries
library(ape)
library(ggtree)
library(ggplot2)


# Reading Data Into Table
D <- read.table("mypca.eigenvec", row.names = 1)

# PCA PLot
pca_plot <- plot(D[,2],D[,3],type="n")
text(D[,2],D[,3],labels = D[,1],cex = 0.5,col = "red")

##=== Estimate phylogeny using IBS distances and UPGMA clustering (1)

#read in IBS distances
M <- read.delim("mydist.mdist", header = F)

#read species names
species <- read.delim("mydist.mdist.id",header = F)[,1]
dimnames(M) <- list(species,species)

#estimate phylogeny, convert to a tree object
tre <- as.phylo(hclust(as.dist(M), method = "average"))

#root the tree at outgroup
tre <- root(tre,outgroup = "Amborella")
tre <- ladderize(tre)

##=== Estimate phylogeny using IBS distances and UPGMA clustering (2)

# Calibrate the tree using fossil evidence
# Amborella angiosperm split approx 180 Mya.
# Assume evolutionary rates to be correlated
calibs <- makeChronosCalib(tre, node = "root", age.max = 180)
time.tre <- chronos(tre, lambda = 1, model = "correlated",
                    calibration = calibs, control = chronos.control())
write.tree(time.tre,file = "SpeciesTree_rooted_calib.tre")

##=== Estimate phylogeny using IBS distances and UPGMA clustering (3)

ftre <- read.tree("SpeciesTree_rooted_calib.tre")
ge <- ggtree(ftre,right = T) + 
  geom_tiplab(size = 3, fontface="italic") + 
  theme_tree2()
ge1 <- revts(ge)
ge1 <- ge1+ scale_x_continuous(breaks= c(-200, -150, -100, -50, 0), 
                               labels = c(-200,-150,-100,-50,0),limits=c(-200,50))
pdf("SNP_Phylogeny.pdf")
print(ge1)
dev.off()
