# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

## DESeq2
# BiocManager::install("airway")
library("airway")
data("airway")
head(airway)
colData(airway)
counts<-assay(airway)

## Preparing a data set for DESeq2:
# Count matrix:
counts<-assay(airway)
## Experiment information - information on case vs. control samples:
treatment <- relevel(airway$dex, "untrt")
## Set up as DESeqDataSet object:

#BiocManager::install("DESeq2")
library(DESeq2)
X = DESeqDataSetFromMatrix(counts,DataFrame(treatment), ~ treatment)

##Filtering
# Check the number of rows in the data:
nrow(X)
## [1] 64102
## Remove all genes that are expressed only in a single replicate:
X <- X[rowSums(counts(X))> 1, ]
nrow(X)
## [1] 29391

## Transform data for visualization purposes
rld <- rlog(X, blind = FALSE)

# Uses experimental design in normalization, which may emphasize
# treatment differences (blind=T doesn’t use)

# Hierachical clustering
sampleDists <- dist(t(assay(rld)))
sampleLabels<-paste(rld$treatment,colnames(rld),sep=":")
plot(hclust(sampleDists),labels=sampleLabels)

## First, open a pdf device:
pdf("myplot.pdf")

# Then make plots:
plotPCA(rld, intgroup = c("treatment"))
plot(hclust(sampleDists),labels=sampleLabels)

# Then close the device:
dev.off()

## Differential expression
# Calculate differential expression:
X.de <- DESeq(X)

# Look at results:
results(X.de)
res<-results(X.de,alpha=0.05)
mcols(res)
summary(res)

# Visualize:
plotMA(res,ylim=c(-8,8))
hist(res$pvalue,100)

## biomaRt: collection of databases with biological information
# BiocManager::install("biomaRt")
library(biomaRt)
listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
listAttributes(ensembl)
listFilters(ensembl)

## GO - gene associations
mapGO <- getBM(attributes=c("ensembl_gene_id","name_1006"), mart=ensembl)
mapGO.red<-mapGO[mapGO[,2]!="",]

## Collect GOs into gene sets
# A ready-made function in piano.
# BiocManager::install("piano")
library(piano)
myGsc<-loadGSC(mapGO.red)

## GO enrichment
# Create input vectors for GSA
pvals <- res$pvalue
directions <- res$log2FoldChange

# The vector items need names, otherwise we don’t know which p-value is which.
names(pvals)=rownames(res)
names(directions)=rownames(res)

# Some p-values are NAs, set them to 1 and NAs in directions to a small value.
pvals[is.na(pvals)]=1
directions[is.na(directions)]=1e-10

# Now we can run GSA
gsa.res<-runGSA(pvals,directions, gsc=GSC.splice, gsSizeLim=c(5,300),
                  nPerm = 1000)

## GO enrichment Part 2

# Get information about GSA results
# Save results as an Excel sheet:
GSAsummaryTable(gsa.res,save=T,file="gsares.xlsx")

# Plot a network (will not work for the airway data, too many enriched GOs):
nw <- networkPlot2(gsa.res,class="non",significance=0.05, physics = FALSE, lay = "2")

