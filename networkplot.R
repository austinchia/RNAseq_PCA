# BiocManager::install("tximport")
# BiocManager::install("rhdf5")

library("tximport")
library("rhdf5")
library(DESeq2)
library(piano)
library(EnhancedVolcano)

# Loading GSC from jarkko's folder
load("~/Desktop/ArabidopsisGO.RData")

# vignette("tximport")
# tximport function
dir <- "/Users/austinchia/Desktop/SRR_files"
samples <- read.table(file.path(dir, "sampletext.txt"),header=TRUE)
files <- file.path(dir, "kallisto", samples$run, "abundance.h5")
names(files) <- paste0(samples$run)

txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene.nosplice)



sampleTable <- data.frame(condition = factor(rep(c("WT", "Spfs2"), each = 3)))
rownames(sampleTable) <- colnames(txi.kallisto$counts)
X <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
X_before <- X
nrow(X)

#Remove genes expressed more than 1 replicate
X_after <- X[rowSums(counts(X))>1, ]
nrow(X)

# DESeq
X.de <- DESeq(X)
results(X.de)
res <- results(X.de, alpha = 0.05)
mcols(res)
summary(res)


X_after.de <- DESeq(X_after)
results(X_after.de)
res2 <- results(X_after.de, alpha = 0.05)
mcols(res2)
summary(res2)

#Visualize
plotMA(res, ylim = c(-8, 8))
hist(res$pvalue)

plotMA(res2, ylim = c(-8, 8))
hist(res$pvalue)


# Volcano Plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

res$padj <- res$padj
res$padj[is.na(res$padj)]=1
names(res$padj)=rownames(res)


resSig <- res[res$padj < 0.05,]
write.table(resSig,".All.txt")
write.table(subset(resSig, log2FoldChange > 0.58), file="./UpRegulated.txt")
write.table(subset(resSig, log2FoldChange < 0.58), file="./DownRegulated.txt")

# Transform data for visualization purposes
rld_after <- rlog(X, blind = FALSE)

# Principal component (PCA) plot of data
plotPCA(rld, intgroup = c("treatment"))

# hierarchical clustering
sampleDists_after <- dist(t(assay(rld_after)))
sampleLabels_after <-paste(rld_after$condition,colnames(rld_after),sep=":")
plot(hclust(sampleDists_after),labels=sampleLabels_after)

## Create input vectors for GSA
pvals <- res$pvalue
directions <- res$log2FoldChange

# The vector items need names, otherwise we don’t know which p-value is which.
names(pvals) <- rownames(res)
names(directions) <- rownames(res)

# Some p-values are NAs, set them to 1 and NAs in directions to a small value.
pvals[is.na(pvals)] <- 1
directions[is.na(directions)] <- 1e-10

# Running GSA
gsa.res <- runGSA(pvals,directions, gsc= GSC, gsSizeLim=c(5,300),
                  nPerm = 1000)
GSAsummaryTable(gsa.res,save=T,file="gsares.xls")

# Plot a network (will not work for the airway data, too many enriched GOs):
nw <- networkPlot2(gsa.res,class="distinct",significance=0.05, direction = "both")
View(nw)

GO:0009733	
GO:0006979	
GO:0071555	
GO:0003824	
GO:0045893	
GO:0005215	

saveRDS(nw, file = "nw.rds")
