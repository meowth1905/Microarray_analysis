#
# Microarray data analysis
# Neethu Raj
#

setwd("C:\\MicroArrayAnalysis")
old.par <- par()

library(oligo)

list.files()
cel.files <- list.celfiles("Data\\E-MTAB-15951", full.names = TRUE)
raw.data <- read.celfiles(cel.files)
raw.data
dim(raw.data)

# Distribution of raw probe intensities per CEL file (per array)
boxplot(raw.data, main="Raw intensities",
        col=c(rep("coral", 3), rep("gray", 3)),
        las = 2,
        cex.axis = 0.7)

# Extract the expression matrix
mat <- exprs(raw.data)

par(mfrow=c(2,3))
for (i in 1:ncol(mat)) {
  hist(mat[,i],
       breaks = 100,
       main = colnames(mat)[i],
       xlab="Intensity",
       col="skyblue",
       cex.main = 0.8)
}


# RMA does three steps: background correction → quantile normalization → probe summarization.
norm.data <- rma(raw.data)

par(old.par)
# Distribution of normalized probe intensities per CEL file (per array)
boxplot(norm.data, main="Normalized intensities",
        col=c(rep("coral", 3), rep("grey", 3)),
        las = 3,
        cex.axis = 0.7)


# matrix of normalized expression values
mat.norm <- exprs(norm.data) 

# Run PCA
pca <- prcomp(t(mat.norm), scale. = TRUE)

# Plot PCA
short.names <- sub(".CEL$", "", colnames(norm.data))
short.names <- gsub("tibia_female", "_", short.names)
short.names <- gsub("weekold", "w", short.names)

plot(pca$x[,1], pca$x[,2],
     pch = 20,
     col = c(rep("orange", 3), rep("gray", 3)),
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of RMA-normalized data",
     bty = "n",
     cex.main = 0.8)

legend("topright",                     
       legend = c("Young female", "Old female"), 
       col = c("orange", "gray"), 
       pch = 20,
       bty = "n",
       cex = 0.7)

text(pca$x[,1], pca$x[,2], labels = short.names, pos = 3, cex = 0.7)

# young females are clustered together

# Perform hierarchical clustering
hc <- hclust(dist(t(mat.norm)))

plot(hc, main="Hierarchical Clustering of Samples", 
     xlab="", sub="", cex=0.8)



group <- factor(c(rep("Young", 3), rep("Old", 3)))
levels(group)

# Create a design matrix for old vs young
# Each row = one sample
# Each column = one experimental group

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

library(limma)
# Fit linear model with limma
fit <- lmFit(mat.norm, design)

# Create a contrast (Old vs Young)
contrast.matrix <- makeContrasts(Old - Young, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#Extract DE results
results <- topTable(fit2, adjust.method = "BH", number = Inf)
head(results)
dim(results)

# Map probes to genes
#
annotation(raw.data)    #  "pd.clariom.s.mouse"
BiocManager::available(pattern="clariomsmouse")

# BiocManager::install("clariomsmousetranscriptcluster.db")
library(clariomsmousetranscriptcluster.db)
library(AnnotationDbi)

columns(clariomsmousetranscriptcluster.db)

results$GeneSymbol <- mapIds(clariomsmousetranscriptcluster.db,
                             keys=rownames(results),
                             column="SYMBOL",
                             keytype="PROBEID",
                             multiVals="first")
head(results)

results.annotated <- results[!is.na(results$GeneSymbol), ]
dim(results.annotated)

# Many probes can map to the same gene symbol.
# treats all probes of the same gene as one group.
# from each group, pick the probe with the smallest p-value, assuming it is the “most significant” measurement for that gene.

table(duplicated(results.annotated$GeneSymbol))
# 20737 FALSE → 20,737 gene symbols are unique at their first appearance
# 258 TRUE → 258 entries are duplicates

library(dplyr)

results.annotated$ProbeID <- rownames(results.annotated)
# this is done because group_by replace rowname with numbers

results.byGene <- results.annotated %>%
  group_by(GeneSymbol) %>%           # group rows by gene symbol
  slice_min(P.Value, n = 1)          # pick the row (probe) with the smallest p-value per gene

results.byGene <- as.data.frame(results.byGene)
rownames(results.byGene) <- results.byGene$ProbeID

dim(results.byGene)

#Filter significant genes
sig.genes <- results.byGene[
  results.byGene$adj.P.Val < 0.05 & abs(results.byGene$logFC) > 1,
]
dim(sig.genes) # 0 8

sig.genes <- results.byGene[
  results.byGene$P.Value < 0.05 & abs(results.byGene$logFC) > 1.5,
]
dim(sig.genes) # 18  8

topgenes <- results.byGene[order(results.byGene$P.Value), ][1:20, ]


# Create x and y
x <- results.byGene$logFC
y <- -log10(results.byGene$P.Val)

# volcano plot
plot(x, y,
     pch=20,
     col = "grey",
     main="Old vs young females",
     xlab="log2 Fold Change",
     ylab="-log10(P-value)")

sigs <- results.byGene$P.Value < 0.05 & abs(results.byGene$logFC) > 1.5
# sum(results.byGene$P.Value < 0.05 & abs(results.byGene$logFC) > 1.5)

points(x[sigs], y[sigs], col="red", pch=20)

legend("bottomleft",                     
       legend=c("Non-significant", "Significant (P.value < 0.05)"),
       col=c("grey", "red"),             
       pch=20,
       bty = "n",
       y.intersp=0.7)                           


# Add a column to mark significance
results.byGene$Significant <- results.byGene$P.Value < 0.05 & abs(results.byGene$logFC) > 1.5

library(ggrepel)
library(ggplot2)

ggplot(results.byGene, aes(x=logFC, y=-log10(P.Value), color=Significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("grey", "red")) +
  xlab("log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Old vs Young females") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="black") +
  geom_text_repel(
    data = subset(results.byGene, Significant),  # only significant genes
    aes(x = logFC, y = -log10(P.Value), label = GeneSymbol),
    color = "forestgreen",
    size = 3
  ) +
  theme_minimal()
 
#==============================
# in progress

# exprs.data from your normalized data
exprs.top <- exprs.data[rownames(exprs.data) %in% top.30$ProbeID, ]

probe2gene <- setNames(top.30$GeneSymbol, top.30$ProbeID)
rownames(exprs.top) <- probe2gene[rownames(exprs.top)]
dim(exprs.top)

exprs.scaled <- t(scale(t(exprs.top)))  # scale by row (gene)

library(pheatmap)
pheatmap(exprs.scaled)
