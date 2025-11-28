#
# Microarray data analysis
# Neethu Raj
#

setwd("C:\\MicroArrayAnalysis")
old.par <- par()

library(oligo)

list.files()
cel.files <- list.celfiles("Data\\female_mouse.data", full.names = TRUE)
raw.data <- read.celfiles(cel.files)
raw.data
dim(raw.data)

# Extract the expression matrix
mat <- exprs(raw.data)

# Distribution of raw probe intensities per CEL file (per array)
boxplot(mat, main="Raw intensities",
        col=c(rep("coral", 3), rep("gray", 3)),
        las = 2,
        cex.axis = 0.7)

# Microarray raw intensities are usually very skewed, often ranging from 0 to tens of thousands.
boxplot(log2(mat), main="Raw intensities",
        col=c(rep("coral", 3), rep("gray", 3)),
        las = 2,
        cex.axis = 0.7)

# Alternatively we can also directly plot raw.data
class(raw.data)

# boxplot() method for ExpressionFeatureSet or ExpressionSet automatically
# applies a log transformation to make the intensity distribution easier to visualize.

boxplot(raw.data, main="Raw intensities",
        col=c(rep("coral", 3), rep("gray", 3)),
        las = 2,
        cex.axis = 0.7)

# Histograms

hist(mat[,1],
     breaks = 100,
     main = colnames(mat)[1],
     xlab="Intensity",
     col="skyblue",
     cex.main = 0.8)

hist(log2(mat[,1]),
     breaks = 100,
     main = colnames(mat)[1],
     xlab="Intensity",
     col="skyblue",
     cex.main = 0.8)

# All samples together
par(mfrow=c(2,3))
for (i in 1:ncol(mat)) {
  hist(mat[,i],
       breaks = 100,
       main = colnames(mat)[i],
       xlab="Intensity",
       col="skyblue",
       cex.main = 0.8)
}
par(old.par)

# rma() stands for Robust Multi-array Average.
# It does three steps: 
# 1. Background correction – reduces noise.
# 2. Quantile normalization – makes distributions across arrays comparable.
# 3. Probe summarization – takes all the probes belonging to a probe set and
# combines their intensities into one value per probe set.

norm.data <- rma(raw.data)
dim(norm.data)

# matrix of normalized expression values
mat.norm <- exprs(norm.data) 

# Distribution of normalized probe intensities per CEL file (per array)
boxplot(mat.norm, main="Normalized intensities",
        col=c(rep("coral", 3), rep("grey", 3)),
        las = 3,
        cex.axis = 0.7)

# Run PCA
pca <- prcomp(t(mat.norm), scale. = TRUE)

# Plot PCA
short.names <- sub(".CEL$", "", colnames(norm.data))
short.names <- gsub("tibia_female","", short.names)
short.names <- gsub("weekold", "w", short.names)
short.names <- gsub("__", "_", short.names)

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

text(pca$x[,1], pca$x[,2], labels = short.names, pos = 3, cex = 0.6)

# young females are clustered together

# Perform hierarchical clustering
hc <- hclust(dist(t(mat.norm)))

plot(hc, main="Hierarchical Clustering of Samples", 
     xlab="", sub="", cex=0.8)

# Create pheno data for samples
# Each row = one sample
# Each column = one experimental group
pheno.df <- data.frame(
  sample = colnames(mat.norm),
  age = factor(c("Young", "Young", "Young", "Old", "Old", "Old"))
  )

# Limma requires a numeric matrix (design) as input for lmFit().
design <- model.matrix(~0 + age, data = pheno.df)
colnames(design) <- c(levels(pheno.df$age))
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

# group_by replace rowname with numbers
# save the probe ids in a column
results.annotated$ProbeID <- rownames(results.annotated)

results.summarized <- results.annotated %>%
  group_by(GeneSymbol) %>%          # group rows by gene symbol
  slice_min(P.Value, n = 1)         # pick one with the smallest p-value per gene

class(results.summarized)

results.summarized <- as.data.frame(results.summarized)
rownames(results.summarized) <- results.summarized$ProbeID

#Filter significant genes
sig.genes <- results.summarized[
  results.summarized$adj.P.Val < 0.05 & abs(results.summarized$logFC) > 1,
]
dim(sig.genes) # 0 8

sig.genes <- results.summarized[
  results.summarized$P.Value < 0.05 & abs(results.summarized$logFC) > 1.5,
]
dim(sig.genes) # 18  8

top.genes <- results.summarized[order(results.summarized$P.Value), ][1:20, ]

# Create x and y
x <- results.summarized$logFC
y <- -log10(results.summarized$P.Value)

# volcano plot
plot(x, y,
     pch=20,
     col = "grey",
     main="Old vs young females",
     xlab="log2 Fold Change",
     ylab="-log10(P value)")

sigs <- results.summarized$P.Value < 0.05 & abs(results.summarized$logFC) > 1.5
points(x[sigs], y[sigs], col="red", pch=20)

legend("bottomleft",                     
       legend=c("Non-significant", "Significant (P.value < 0.05)"),
       col=c("grey", "red"),             
       pch=20,
       bty = "n",
       y.intersp=0.7,
       cex = 0.8)                           

# Add a column to mark significance
results.summarized$Gene <- ifelse(
  results.summarized$P.Value < 0.05 & abs(results.summarized$logFC) > 1.5,
  "Significant",
  "Non-significant"
)

library(ggrepel)
library(ggplot2)

ggplot(results.summarized, aes(x=logFC, y=-log10(P.Value), color=Gene)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("grey", "red")) +
  xlab("log2 Fold Change") +
  ylab("-log10(P.Value)") +
  ggtitle("Old vs Young females") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="black") +
  geom_text_repel(
    data = sig.genes,  # keep rows where SigGene=T
    aes(x = logFC, y = -log10(P.Value), label = GeneSymbol),
    color = "forestgreen",
    size = 3
  ) +
  theme_minimal()

head(top.genes)

# Get intensity values for the top 20 genes
mat.top <- mat.norm[rownames(mat.norm) %in% rownames(top.genes), ]

dim(mat.top)
head(mat.top)

heatmap(mat.top, cexRow = 0.7, cexCol=0.7)

colnames(mat.top)
colnames(mat.top) <- short.names

gene.labels <- top.genes[rownames(mat.top), "GeneSymbol"]
rownames(mat.top) <- gene.labels

heatmap(mat.top, cexRow = 0.7, cexCol=0.7)

