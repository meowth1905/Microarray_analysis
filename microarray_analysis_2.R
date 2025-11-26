#
# Microarray data analysis
# Neethu Raj
#

# -------------------------------------------------------------
# Dataset description:
# -------------------------------------------------------------

# Format:
# - Rows: Each row corresponds to a single probeset (ID_REF), e.g., "1007_s_at".
#   Note: The first 66 rows contain metadata / experiment descriptions 
#     (e.g., who performed the experiment, contact info, array details) and not actual expression values.
# - Columns: 
#       - The first column is "ID_REF" (probeset identifier).
#       - Subsequent columns are samples (e.g., GSM188013, GSM188014, ...).
#       - Each value represents raw probe intensity for that probeset in that sample.
#
# Example:
# ID_REF      GSM188013  GSM188014  GSM188016 ...
# 1007_s_at   15630.2    17048.8    13667.5   ...
# 1053_at     3614.4     3563.22    2604.65   ...
#
# -------------------------------------------------------------

setwd("C:\\MicroArrayAnalysis")
list.files()

raw.data <- read.csv(file = "Assignments04_Data/GSE7765-GPL96-series-matrix.txt",
                       sep = "\t",
                       header = TRUE,
                       row.names = 1,
                       skip = 66) # skip first 66 lines

dim(raw.data)

boxplot(raw.data, las = 2, cex.axis = 0.7)

# boxplot works on columns
# las = 2 perpendicular (vertical for x-axis)
# cex.axis = 0.7 scales the size of the axis labels to 70% of the default.

hist(raw.data[,1], breaks = 100)

# probe intensities vary exponentially with target abundance
# a probe gives more signal when there is more target bound to it.
# relationship between “how much target is present” and “how bright the signal is” is often non-linear
# small increases in target can cause big increases in signal.
# PCR can make this even more extreme (PCR amplifies target exponentially).

# this means raw intensities can span several orders of magnitude (from very dim to extremely bright).
# log2 transform compresses the dynamic range
# microarrays are always log transformed.

log.data <- apply(raw.data, 2, log2)
dim(log.data)

head(raw.data,5)
head(log.data,5)

hist(log.data[,1], breaks = 100)

boxplot(log.data, las = 2, cex.axis = 0.7)

# BiocManager::install("preprocessCore")
library(preprocessCore)

log.norm.data <- normalize.quantiles(log.data)

boxplot(log.norm.data, las = 2, cex.axis = 0.7)

rownames(log.norm.data) <- rownames(log.data)
colnames(log.norm.data) <- colnames(log.data)

boxplot(log.norm.data, las = 2, cex.axis = 0.7)

# how Quantile normalization works?
# sort each sample from smallest to largest.
# compute the average at each rank across samples. For example, average all the smallest values, all the second smallest, etc.
# replace the original values with these averages, but put them back in the original order of each sample.

# after quantile normalization
# each sample has the same distribution of values.
# means, medians, percentiles — all become the same.
# but the relative ranking within a sample (which spot was originally higher or lower) is preserved.
#
# actual intensity values get “standardized” to the same shape, but their positions within each sample stay the same.

# each probe gives a signal depending on how much target binds.
# more target → higher intensity.
# if one sample accidentally got more RNA/DNA during pipetting, it would have higher overall intensities.
# quantile normalization removes that technical bias, so samples become comparable.
#
# if a sample really does have more biological target
# e.g., one tissue truly has higher overall expression, 
# normalization will also remove or reduce that real difference.
#
# normalization assumes that most differences between samples are technical, not biological.
# assumes only a subset of genes change across samples.
# can distort data if most genes truly differ in abundance across samples.

# How different samples are to each other
# measure distance between each sample

clusters <- hclust(dist(t(log.norm.data)))
plot(clusters)

# dist() works on rows 
# hclust unsupervised hierarchical clustering
# takes a distance matrix and builds a tree showing how samples cluster together.

# Remove rows (probes) with zero variance
log.norm.data.filtered <- log.norm.data[apply(log.norm.data, 1, var) != 0, ]

# Transpose so samples are rows (needed for PCA)
pca <- prcomp(t(log.norm.data.filtered), scale. = TRUE)

plot(pca$x[,1], pca$x[,2],
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of samples",
     type = 'n')

text(pca$x[,1], pca$x[,2],
     labels = colnames(log.norm.data.filtered),
     cex = 0.6) 


# PCA plot
plot(pca$x[,1], pca$x[,2],
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of samples",
     pch = 19)


# Statistical tests compare predefined groups.
# All we have here is raw intensity values and
# no information about which sample belong to which group
# without knowing sample groups, differential expression analysis is not possible.

# Instead look at genes (probe set IDs) that are highly variable across samples

# For each row in log.norm.data calculate var(row) / mean(row)
scaled.variation <- apply(log.norm.data, 1, function(x){
  var(x)/mean(x)})

length(scaled.variation)

plot(scaled.variation)
# most have low scaled.variation
# we can see probes are targeting something that is equally expressed across samples

# highly variable probes
hvp <- names(which(scaled.variation > 1))
hvp

heatmap(log.norm.data[hvp,], cexRow = 0.7, cexCol=0.7)
# can see the same samples that clustered together in hclust()
# also has similar expression patterns in the heatmap

# How to get gene names from probes

#library(GEOquery) 
gse <- getGEO('GSE7765')[[1]] # from 2nd line of data file !Series_geo_accession	"GSE7765" 

platform <- annotation(gse)
platform

# search in 
# https://bioconductor.org/packages/release/data/annotation/


library(AnnotationDbi)

BiocManager::install("hgu133a.db")
library(hgu133a.db)

# map your highly variable probes (hvp) to gene symbols
hvg <- mapIds(hgu133a.db,
                      keys = hvp,
                      column = "SYMBOL",
                      keytype = "PROBEID",
                      multiVals = "first")

as.character(hvg)