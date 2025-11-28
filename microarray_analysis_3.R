#
# Microarray data analysis
# Neethu Raj
#

setwd("C:\\MicroArrayAnalysis")
old.par <- par()

library(oligo)

male.files   <- list.celfiles("Data/male_mouse.data", full.names = TRUE)
female.files <- list.celfiles("Data/female_mouse.data", full.names = TRUE)
all.files <- c(male.files, female.files)

raw.data <- read.celfiles(all.files)
dim(raw.data)

colnames(raw.data)

short.names <- sub(".CEL$", "", colnames(raw.data))
short.names <- gsub("weekold", "w", short.names)
colnames(raw.data) <- short.names

# Extract the expression matrix
mat <- exprs(raw.data)

# Distribution of raw probe intensities per CEL file (per array)
boxplot(mat, main="Raw intensities",
        col=c(rep("firebrick", 6), rep("skyblue", 6)),
        las = 2,
        cex.axis = 0.7)

# Microarray raw intensities are usually very skewed, often ranging from 0 to tens of thousands.
boxplot(log2(mat), main="Raw intensities",
        col=c(rep("coral1", 3), rep("coral3", 3), rep("skyblue", 3), rep("steelblue", 3)),
        las = 2,
        cex.axis = 0.7)

# Alternatively we can also directly plot raw.data
class(raw.data)

# boxplot() method for ExpressionFeatureSet or ExpressionSet automatically
# applies a log transformation to make the intensity distribution easier to visualize.

boxplot(raw.data, main="Raw intensities",
        col=c(rep("coral1", 3), rep("coral3", 3), rep("skyblue", 3), rep("steelblue", 3)),
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
        col=c(rep("coral1", 3), rep("coral3", 3), rep("skyblue", 3), rep("steelblue", 3)),
        las = 3,
        cex.axis = 0.7)
# OR

boxplot(norm.data, main="Normalized intensities",
        col=c(rep("coral1", 3), rep("coral3", 3), rep("skyblue", 3), rep("steelblue", 3)),
        las = 3,
        cex.axis = 0.7)

# Run PCA
pca <- prcomp(t(mat.norm), scale. = TRUE)

# Plot PCA
plot(pca$x[,1], pca$x[,2],
     pch = 20,
     col=c(rep("coral1", 3), rep("coral3", 3), rep("skyblue", 3), rep("steelblue", 3)),
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA of RMA-normalized data",
     bty = "n",
     cex.main = 0.8)

short.names <- sub(".CEL$", "", colnames(raw.data))
short.names <- gsub("^M[0-9]{2}_tibia_", "", short.names)
short.names <- gsub("weekold", "w", short.names)

text(pca$x[,1], pca$x[,2], labels = short.names, pos = 3, cex = 0.6)


# Perform hierarchical clustering
hc <- hclust(dist(t(mat.norm)))

plot(hc, main="Hierarchical Clustering of Samples", 
     xlab="", sub="", cex=0.8)

# Create pheno data for samples
# Each row = one sample
# Each column = one experimental group
pheno.df <- data.frame(
  sample = colnames(mat.norm),
  age = factor(c(rep("Young",3), rep("Old",3), rep("Young",3),  rep("Old",3))),
  gender = factor(c(rep("Male",6), rep("Female",6)))
)
# Check sample order matches between pheno.df and expression matrix:
all(pheno.df$sample == colnames(mat.norm))

# Limma requires a numeric matrix (design) as input for lmFit().
design <- model.matrix(~ age + gender, data = pheno.df)
design

levels(pheno.df$age)    # Reference = Old
levels(pheno.df$gender) # Reference = Female

# If you want Young and Female to be reference levels
# pheno.df$age <- relevel(pheno.df$age, ref="Young")
# pheno.df$gender <- relevel(pheno.df$gender, ref="Female")

library(limma)

# Fit linear model with limma
fit <- lmFit(mat.norm, design)
fit <- eBayes(fit)

# Extract DE results
# limma performs an F-test asking to answer 
# Is this gene affected by age OR gender (or both)?
res <- topTable(fit, coef = c("ageYoung", "genderMale"),
                    adjust.method = "BH", number = Inf)

head(res)

# Map probes to genes
#
annotation(raw.data)    #  "pd.clariom.s.mouse"
BiocManager::available(pattern="clariomsmouse")

# BiocManager::install("clariomsmousetranscriptcluster.db")
library(clariomsmousetranscriptcluster.db)
library(AnnotationDbi)

# Check which keys are available for this db
columns(clariomsmousetranscriptcluster.db)

res$GeneSymbol <- mapIds(clariomsmousetranscriptcluster.db,
                             keys=rownames(res),
                             column="SYMBOL",
                             keytype="PROBEID",
                             multiVals="first")
head(res)

res.annotated <- res[!is.na(res$GeneSymbol), ]
dim(res.annotated)

# Many probes can map to the same gene symbol.
# treats all probes of the same gene as one group.
# from each group, pick the probe with the smallest p-value, assuming it is the “most significant” measurement for that gene.

table(duplicated(res.annotated$GeneSymbol))
# 20737 FALSE → 20,737 gene symbols are unique at their first appearance
# 258 TRUE → 258 entries are duplicates

library(dplyr)

# group_by replace rowname with numbers
# save the probe ids in a column
res.annotated$ProbeID <- rownames(res.annotated)

res.summarized <- res.annotated %>%
  group_by(GeneSymbol) %>%          # group rows by gene symbol
  slice_min(P.Value, n = 1)         # pick one with the smallest p-value per gene

class(res.summarized)

res.summarized <- as.data.frame(res.summarized)
rownames(res.summarized) <- res.summarized$ProbeID

#Filter significant genes
# genes significantly affected by age, gender, or both.
sig.genes <- res.summarized[res.summarized$adj.P.Val < 0.05, ]
dim(sig.genes) 

# Get intensity values for the sig genes
mat.sig <- mat.norm[rownames(mat.norm) %in% rownames(sig.genes), ]
dim(mat.sig)

heatmap(mat.sig, cexRow = 0.7, cexCol=0.7)

colnames(mat.sig)
colnames(mat.sig) <- short.names

gene.labels <- sig.genes[rownames(mat.sig), "GeneSymbol"]
rownames(mat.sig) <- gene.labels

heatmap(mat.sig, cexRow = 0.7, cexCol=0.7)

#
#

# Age effect:
# ageYoung: how much higher or lower Young is compared to Old.
# genes differentially expressed due to age alone, while controlling for gender in the model.
res.age <- topTable(fit, coef="ageYoung", adjust.method="BH", number=Inf)

res.age$GeneSymbol <- mapIds(clariomsmousetranscriptcluster.db,
                         keys=rownames(res.age),
                         column="SYMBOL",
                         keytype="PROBEID",
                         multiVals="first")
head(res.age)
dim(res.age)

res.age <- res.age[!is.na(res.age$GeneSymbol), ]
dim(res.age)

table(duplicated(res.age$GeneSymbol))

res.age$ProbeID <- rownames(res.age)

res.age <- res.age %>%
  group_by(GeneSymbol) %>%         
  slice_min(P.Value, n = 1)         

class(res.age)

res.age <- as.data.frame(res.age)
rownames(res.age) <- res.age$ProbeID

#Filter significant genes
# genes significantly affected by age.
age.genes <- res.age[res.age$P.Value < 0.05 & abs(res.age$logFC) > 1, ]
dim(age.genes)  

# Get intensity values 
mat.age <- mat.norm[rownames(mat.norm) %in% rownames(age.genes), ]
dim(mat.age)

heatmap(mat.age, cexRow = 0.7, cexCol=0.7)

colnames(mat.age)
colnames(mat.age) <- short.names

gene.labels <- age.genes[rownames(mat.age), "GeneSymbol"]
rownames(mat.age) <- gene.labels

heatmap(mat.age, cexRow = 0.7, cexCol=0.7)
#
#
# Gender effect:
# genderMale: how much higher or lower Male is compared to Female

res.gender <- topTable(fit, coef="genderMale", adjust.method="BH", number=Inf)

res.gender$GeneSymbol <- mapIds(clariomsmousetranscriptcluster.db,
                             keys=rownames(res.gender),
                             column="SYMBOL",
                             keytype="PROBEID",
                             multiVals="first")
head(res.gender)
dim(res.gender)

res.gender <- res.gender[!is.na(res.gender$GeneSymbol), ]
dim(res.gender)

table(duplicated(res.gender$GeneSymbol))

res.gender$ProbeID <- rownames(res.gender)

res.gender <- res.gender %>%
  group_by(GeneSymbol) %>%         
  slice_min(P.Value, n = 1)         

class(res.gender)

res.gender <- as.data.frame(res.gender)
rownames(res.gender) <- res.gender$ProbeID

#Filter significant genes
# genes significantly affected by gender
gender.genes <- res.gender[res.gender$P.Value < 0.05 & abs(res.gender$logFC) > 1, ]
dim(gender.genes)  

# Get intensity values 
mat.gender <- mat.norm[rownames(mat.norm) %in% rownames(gender.genes), ]
dim(mat.gender)

heatmap(mat.gender, cexRow = 0.7, cexCol=0.7)

colnames(mat.gender)
colnames(mat.gender) <- short.names

gene.labels <- gender.genes[rownames(mat.gender), "GeneSymbol"]
rownames(mat.gender) <- gene.labels

heatmap(mat.gender, cexRow = 0.7, cexCol=0.7)
