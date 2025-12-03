setwd("C:\\Neethu\\MicroArrayAnalysis")

list.files()
pheno.data <- read.csv("assig6_data\\arrays.txt", sep = "\t", header = TRUE)
dim(pheno.data)

# HT: Hypothalamus
# GF: Gonadal fat

array.data <- read.csv("assig6_data\\arraydata.txt", sep = "\t", header = TRUE)
dim(array.data)

array.num <- array.data[, 2:17]
boxplot(array.num, las = 2, cex.axis = 0.7)

array.log <- apply(array.num, 2, log2)
boxplot(array.log, las = 2, cex.axis = 0.7)

library(preprocessCore)
array.norm <- normalize.quantiles(array.log)
colnames(array.norm) <- colnames(array.log)
rownames(array.norm) <- rownames(array.log)

boxplot(array.norm, las = 2, cex.axis = 0.7)

old.par <- par()

par(mfrow=c(1,2))
boxplot(array.log, las = 2, cex.axis = 0.7, main="Before normalization")
boxplot(array.norm, las = 2, cex.axis = 0.7, main="After normalization")

par(old.par)

corr <- cor(array.norm, method = "pearson")
heatmap(corr)

library(pheatmap)
pheatmap(corr)

dim(corr)

par(mfrow=c(2,3))
for (i in 1:6){
  barplot(corr[i, -1], # exclude self-correlation
        main=colnames(corr)[i],
        col="skyblue", las=2)
}

par(mfrow=c(2,3))
for (i in 7:12){
  barplot(corr[i, -1], # exclude self-correlation
          main=colnames(corr)[i],
          col="skyblue", las=2)
}

par(mfrow=c(2,2))
for (i in 13:16){
  barplot(corr[i, -1], # exclude self-correlation
          main=colnames(corr)[i],
          col="skyblue", las=2)
}

corr.spear <- cor(array.norm, method = "spearman")
heatmap(corr.spear)

HT <- array.norm[, grep("^HT", colnames(array.norm))]
GF <- array.norm[, grep("^GF", colnames(array.norm))]

# or

HT.samples <- pheno.data[pheno.data[,"Tissue"] == "HT", "AtlasID"]
HT <- array.norm[,HT]

GF.samples <- pheno.data[pheno.data[,"Tissue"] == "GF", "AtlasID"]
GF <- array.norm[,GF]


apply(HT, 1, mean)
apply(HT, 1, sd)

apply(GF, 1, mean)
apply(GF, 1, sd)

t.test(HT, GF)

results <- NULL
for (i in rownames(array.norm)){
  ht <- as.numeric(HT[i,])
  gf <- as.numeric(GF[i,])
  values <- c(mean(ht), mean(gf), sd(ht), sd(gf), t.test(ht,gf)$p.value)
  results <- rbind(results, values)
}


results <- NULL
for (i in rownames(array.norm)){
  ht <- HT[i,]
  gf <- GF[i,]
  values <- c(mean(ht), mean(gf), sd(ht), sd(gf), t.test(ht,gf)$p.value)
  results <- rbind(results, values)
}
results <- as.data.frame(results)

colnames(results) <- c("meanHT","meanGF","sdHT","sdGF","pValues")
rownames(results) <- rownames(array.norm)

# p-value < 0.05
# probe’s intensities are significantly different between HT & GF samples
results$sig <- results$pValues < 0.05

?p.adjust

results$pAdjusted <- p.adjust(results$pValues, method = "BH")
results$sigAdjusted <- results$pAdjusted < 0.05

table(results$sig)
table(results$sigAdjusted)

