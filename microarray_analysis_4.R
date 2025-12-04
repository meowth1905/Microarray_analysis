setwd("C:\\Neethu\\Data_Analysis_With_R")

list.files()
pheno.data <- read.csv("assig6_data\\arrays.txt", sep = "\t", header = TRUE,
                       row.names = 2)
dim(pheno.data)

array.data <- read.csv("assig6_data\\arraydata.txt", sep = "\t", header = TRUE)
dim(array.data)

# Make a Boxplot
boxplot(array.data[, rownames(pheno.data)], las = 2, cex.axis = 0.7)
# intensity data is not normal (as expected)

# Log2 transformation of array data
array.log <- log2(array.data[, rownames(pheno.data)])

# Make a Boxplot
boxplot(array.log, las = 2, cex.axis = 0.7)

# Normalization
library(preprocessCore)

class(array.log)

array.norm <- normalize.quantiles(as.matrix(array.log))
colnames(array.norm) <- colnames(array.log)
rownames(array.norm) <- rownames(array.log)

# Make a Boxplot
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


# Fitting a linear model
# 2 entries (HT and GF) for each Individual (here mouse)
# example: pheno.data[pheno.data$Individual=="513",]
# Linear model we assume observations are independent.
# Not the case here
# So we split data


# linear model for HT Tissue
#
pheno.data[,"Tissue"] == "HT"
# or
pheno.data$Tissue == "HT"

which(pheno.data$Tissue == "HT")
HT.samples <- rownames(pheno.data)[which(pheno.data$Tissue == "HT")]

# expression values of 1st probe for HT samples
array.norm[1, HT.samples]

# strains of HT samples
# there are 3 different strains: "B6N"  "BFMI" "F1" 
pheno.data[HT.samples, "Strain"]    

row.number <- 500
expr <- array.norm[row.number, HT.samples]
strain <- factor(pheno.data[HT.samples, "Strain"])
lm(expr ~ strain)
summary(lm(expr ~ strain))

# Or use anova to get p value
anova(lm(expr ~ strain))
anova(lm(expr ~ strain))["strain", "Pr(>F)"]

p.vals <- c()
for (i in 1: nrow(array.norm)){
  expr <- array.norm[i, HT.samples]
  strain <- factor(pheno.data[HT.samples, "Strain"])
  p.val <- anova(lm(expr ~ strain))["strain", "Pr(>F)"]
  p.vals <- c(p.vals, p.val)
}

# Adjusting p values for multiple testing
p.adjusted <- p.adjust(p.vals, method = "fdr")

# How many below 0.05
length(which(p.adjusted < 0.05))

# so in 1198 probes one of the strain is different from the other two

# 
# To find which probes are different between
# BFMI and B6n
# F1 and BFMI
# PostHoc tests for figuring which of the strains are really different

which(p.adjusted < 0.05)
array.norm[88, HT.samples]

strain <- pheno.data[HT.samples, "Strain"]
array.norm[88, HT.samples][which(strain == "F1")]

n.tests <- length(which(p.adjusted < 0.05))
bfmi.b6n <- 0; bfmi.f1 <- 0
b6n.probes <- c(); f1.probes <- 0

for (i in which(p.adjusted < 0.05)){
  expr <- array.norm[i, HT.samples]
  strain <- factor(pheno.data[HT.samples, "Strain"])
  
  bfmi <- expr[which(strain == "BFMI")]
  f1 <- expr[which(strain == "F1")]
  b6n <- expr[which(strain == "B6N")]
  
  t1 <- t.test(bfmi, b6n)
  t2<- t.test(bfmi, f1)
  
  if (p.adjust(t1$p.value, "fdr", n = n.tests) < 0.05){
    bfmi.b6n <- bfmi.b6n + 1
    b6n.probes <- c(b6n.probes, i)
  }
  
  if (p.adjust(t2$p.value, "fdr", n = n.tests) < 0.05){
    bfmi.f1 <- bfmi.f1 + 1
    f1.probes <- c(f1.probes, i)
  }

}

bfmi.b6n
# 96 probes different between bfmi and b6n 

bfmi.f1
# 40 probes different between bfmi and b6n

# Instead of copying above code and replacing HT with GF

linear.model <- function(normalized.array.data, samples){
  p.vals <- c()
  for (i in 1: nrow(normalized.array.data)){
    expr <- normalized.array.data[i, samples]
    strain <- factor(pheno.data[samples, "Strain"])
    p.val <- anova(lm(expr ~ strain))["strain", "Pr(>F)"]
    p.vals <- c(p.vals, p.val)
  }
  return(p.vals)
}

posthoc <- function(p.adjusted, samples, pheno.data, normalized.array.data){
  
  bfmi.b6n <- 0; bfmi.f1 <- 0; b6n.f1 <- 0
  b6n.probes <- c(); f1.probes <- 0; b6.f1.probes <- c()
  
  for (i in which(p.adjusted < 0.05)){
    expr <- normalized.array.data[i, samples]
    strain <- factor(pheno.data[samples, "Strain"])
    
    bfmi <- expr[which(strain == "BFMI")]
    f1 <- expr[which(strain == "F1")]
    b6n <- expr[which(strain == "B6N")]
    
    t1 <- t.test(bfmi, b6n)
    t2 <- t.test(bfmi, f1)
    t3 <- t.test(b6n, f1)
    
    n.tests <- length(which(p.adjusted < 0.05))
    
    if (p.adjust(t1$p.value, "fdr", n = n.tests) < 0.05){
      bfmi.b6n <- bfmi.b6n + 1
      b6n.probes <- c(b6n.probes, i)
    }
    
    if (p.adjust(t2$p.value, "fdr", n = n.tests) < 0.05){
      bfmi.f1 <- bfmi.f1 + 1
      f1.probes <- c(f1.probes, i)
    }
    
    if (p.adjust(t3$p.value, "fdr", n = n.tests) < 0.05){
      b6n.f1 <- b6n.f1 + 1
      b6.f1.probes <- c(b6.f1.probes, i)
    }

  }
  return(list("bfmi vs b6n"=bfmi.b6n, "bfmi vs f1"=bfmi.f1, "b6n vs f1"=b6n.f1,
              "probe index bfmi vs b6n"=b6n.probes, "probe index bfmi vs f1"=f1.probes,
              "probe index b6n vs f1"=b6.f1.probes))
}


p.adj.HT <- p.adjust(linear.model(array.norm, HT.samples), method = "fdr")
length(which(p.adj.HT < 0.05))

result.HT <- posthoc(p.adj.HT, HT.samples, pheno.data, array.norm)
result.HT


# GF samples
GF.samples <- rownames(pheno.data)[which(pheno.data$Tissue == "GF")]

p.adj.GF <- p.adjust(linear.model(array.norm, GF.samples), method = "fdr")
length(which(p.adj.GF < 0.05))

result.GF <- posthoc(p.adj.GF, GF.samples, pheno.data, array.norm)
result.GF
