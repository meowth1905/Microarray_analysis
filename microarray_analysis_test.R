# BiocManager::install("GEOquery")
# BiocManager::install("hgu133plus2.db")



library(affy)
library(oligo)
library(limma)
library(hgu133plus2.db)
library(pd.clariom.s.mouse)
library(org.Hs.eg.db)

sessionInfo()

# To get information on a bioconductor package
browseVignettes("affy")

# The ExpressionSet class
# commonly used to store microarray experiment data.
# A typical ExpressionSet class object contains the following data:
# 1. assayData: Numeric matrix of raw or normalized intensities, Rows = probes (or probe sets), columns = samples.
# 2. phenoData: Metadata about samples, Rows = samples, columns = sample attributes.
# 3. featureData: - probe-to-gene mapping (gene symbols, Entrez IDs, etc.), usually populated by loading annotation packages (hgu133a2.db, ..)
# 4. annotation: - a character string naming the platform (e.g., "hgu133a2")

# GEOquery: for reading files directly from GEO based on the accession
browseVignettes('GEOquery')
library(GEOquery) 

# GSM820817 	DKAT MEGM rep1
gsm <- getGEO('GSM820817')

## What methods are available to call on this class of data?
methods(class=class(gsm))

# Meta(gsm)    information about experiment.
head(Meta(gsm))

# $data_row_count
# [1] "54675"

Columns(gsm)  # 2 columns, ID_REF (probe ID) and intensity column

# Table(gsm)  # get the intensity values
head(Table(gsm))

# To get all samples under GSE33146:
# (
# GSM820817 	DKAT MEGM rep1
# GSM820818 	DKAT MEGM rep2
# GSM820819 	DKAT MEGM rep3
# GSM820820 	DKAT SCGM day14 rep1
# GSM820821 	DKAT SCGM day14 rep2
# GSM820822 	DKAT SCGM day14 rep3
# )


# processed expression matrix,
# ready for downstream analysis 
# preprocessing choices were made by the original authors.

gse <- getGEO('GSE33146')

class(gse)       # a list
# Information associated with a microarray experiment.
gse[[1]]

class(gse[[1]]) # ExpressionSet object

# we only need this


gse <- gse[[1]]
exprs(gse)      # matrix of expression values
pData(gse)      # sample annotations
fData(gse)      # feature annotations

# Download .CEL files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33146

setwd("C:\\MicroArrayAnalysis")

# Read raw probe intensities

library(oligo)

list.files()

cel.files <- list.celfiles("GSE33146_RAW", listGzipped=TRUE, full.names=TRUE)
cel.files

raw.data <- read.celfiles(cel.files)

raw.data
# assayData: 1354896 features, 6 samples 

image(raw.data[,1])
image(raw.data[,6])



