# script to identify deifferential gene expression using DESeq2 package 

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

#Step 1: preparing count data 

# read count matrix
counts_data <- read.csv('counts_data.csv')
head(counts_data)
rownames(counts_data)

# read sample info
colData <- read.csv('sample_info.csv')
head(colData)

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2. construct DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ dexamethasone
                       )

# pre-filtering removing rows with low gene counts
# keeping rows that have at least 10 reads total 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the factor level 
dds$dexamethasone <- relevel(dds$dexamethasone, ref="untreated")

# NOTE: collapse technical replicates

# Step3. run DESeq
dds <- DESeq(dds)
res <- results(dds)

# explore results 
summary(res)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# e.g. treated_4hrs, treated_8hrs, untreated
results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))


# MA plot 
plotMA(res)



