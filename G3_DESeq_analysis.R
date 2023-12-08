# 07/18/2023
# Complete DESeq analysis on G3 count matrix from featureCounts

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")
library("tximport")
library("readr")
BiocManager::install("tximportData")
library("tximportData")
library(tidyverse)
# BiocManager::install("pasilla")
library(pasilla)

# Count matrix input
g3coldata <- read.delim("/Users/colleenahern/Documents/EBTNG3/g3sampleinfo.txt")
colnames(g3coldata) <- gsub(".*\\.", "", colnames(g3coldata))
g3coldata$groupName <- gsub(".*\\_", "", g3coldata$groupName)
g3cts <- read.delim("/Users/colleenahern/Documents/EBTNG3/G3_substrate_counts.txt")
rownames(g3cts) <- g3cts[,1]
g3cts <- g3cts[,-(1)]
g3coldata
rownames(g3coldata) <- g3coldata[,1]
g3coldata <- g3coldata[,-(1)]
g3cts
g3coldata

# Examine the count matrix and column data to see if they are consistent in terms of sample order
head(g3cts, 2)
g3coldata

# Rearrange
all(rownames(g3coldata) %in% colnames(g3cts))
all(rownames(g3coldata) == colnames(g3cts))

g3cts <- g3cts[, rownames(g3coldata)]
all(rownames(g3coldata) == colnames(g3cts))

# With the count matrix and the sample info, construct a DESeqDataSet:
library("DESeq2")
g3dds <- DESeqDataSetFromMatrix(countData = g3cts,
                              colData = g3coldata,
                              design = ~ groupName)
g3dds

# Pre-filter

keep <- rowSums(counts(g3dds)) > 0
g3dds <- g3dds[keep,]   # reduced number of genes from 23663 to 19384

# Set the factor levels 
g3dds$groupName <- factor(g3dds$groupName, levels = c("RCG","FP","glucose", "cellobiose"))

# Differential expression analysis: compare each substrate to glucose

g3dds <- DESeq(g3dds)
g3res <- results(g3dds)
g3res
g3res <- results(g3dds, contrast=c("groupName","RCG","FP"))
g3res

# Now compare all of the different the G3 substrates to each other: 6 in total
# RCG vs FP
g3res_RCG_FP <- results(g3dds, contrast=c("groupName","RCG","FP"))
g3res_RCG_FP
g3resmat_RCG_FP <- as.matrix(g3res_RCG_FP)
g3resmat_RCG_FP
colnames(g3resmat_RCG_FP) <- paste("RCG_v_FP", colnames(g3resmat_RCG_FP), sep = "_")
g3resmat_RCG_FP <- as.data.frame(g3resmat_RCG_FP)
sum(g3resmat_RCG_FP$RCG_v_FP_padj < 0.05, na.rm=TRUE)

# RCG vs glucose
g3res_RCG_glu <- results(g3dds, contrast=c("groupName","RCG","glucose"))
g3res_RCG_glu
g3resmat_RCG_glu <- as.matrix(g3res_RCG_glu)
g3resmat_RCG_glu
colnames(g3resmat_RCG_glu) <- paste("RCG_v_glucose", colnames(g3resmat_RCG_glu), sep = "_")

# RCG vs cellobiose
g3res_RCG_cel <- results(g3dds, contrast=c("groupName","RCG","cellobiose"))
g3res_RCG_cel
g3resmat_RCG_cel <- as.matrix(g3res_RCG_cel)
g3resmat_RCG_cel
colnames(g3resmat_RCG_cel) <- paste("RCG_v_cellobiose", colnames(g3resmat_RCG_cel), sep = "_")

# FP vs glucose
g3res_FP_glu <- results(g3dds, contrast=c("groupName","FP","glucose"))
g3res_FP_glu
g3resmat_FP_glu <- as.matrix(g3res_FP_glu)
g3resmat_FP_glu
colnames(g3resmat_FP_glu) <- paste("FP_v_glucose", colnames(g3resmat_FP_glu), sep = "_")

# FP vs cellobiose
g3res_FP_cel <- results(g3dds, contrast=c("groupName","FP","cellobiose"))
g3res_FP_cel
g3resmat_FP_cel <- as.matrix(g3res_FP_cel)
g3resmat_FP_cel
colnames(g3resmat_FP_cel) <- paste("FP_v_cellobiose", colnames(g3resmat_FP_cel), sep = "_")

# glucose vs cellobiose
g3res_glu_cel <- results(g3dds, contrast=c("groupName","glucose","cellobiose"))
g3res_glu_cel
g3resmat_glu_cel <- as.matrix(g3res_glu_cel)
g3resmat_glu_cel
colnames(g3resmat_glu_cel) <- paste("glucose_v_cellobiose", colnames(g3resmat_glu_cel), sep = "_")

# Confirm that the order of genes in all 6 matrices is the same
all(rownames(g3resmat_RCG_FP) == rownames(g3resmat_RCG_glu))
all(rownames(g3resmat_RCG_glu) == rownames(g3resmat_RCG_cel))
all(rownames(g3resmat_RCG_cel) == rownames(g3resmat_FP_glu))
all(rownames(g3resmat_FP_glu) == rownames(g3resmat_FP_cel))
all(rownames(g3resmat_FP_cel) == rownames(g3resmat_glu_cel))

CBA_DESeq <- cbind(g3resmat_RCG_FP, g3resmat_RCG_glu)
CBA_DESeq <- CBA_DESeq[, c(2,6,8,12)]
CBA_DESeq <- as.data.frame(CBA_DESeq)
CBA_DESeq$RCG_v_FP_Sig <- ifelse(CBA_DESeq$RCG_v_FP_padj <= 0.05,"TRUE","FALSE")
CBA_DESeq <- CBA_DESeq[, c(1,2,5,3,4)]
CBA_DESeq$RCG_v_glucose_Sig <- ifelse(CBA_DESeq$RCG_v_glucose_padj <= 0.05,"TRUE","FALSE")

CBA_DESeq <- cbind(CBA_DESeq, g3resmat_RCG_cel)
CBA_DESeq <- CBA_DESeq[, c(1:6,8,12)]
CBA_DESeq$RCG_v_cellobiose_Sig <- ifelse(CBA_DESeq$RCG_v_cellobiose_padj <= 0.05,"TRUE","FALSE")

CBA_DESeq <- cbind(CBA_DESeq, g3resmat_FP_glu)
CBA_DESeq <- CBA_DESeq[, c(1:9,11,15)]
CBA_DESeq$FP_v_glucose_Sig <- ifelse(CBA_DESeq$FP_v_glucose_padj <= 0.05,"TRUE","FALSE")

CBA_DESeq <- cbind(CBA_DESeq, g3resmat_FP_cel)
CBA_DESeq <- CBA_DESeq[, c(1:12,14,18)]
CBA_DESeq$FP_v_cellobiose_Sig <- ifelse(CBA_DESeq$FP_v_cellobiose_padj <= 0.05,"TRUE","FALSE")

CBA_DESeq <- cbind(CBA_DESeq, g3resmat_glu_cel)
CBA_DESeq <- CBA_DESeq[, c(1:15,17,21)]
CBA_DESeq$glucose_v_cellobiose_Sig <- ifelse(CBA_DESeq$glucose_v_cellobiose_padj <= 0.05,"TRUE","FALSE")

CBA_DESeq$GeneID <- rownames(CBA_DESeq)
CBA_DESeq <- CBA_DESeq[,c(19,1:18)]
  
write_csv(CBA_DESeq, "/Users/colleenahern/Documents/EBTNG3/CBA_DESeq_final.csv")


