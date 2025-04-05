# Colleen Ahern
# MDS plot from G3 featureCounts count matrix

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))

library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
BiocManager::install("org.Mm.eg.db")

# Import the count data and sample info
G3seqdata <- read.delim("/Users/colleenahern/Documents/EBTNG3/G3_substrate_counts.txt", stringsAsFactors = FALSE)
G3sampleinfo <- read.delim("/Users/colleenahern/Documents/EBTNG3/g3sampleinfo.txt", stringsAsFactors = TRUE)
colnames(G3sampleinfo) <- gsub(".*\\.", "", colnames(G3sampleinfo))
head(G3seqdata)
dim(G3seqdata)

G3sampleinfo

# Remove first column from G3seqdata
G3countdata <- G3seqdata[,-(1)]
# Look at the output
head(G3countdata)

# Store GeneID as rownames
rownames(G3countdata) <- G3seqdata[,1]
head(G3countdata)

# Make sure that the column names are the same as libraryName in the G3sampleinfo file
colnames(G3countdata) # using substr, you can change the characters of the colnames
table(colnames(G3countdata)==G3sampleinfo$libraryName)
G3countdata <- G3countdata[,as.character(G3sampleinfo$libraryName)]
table(colnames(G3countdata)==G3sampleinfo$libraryName) # Now they match

G3sampleinfo$groupName <- gsub(".*\\_", "", G3sampleinfo$groupName)
G3sampleinfo

## TPM cutoff

G3TPM <- read.delim("/Users/colleenahern/Documents/EBTNG3/tpm_counts.txt", header=T, check.names = FALSE)
G3TPM$keep <- 0

for (i in 1:nrow(G3TPM)) {
  G3TPM$keep[i] <- ifelse(mean(c(G3TPM$HSOHC[i],G3TPM$HSOHG[i],G3TPM$HSOHH[i])) < 2 && mean(c(G3TPM$HSOHN[i],G3TPM$HWOOB[i],G3TPM$HSOHP[i],G3TPM$HSOHS[i])) < 2 && mean(c(G3TPM$HSOHT[i],G3TPM$HSOHU[i],G3TPM$HSOHW[i],G3TPM$HSOHX[i])) < 2 && mean(c(G3TPM$HSOHY[i],G3TPM$HSOHZ[i],G3TPM$HSONA[i],G3TPM$HSONB[i])) < 2, "remove", "keep") 
}

G3TPM <- G3TPM[G3TPM$keep == "keep",] # reduces genes from 23664 to 11749
G3TPM <- G3TPM[, -17]

G3countdata_tpmfilt <- G3countdata[rownames(G3countdata) %in% G3TPM$Geneid,]

## Convert counts to DGEList object
G3y_tpmfilt <- DGEList(G3countdata_tpmfilt)
G3y_tpmfilt

group <- paste(G3sampleinfo$groupName)
group

# Convert to factor
group <- factor(group)
group

# Add the group information into the DGEList
G3y_tpmfilt$samples$group <- group
G3y_tpmfilt$samples

G3y_tpmfilt$samples$lib.size  # the TPM filter has slightly reduced the library size

# Apply normalisation to DGEList object
G3y_tpmfilt_norm <- calcNormFactors(G3y_tpmfilt)

# This will update the normalization factors in the DGEList object (their default values are 1). 
G3y_tpmfilt_norm$samples
levels(as.factor(G3sampleinfo$groupName))
col.cell <- c("blue", "green", "purple","orange")[as.factor(G3sampleinfo$groupName)]
data.frame(G3sampleinfo$groupName,col.cell)
plotMDS(G3y_tpmfilt_norm,col=col.cell, pch=16)
legend("topleft",fill=c("blue", "green", "purple","orange"),legend=c("Cellobiose", "Filter paper", "Glucose", "Reed canary grass"), title = expression(bold("Substrate")))
title("Substrate type")


