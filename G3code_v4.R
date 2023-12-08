# 06/12/2023 G3 transcriptomics analysis: MDS

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

# Convert counts to DGEList object
# This is an object used by edgeR to store count data. 
# It has a number of slots for storing various parameters about the data.

G3y <- DGEList(G3countdata)
G3y

# See what slots are stored in G3y
names(G3y)

# Library size information is stored in the samples slot
G3y$samples

group <- paste(G3sampleinfo$groupName)
group

# Convert to factor
group <- factor(group)
group

# Add the group information into the DGEList
G3y$samples$group <- group
G3y$samples

## Make an MDS plot of the unnormalized counts/no TPM filter (for comparison)

G3y$samples$lib.size

# Plot the library sizes as a barplot to see whether there are any major discrepancies between the samples 
# The names argument tells the barplot to use the sample names on the x-axis
# The last argument rotates the axis names
barplot(G3y$samples$lib.size,names=colnames(G3y),las=2)

# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labeling if we want
barplot(G3y$samples$lib.size/1e06, names=colnames(G3y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

## Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log 
## the counts. Use box plots to check the distribution of the read counts on the log2 scale. Use the cpm function to
## get log2 counts per million, which are corrected for the different library sizes. The cpm function also adds a small 
## offset to avoid taking log of zero.

# Get log2 counts per million
G3logcounts <- cpm(G3y,log=TRUE)

# Check distributions of samples using boxplots
boxplot(G3logcounts, xlab="", ylab="Log2 counts per million",las=2)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(G3logcounts),col="blue")
title("Boxplots of logCPMs \n (no TPM filter/TMM normalization)")

## From the boxplots we see that overall the density distributions of raw log-intensities are not identical but still
## not very different. If a sample is really far above or below the blue horizontal line we may need to investigate 
## that sample further. Another kind of QC plot that is helpful in checking for dodgy samples is a relative log 
## expression (RLE) plot, which can be generated with plotRLE from the EDASeq package.

plotMDS(G3y)

# Let's set up colour schemes for groupName
# How many cell types and in what order are they stored?
levels(as.factor(G3sampleinfo$groupName))

col.cell <- c("blue", "green", "purple","orange")[as.factor(G3sampleinfo$groupName)]
data.frame(G3sampleinfo$groupName,col.cell)

# Redo the MDS with cell type colouring
plotMDS(G3y,col=col.cell, pch=16)

# Add a legend to the plot so we know which colours correspond to which cell type
# legend("topleft",fill=c("blue", "green", "purple","orange"),legend=levels(as.factor(G3sampleinfo$groupName)))
legend("topleft",fill=c("blue", "green", "purple","orange"),legend=c("Cellobiose", "Filter paper", "Glucose", "Reed canary grass"), title = expression(bold("Substrate")))

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Now repeat this process with a TPM cutoff

G3TPM <- read.delim("/Users/colleenahern/Documents/EBTNG3/tpm_counts.txt", header=T, check.names = FALSE)
G3TPM$keep <- 0

for (i in 1:nrow(G3TPM)) {
  G3TPM$keep[i] <- ifelse(mean(c(G3TPM$HSOHC[i],G3TPM$HSOHG[i],G3TPM$HSOHH[i])) < 2 && mean(c(G3TPM$HSOHN[i],G3TPM$HWOOB[i],G3TPM$HSOHP[i],G3TPM$HSOHS[i])) < 2 && mean(c(G3TPM$HSOHT[i],G3TPM$HSOHU[i],G3TPM$HSOHW[i],G3TPM$HSOHX[i])) < 2 && mean(c(G3TPM$HSOHY[i],G3TPM$HSOHZ[i],G3TPM$HSONA[i],G3TPM$HSONB[i])) < 2, "remove", "keep") 
}

G3TPM <- G3TPM[G3TPM$keep == "keep",] # reduces genes from 23664 to 11749
G3TPM <- G3TPM[, -17]

G3countdata_tpmfilt <- G3countdata[rownames(G3countdata) %in% G3TPM$Geneid,]

# remake DGEList object
G3y_tpmfilt <- DGEList(G3countdata_tpmfilt)
G3y_tpmfilt

# See what slots are stored in G3y
names(G3y_tpmfilt)

# Library size information is stored in the samples slot
G3y_tpmfilt$samples

group <- paste(G3sampleinfo$groupName)
group

# Convert to factor
group <- factor(group)
group

# Add the group information into the DGEList
G3y_tpmfilt$samples$group <- group
G3y_tpmfilt$samples

G3y_tpmfilt$samples$lib.size
G3y$samples$lib.size.  # the TPM filter has slightly reduced the library size

# Plot the library sizes as a barplot to see whether there are any major discrepancies between the samples 
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(G3y_tpmfilt$samples$lib.size,names=colnames(G3y_tpmfilt),las=2)

# Add a title to the plot
title("Barplot of library sizes")

# Plot the library sizes with and without the TPM cutoff
par(mfrow=c(1,2))
barplot(G3y$samples$lib.size,names=colnames(G3y),las=2)
title("Barplot of library sizes (no TPM cutoff, no TMM normalization)")
barplot(G3y_tpmfilt$samples$lib.size,names=colnames(G3y_tpmfilt),las=2)
title("Barplot of library sizes (TPM cutoff, no TMM normalization)")  # they pretty such look the same

# Get log2 counts per million
G3logcounts_tpmfilt <- cpm(G3y_tpmfilt,log=TRUE)

# Check distributions of samples using boxplots
boxplot(G3logcounts_tpmfilt, xlab="", ylab="Log2 counts per million (TPM filt, no TOMM norm)",las=2)

# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(G3logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

plotMDS(G3y)
levels(as.factor(G3sampleinfo$groupName))

col.cell <- c("blue", "green", "purple","orange")[as.factor(G3sampleinfo$groupName)]
data.frame(G3sampleinfo$groupName,col.cell)

# Redo the MDS with substrate coloring
plotMDS(G3y,col=col.cell, pch=16)

# Add a legend to the plot so we know which colours correspond to which substrate
# legend("topleft",fill=c("blue", "green", "purple","orange"),legend=levels(as.factor(G3sampleinfo$groupName)))
legend("topleft",fill=c("blue", "green", "purple","orange"),legend=c("Cellobiose", "Filter paper", "Glucose", "Reed canary grass"))

# Add a title
title("Substrate type")

# Apply normalisation to DGEList object
G3y <- calcNormFactors(G3y)

# This will update the normalisation factors in the DGEList object (their default values are 1). 
G3y$samples

# If we plot mean difference plots using the plotMDS function for these samples, we should be able to see the 
# composition bias problem. We will use the logcounts, which have been normalised for library size, but not for 
# composition bias.

par(mfrow=c(1,2))
plotMD(G3logcounts,column = 7)
abline(h=0,col="grey")
plotMD(G3logcounts,column = 11)
abline(h=0,col="grey")

# The mean-difference plots show average expression (mean: x-axis) against log-fold-changes (difference: y-axis). 
# Because our DGEList object contains the normalisation factors, if we redo these plots using y, we should see the 
# composition bias problem has been solved.

par(mfrow=c(1,2))
plotMD(G3y,column = 7)
abline(h=0,col="grey")
plotMD(G3y,column = 11)
abline(h=0,col="grey")

G3normlogcounts <- cpm(G3y,log=TRUE)

par(mfrow=c(1,2))

boxplot(G3logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(G3logcounts),col="blue")

boxplot(G3normlogcounts, xlab="", ylab="Log2 counts per million",las=2,main="TMM transformed logCPM")
abline(h=median(G3normlogcounts),col="blue")

## Replot the MDS to see if the TMM normalization changed the MDS plot

plotMDS(G3y)
levels(as.factor(G3sampleinfo$groupName))
col.cell <- c("blue", "green", "purple","orange")[as.factor(G3sampleinfo$groupName)]
data.frame(G3sampleinfo$groupName,col.cell)
plotMDS(G3y,col=col.cell, pch=16)
legend("topleft",fill=c("blue", "green", "purple","orange"),legend=c("Cellobiose", "Filter paper", "Glucose", "Reed canary grass"), title = expression(bold("Substrate")))
title("Substrate type")
