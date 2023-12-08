# 06/12/2023 Use https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html tutorial to look at Elaina's transcriptomics data

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
BiocManager::install("org.Mm.eg.db", force = TRUE)

G3seqdata <- read.delim("/Users/colleenahern/Documents/EBTNG3/trabundg32/gene_count_matrix.txt")
G3sampleinfo <- read.delim("/Users/colleenahern/Documents/EBTNG3/g3sampleinfo_v2.txt")
head(G3seqdata)
dim(G3seqdata)

# Remove first column from G3seqdata
G3countdata <- G3seqdata[,-(1)]
# Look at the output
head(G3countdata)

# Store GeneID as rownames
rownames(G3countdata) <- G3seqdata[,1]
head(G3countdata)

colnames(G3countdata) # using substr, you can change the characters of the colnames
table(colnames(G3countdata)==G3sampleinfo$X1.libraryName)
colnames(G3countdata) <- gsub(".*\\X", "", colnames(G3countdata))
G3countdata <- G3countdata[, c(11,3,2,8,15,10,7,1,14,4,9,6,13,5,12)]
G3sampleinfo$libraryName <- gsub("-", ".", G3sampleinfo$libraryName)
table(colnames(G3countdata)==G3sampleinfo$libraryName)

colnames(G3sampleinfo) <- gsub(".*\\.", "", colnames(G3sampleinfo))
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

# Adding annotation - I don't think this is correct, so ignore
columns(org.Mm.eg.db)
rownames(G3y$counts) <- gsub(".*\\_", "", rownames(G3y$counts))

G3ann <- select(org.Mm.eg.db,keys=rownames(G3y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
head(G3ann)
table(G3ann$ENTREZID==rownames(G3y$counts))
# G3y$genes <- G3ann

# Filtering lowly expressed genes
# Obtain CPMs
G3CPM <- cpm(G3countdata)
# Have a look at the output
head(G3CPM)

# Which values in G3CPM are greater than 0.5? .3
thresh <- G3CPM > .3
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 10101 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(G3CPM[,1],G3countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(G3CPM[,1],G3countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)   # In this case, a cpm of ~1.2 corresponds to a count of ~10, so use a cpm cutoff of 1.2

# filter the DGEList object
G3y <- G3y[keep, keep.lib.sizes=FALSE]

# Quality control: Now that we have got rid of the lowly expressed genes and have our counts stored in a DGEList object, 
# we can look at a few different plots to check that the data is good quality, and that the samples are as we would expect.

G3y$samples$lib.size

# Plot the library sizes as a barplot to see whether there are any major discrepancies between the samples 
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(G3y$samples$lib.size,names=colnames(G3y),las=2)
# Add a title to the plot
title("Barplot of library sizes")

# we can also adjust the labelling if we want
barplot(G3y$samples$lib.size/1e06, names=colnames(G3y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. 
# Use box plots to check the distribution of the read counts on the log2 scale. 
# Use the cpm function to get log2 counts per million, which are corrected for the different library sizes. 
# The cpm function also adds a small offset to avoid taking log of zero.

# Get log2 counts per million
G3logcounts <- cpm(G3y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(G3logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(G3logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# From the boxplots we see that overall the density distributions of raw log-intensities are not identical but still not very different. 
# If a sample is really far above or below the blue horizontal line we may need to investigate that sample further. 
# Another kind of QC plot that is helpful in checking for dodgy samples is a relative log expression (RLE) plot, which can be generated with plotRLE from the EDASeq package.

plotMDS(G3y)


# Let's set up colour schemes for groupName
# How many cell types and in what order are they stored?
levels(as.factor(G3sampleinfo$groupName))

col.cell <- c("blue", "green", "purple","orange")[as.factor(G3sampleinfo$groupName)]
data.frame(G3sampleinfo$groupName,col.cell)

# Redo the MDS with cell type colouring
plotMDS(G3y,col=col.cell, pch=16)
# Add a legend to the plot so we know which colours correspond to which cell type
legend("bottomright",fill=c("blue", "green", "purple","orange"),legend=levels(as.factor(G3sampleinfo$groupName)))
s# Add a title
title("G3 Substrate type")

# Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
G3var_genes <- apply(G3logcounts, 1, var)
head(G3var_genes)

# Get the gene names for the top 500 most variable genes
G3select_var <- names(sort(G3var_genes, decreasing=TRUE))[1:500]
head(G3select_var)

# Subset logcounts matrix
G3highly_variable_lcpm <- G3logcounts[G3select_var,]
dim(G3highly_variable_lcpm)
head(G3highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("blue", "green", "purple","orange")[as.factor(G3sampleinfo$groupName)]

# Plot the heatmap
heatmap.2(G3highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")


