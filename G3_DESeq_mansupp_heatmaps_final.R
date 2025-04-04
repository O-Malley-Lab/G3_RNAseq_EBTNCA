# Colleen Ahern
# Final manuscript heatmap from my DESeq2 analysis (G3_DESeq_analysis03282025.R)
# all genes with CAZymes labeled with dbCAN annotation; tpm filter applied and padjmin ranked
# Also includes supplement heatmaps: all genes and CAZyme-only

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("pasilla")

library(tximport)
library(readr)
library(tximportData)
library(tidyverse)
library(pasilla)

# Heatmap for all genes

# Apply tpm cutoff to the gene count matrix
# Import tpm gene count matrix
tpm <- read.delim("/Users/colleenahern/Documents/EBTNG3/tpm_counts.txt", header=T, check.names = FALSE)

# Import/reformat sample information 
asm <- read.delim("/Users/colleenahern/Documents/EBTNG3/asm_notes2", header=T, check.names = FALSE)
asm <- asm[,c(1,5)]
colnames(asm) <- c("libraryName", "sampleName")
asm$sampleName <- gsub("G3_", "", asm$sampleName)
colnames(tpm)[2:16] == asm$libraryName

a <-  tpm[asm$libraryName]
a
tpm <- tpm[1]
tpm[2:16] <- a

# Make sure the sample names in the count matrix and sample info match
colnames(tpm)[2:16] == asm$libraryName

# Replace library names with sample names
colnames(tpm)[2:16] <- asm$sampleName

tpm$rem <- 0

# Apply tpm cutoff: the average TPM of the  biological replicates for any of the four substrates must be >2
for (i in 1:nrow(tpm)) {
  tpm$rem[i] <- ifelse(mean(c(tpm$RCG2[i],tpm$RCG3[i],tpm$RCG4[i])) < 2 && mean(c(tpm$FP1[i],tpm$FP2[i],tpm$FP3[i],tpm$FP4[i])) < 2 && mean(c(tpm$glucose1[i],tpm$glucose2[i],tpm$glucose3[i],tpm$glucose4[i])) < 2 && mean(c(tpm$cellobiose1[i],tpm$cellobiose2[i],tpm$cellobiose3[i],tpm$cellobiose4[i])) < 2, "REMOVE", "KEEP")
}

tpm <- tpm[(tpm$rem=="KEEP"),] # reduces from 23664 genes to 11749
tpm <- tpm[,-c(17)]

# Import list of differentially regulated genes for heatmap
# Remove genes with a p-value = NA for all three comparitive conditions

vary_43 <- read.delim("/Users/colleenahern/Documents/EBTNG3/Finalized_codes/Tables/CBA_DESeq_final_03282025.txt", header=T, check.names = FALSE) # My DESeq
vary_43 <- vary_43[,c(1,5,6,11,12,17,18)] 

vary_43$rem <- 0
for (i in 1:nrow(vary_43)) {
  vary_43$rem[i] <- ifelse(is.na(vary_43$RCG_v_glucose_padj[i]) == TRUE && is.na(vary_43$FP_v_glucose_padj[i]) == TRUE && is.na(vary_43$cellobiose_v_glucose_padj[i]) == TRUE, "REMOVE", "KEEP")
}

vary_43 <- vary_43[(vary_43$rem=="KEEP"),] # reduces from 19384 to 18247
vary_43 <- vary_43[,-8]
names(vary_43)[names(vary_43) == 'GeneID'] <- 'Geneid'

# Merge tpm matrix with vary_43 matrix - to make a matrix that contains only genes that pass both of those filters
a <- merge(vary_43,tpm,"Geneid") # 11,736 genes
vary_43 <- a[,c(1:7)]

# Create a minimum p value column that takes the minimum p value from any of the three comparative conditions for each gene
# Then sort the matrix from lowest padjmin to largest
vary_43$padjmin <- apply(vary_43[,c("RCG_v_glucose_padj","FP_v_glucose_padj","cellobiose_v_glucose_padj")],1,min, na.rm = TRUE)
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]
colnames(vary_43) <- c("geneID", "RCG", "RCG_glucose_padj", "FP", "FP_glucose_padj", "Cellobiose", "Cellobiose_glucose_padj", "padjmin")

# Import list of JGI gene annotations
my_gene_col_44 <- read.delim("/Users/colleenahern/Documents/EBTNG3/Neocon1_GeneCatalog_20220615_2_edit.txt", check.names = FALSE)
my_gene_col_44 <- my_gene_col_44[,c(1,3,4)]
my_gene_col_44

# Add annotations to the data frame and reorder padjmin
vary_43 <- merge(vary_43, my_gene_col_44, "geneID")
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]

for (i in 1:nrow(vary_43)) {
  vary_43$name[i] <- paste(vary_43$product_name[i], " ", "(proteinID = ", vary_43$proteinID[i], ")", sep = "")
}

write_csv(vary_43, "/Users/colleenahern/Documents/EBTNG3/Finalized_codes/Tables/vary_43.csv")
vary_43 <- read_csv("/Users/colleenahern/Documents/EBTNG3/Finalized_codes/Tables/vary_43.csv")

# replace the JGI name of any CAZyme with the CAZyme name
cazlist <- read.delim("/Users/colleenahern/Documents/EBTNG3/Updated_CAZyme_list.txt", header = T, check.names = FALSE)
cazlist$proteinId <- gsub("^.{0,8}", "", cazlist$proteinId)
for (i in 1:nrow(cazlist)) {
  cazlist$name[i] <- ifelse(cazlist$annotation[i]=="DOC2", "Dockerin", cazlist$name[i])
}

library(tidyverse)
library(plyr)
library(dplyr)

cazcomb <- ddply(cazlist,.(proteinId,annotation,name),nrow)
cazcomb$name2 <- 0
cazcomb$annotation2 <- 0

for (i in 1:nrow(cazcomb)) {
  cazcomb$name2[i] <- ifelse(cazcomb$V1[i] > 1, paste(cazcomb$name[i], cazcomb$V1[i], sep = " x"), cazcomb$name[i])
  cazcomb$annotation2[i] <- ifelse(cazcomb$V1[i] > 1, paste(cazcomb$annotation[i], cazcomb$V1[i], sep = " x"), cazcomb$annotation[i])
}

detach(package:plyr)

cazcomb2 <- cazcomb %>%
  group_by(proteinId) %>%
  mutate(annotation2 = paste0(annotation2, collapse = ", ")) %>%
  mutate(name2 = paste0(name2, collapse = ", "))

cazcomb2 <- cazcomb2[!duplicated(cazcomb2$proteinId),]
cazcomb2 <- cazcomb2[, -c(2,3,4)]
cazcomb2 <- cazcomb2[, c(1,3,2)]

for (i in 1:nrow(cazcomb2)) {
  cazcomb2$name3[i] <- paste("(proteinID = ", cazcomb2$proteinId[i], ")", sep = "")
  cazcomb2$namefin[i] <- paste("CAZyme: ", cazcomb2$name2[i], cazcomb2$name3[i], sep = " ")
}

cazcomb2 <- cazcomb2[, -c(3,4)]
colnames(cazcomb2) <- c("proteinID", "CAZ_annotation", "CAZ_name")
write_csv(cazcomb2, "/Users/colleenahern/Documents/EBTNG3/Finalized_codes/Tables/cazcomb2.csv")
cazcomb2 <- read_csv("/Users/colleenahern/Documents/EBTNG3/Finalized_codes/Tables/cazcomb2.csv")
cazcomb2 <- cazcomb2[,-1]
fullcazgenelist <- merge(vary_43,cazcomb2, "proteinID", all.x = TRUE)
fullcazgenelist <- fullcazgenelist[,c(2,1,10,11,12,13,3:9)]

fullcazgenelist$combname <- 0

for (i in 1:nrow(fullcazgenelist)) {
  fullcazgenelist$combname[i] <- ifelse(is.na(fullcazgenelist$CAZ_name[i]) == TRUE, fullcazgenelist$name[i], fullcazgenelist$CAZ_name[i])
}

fullcazgenelist <- fullcazgenelist[,c(1,2,14,7:13)]
fullcazgenelist <- fullcazgenelist[order(fullcazgenelist$padjmin,decreasing=FALSE),]
fullcazgenelist <- fullcazgenelist[,c(1:3,8,9,6,7,4,5,10)]

newnames <- lapply(
  colnames(fullcazgenelist[,c(4,6,8)]),
  function(x) bquote(.(x)))
newnames
head(fullcazgenelist)
nrow(fullcazgenelist)
rownum_from <- 1
rownum_to <- 100

vary_44 <- fullcazgenelist[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,3]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

g3_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(4,6,8)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,3]), border_color=NA,breaks=seq(-4, 6, length.out=300))


#####################################################################################################################################
# Same heatmap but ranking based on the lowest p-value taken from the condition of reed canary grass (RCG) compared to glucose

# vary_43 from line 69
# Then sort the matrix from lowest padj to largest of the RCG vs. glucose padj values

vary_43_rcg <- vary_43[order(vary_43$RCG_v_glucose_padj,decreasing=FALSE),]
colnames(vary_43_rcg) <- c("geneID", "RCG", "RCG_glucose_padj", "FP", "FP_glucose_padj", "Cellobiose", "Cellobiose_glucose_padj")

# Import list of JGI gene annotations
my_gene_col_44 <- read.delim("/Users/colleenahern/Documents/EBTNG3/Neocon1_GeneCatalog_20220615_2_edit.txt", check.names = FALSE)
my_gene_col_44 <- my_gene_col_44[,c(1,3,4)]
my_gene_col_44

# Add annotations to the data frame and reorder padjmin
vary_43_rcg <- merge(vary_43_rcg, my_gene_col_44, "geneID")
vary_43_rcg <- vary_43_rcg[order(vary_43_rcg$RCG_glucose_padj,decreasing=FALSE),]

for (i in 1:nrow(vary_43_rcg)) {
  vary_43_rcg$name[i] <- paste(vary_43_rcg$product_name[i], " ", "(proteinID = ", vary_43_rcg$proteinID[i], ")", sep = "")
}

# replace the JGI name of any CAZyme with the CAZyme name
cazlist <- read.delim("/Users/colleenahern/Documents/EBTNG3/Updated_CAZyme_list.txt", header = T, check.names = FALSE)
cazlist$proteinId <- gsub("^.{0,8}", "", cazlist$proteinId)
for (i in 1:nrow(cazlist)) {
  cazlist$name[i] <- ifelse(cazlist$annotation[i]=="DOC2", "Dockerin", cazlist$name[i])
}

library(tidyverse)
library(plyr)
library(dplyr)

cazcomb <- ddply(cazlist,.(proteinId,annotation,name),nrow)
cazcomb$name2 <- 0
cazcomb$annotation2 <- 0

for (i in 1:nrow(cazcomb)) {
  cazcomb$name2[i] <- ifelse(cazcomb$V1[i] > 1, paste(cazcomb$name[i], cazcomb$V1[i], sep = " x"), cazcomb$name[i])
  cazcomb$annotation2[i] <- ifelse(cazcomb$V1[i] > 1, paste(cazcomb$annotation[i], cazcomb$V1[i], sep = " x"), cazcomb$annotation[i])
}

detach(package:plyr)

cazcomb2 <- cazcomb %>%
  group_by(proteinId) %>%
  mutate(annotation2 = paste0(annotation2, collapse = ", ")) %>%
  mutate(name2 = paste0(name2, collapse = ", "))

cazcomb2 <- cazcomb2[!duplicated(cazcomb2$proteinId),]
cazcomb2 <- cazcomb2[, -c(2,3,4)]
cazcomb2 <- cazcomb2[, c(1,3,2)]

for (i in 1:nrow(cazcomb2)) {
  cazcomb2$name3[i] <- paste("(proteinID = ", cazcomb2$proteinId[i], ")", sep = "")
  cazcomb2$namefin[i] <- paste("CAZyme: ", cazcomb2$name2[i], cazcomb2$name3[i], sep = " ")
}

cazcomb2 <- cazcomb2[, -c(3,4)]
colnames(cazcomb2) <- c("proteinID", "CAZ_annotation", "CAZ_name")
fullcazgenelist_rcg <- merge(vary_43_rcg,cazcomb2, "proteinID", all.x = TRUE)
fullcazgenelist_rcg <- fullcazgenelist_rcg[,c(2,1,9:12,7,8,5,6,3,4)]

fullcazgenelist_rcg$combname <- 0

for (i in 1:nrow(fullcazgenelist_rcg)) {
  fullcazgenelist_rcg$combname[i] <- ifelse(is.na(fullcazgenelist_rcg$CAZ_name[i]) == TRUE, fullcazgenelist_rcg$name[i], fullcazgenelist_rcg$CAZ_name[i])
}

fullcazgenelist_rcg <- fullcazgenelist_rcg[,c(1,2,13,7:12)]
fullcazgenelist_rcg <- fullcazgenelist_rcg[order(fullcazgenelist_rcg$RCG_glucose_padj,decreasing=FALSE),]

newnames <- lapply(
  colnames(fullcazgenelist_rcg[,c(4,6,8)]),
  function(x) bquote(.(x)))
newnames
head(fullcazgenelist_rcg)
nrow(fullcazgenelist_rcg)
rownum_from <- 1
rownum_to <- 100

vary_44 <- fullcazgenelist_rcg[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,3]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

g3_rcg_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(4,6,8)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,3]), border_color=NA,breaks=seq(-4, 6, length.out=300))


####################################################################################################################################
# Same heatmap but ranking based on the lowest p-value taken from the condition of cellobiose compared to glucose

fullcazgenelist_cb <- fullcazgenelist_rcg[order(fullcazgenelist_rcg$Cellobiose_glucose_padj,decreasing=FALSE),]

newnames <- lapply(
  colnames(fullcazgenelist_cb[,c(4,6,8)]),
  function(x) bquote(.(x)))
newnames
head(fullcazgenelist_cb)
nrow(fullcazgenelist_cb)
rownum_from <- 1
rownum_to <- 100

vary_44 <- fullcazgenelist_cb[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,3]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

g3_cb_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(4,6,8)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                    fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,3]), border_color=NA,breaks=seq(-4, 6, length.out=300))


####################################################################################################################################
# Same heatmap but ranking based on the lowest p-value taken from the condition of cellobiose compared to glucose

fullcazgenelist_fp <- fullcazgenelist_rcg[order(fullcazgenelist_rcg$FP_glucose_padj,decreasing=FALSE),]

newnames <- lapply(
  colnames(fullcazgenelist_fp[,c(4,6,8)]),
  function(x) bquote(.(x)))
newnames
head(fullcazgenelist_fp)
nrow(fullcazgenelist_fp)
rownum_from <- 1
rownum_to <- 100

vary_44 <- fullcazgenelist_fp[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,3]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

g3_fp_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(4,6,8)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                   fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,3]), border_color=NA,breaks=seq(-4, 6, length.out=300))


#####################################################################################################################################
# Heatmap for CAZYmes only

cazcomb2 <- read_csv("/Users/colleenahern/Documents/EBTNG3/Finalized_codes/Tables/cazcomb2.csv")
vary_43 <- read_csv("/Users/colleenahern/Documents/EBTNG3/Finalized_codes/Tables/vary_43.csv")

cazgenelist <- merge(vary_43,cazcomb2, "proteinID")
cazgenelist$CAZ_name2 <- gsub("CAZyme: ", "", cazgenelist$CAZ_name)
cazgenelist <- cazgenelist[,c(2,1,10:14,7,8,5,6,3,4,9)]

cazgenelist <- cazgenelist[order(cazgenelist$padjmin,decreasing=FALSE),]

newnames <- lapply(
  colnames(cazgenelist[,c(8,10,12)]),
  function(x) bquote(.(x)))
newnames
head(cazgenelist)
nrow(cazgenelist)
rownum_from <- 1
rownum_to <- 100

vary_44 <- cazgenelist[c(rownum_from:rownum_to),]

ann_cols <- vary_44[,7]
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

g3_caz_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(8,10,12)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,7]), border_color=NA,breaks=seq(-4, 6, length.out=300))



















