# 11/14/2023
# final proofread of complete DESeq analysis on G3 count matrix 
# Manuscript heatmap: all genes with CAZymes labeled with dbCAN annotation; tpm filter applied)

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
vary_43 <- read.delim("/Users/colleenahern/Documents/EBTNG3/DGE_summary_1304689.txt", header=T, check.names = FALSE) #Daniel
vary_43 <- vary_43[,c(1,2,3,8,9,17,18)] 
vary_43[vary_43=="  NA"] <- NA   

vary_43$rem <- 0
for (i in 1:nrow(vary_43)) {
  vary_43$rem[i] <- ifelse(is.na(vary_43$RCG_glucose_Padj[i]) == TRUE && is.na(vary_43$FP_glucose_Padj[i]) == TRUE && is.na(vary_43$cellobiose_glucose_Padj[i]) == TRUE, "REMOVE", "KEEP")
}

vary_43 <- vary_43[(vary_43$rem=="KEEP"),]
vary_43 <- vary_43[,-8]

# Merge tpm matrix with vary_43 matrix - to make a matrix that contains only genes that pass both of those filters
a <- merge(vary_43,tpm,"Geneid")
vary_43 <- a[,c(1:7)]

# tpm$huh <- 0
# for (i in 1:nrow(tpm)) {
# tpm$huh[i] <- ifelse(tpm$Geneid[i] %in% vary_43e$Geneid == TRUE, "IN", "OUT")
# }

vary_43$RCG_glucose_Padj <- as.numeric(gsub('E', 'e', vary_43$RCG_glucose_Padj))
vary_43$cellobiose_glucose_Padj <- as.numeric(gsub('E', 'e', vary_43$cellobiose_glucose_Padj))
vary_43$FP_glucose_Padj <- as.numeric(gsub('E', 'e', vary_43$FP_glucose_Padj))

# Create a minimum p value column that takes the minimum p value from any of the three comparitive conditions for each gene
# Then sort the matrix from lowest padjmin to largest
vary_43$padjmin <- apply(vary_43[,c("cellobiose_glucose_Padj","FP_glucose_Padj","RCG_glucose_Padj")],1,min, na.rm = TRUE) #Daniel
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]
colnames(vary_43) <- c("geneID", "Cellobiose", "cellobiose_glucose_Padj", "FP", "FP_glucose_Padj", "RCG", "RCG_glucose_Padj", "padjmin")

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
fullcazgenelist <- merge(vary_43,cazcomb2, "proteinID", all.x = TRUE)
fullcazgenelist <- fullcazgenelist[,c(2,1,10,11,12,13,3:9)]

fullcazgenelist$combname <- 0

for (i in 1:nrow(fullcazgenelist)) {
  fullcazgenelist$combname[i] <- ifelse(is.na(fullcazgenelist$CAZ_name[i]) == TRUE, fullcazgenelist$name[i], fullcazgenelist$CAZ_name[i])
}

fullcazgenelist <- fullcazgenelist[,c(1,2,14,7:13)]
fullcazgenelist <- fullcazgenelist[order(fullcazgenelist$padjmin,decreasing=FALSE),]

newnames <- lapply(
  colnames(fullcazgenelist[,c(4,6,8)]),
  function(x) bquote(.(x)))
head(fullcazgenelist)
nrow(fullcazgenelist)
#quantile(rowSums(vary_43))
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

#ann_colors = list("GH Family"=c("GH1"="#00CC99","GH10"="#984EA3" ,"GH11"="#FF7F00","GH114"="#FFFF33","GH115"="#4DAF4A","GH13_28"="#DECBE4","GH16"="#996633","GH18"="#999999","GH2"="#E5D8BD","GH24"="#80B1D3","GH26"="#FDB462","GH3"="#B3DE69","GH30_5"="#FCCDE5","GH31"="#E41A1C","GH32"="#CCEBC5","GH39"="#F781BF","GH43"="#FB8072","GH45"="#FFCC00","GH48"="#7FC97F","GH5_1"="#666666","GH5_4"="#999999","GH5_5"="#F0027F","GH5_7"="#FF9933","GH6"="#FB8072","GH8"="#386CB0","GH9"="#FDCDAC","GH95"="#984EA3"),"Multiple Annotations"=c("Y"="#1F78B4","N"="#A6CEE3"))

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
#GH_regulation_heatmap<-pheatmap(data_matrix_44, annotation_row = my_gene_col_final_S4timecourse, annotation_colors = ann_colors, col=my_palette, cellwidth=55,cellheight=8, cluster_rows=TRUE, cluster_cols=F, 
#  fontsize=14,show_rownames=F, fontsize_row=8, fontsize_col=15, labels_col = as.expression(newnames), border_color=NA,breaks=seq(-6, 6, length.out=300))

g3_regulation_heatmap<-pheatmap(data.matrix(vary_44[,c(4,6,8)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(vary_44[,3]), border_color=NA,breaks=seq(-4, 6, length.out=300))



save_pheatmap_png <- function(x, filename, width=1200, height=2000, res = 1500) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}



save_pheatmap_png(g3_regulation_heatmap, "g3_regulation_heatmap.png")


# CAZYmes only
genelist <- read.delim("/Users/colleenahern/Documents/EBTNG3/Neocon1_GeneCatalog_20220615_2_edit.txt", header = T, check.names = FALSE)
genelist <- genelist[-c(2)]
cazlist <- read.delim("/Users/colleenahern/Documents/EBTNG3/Updated_CAZyme_list.txt", header = T, check.names = FALSE)
cazlist$proteinId <- gsub("^.{0,8}", "", cazlist$proteinId)
for (i in 1:nrow(cazlist)) {
  cazlist$name[i] <- ifelse(cazlist$annotation[i]=="DOC2", "Dockerin", cazlist$name[i])
}

library(tidyverse)
library(plyr)
cazcomb <- ddply(cazlist,.(proteinId,annotation,name),nrow)
cazcomb$name2 <- 0
cazcomb$annotation2 <- 0
for (i in 1:nrow(cazcomb)) {
  cazcomb$name2[i] <- ifelse(cazcomb$V1[i] > 1, paste(cazcomb$name[i], cazcomb$V1[i], sep = " x"), cazcomb$name[i])
  cazcomb$annotation2[i] <- ifelse(cazcomb$V1[i] > 1, paste(cazcomb$annotation[i], cazcomb$V1[i], sep = " x"), cazcomb$annotation[i])
}
cazcomb2 <- cazcomb %>%
  group_by(proteinId) %>%
  mutate(annotation2 = paste0(annotation2, collapse = ", ")) %>%
  mutate(name2 = paste0(name2, collapse = ", ")
  ) 
cazcomb2 <- cazcomb2[!duplicated(cazcomb2$proteinId),]
cazcomb2 <- cazcomb2[, -c(2,3,4)]
cazcomb2 <- cazcomb2[, c(1,3,2)]
write.csv(cazcomb2, "/Users/colleenahern/Documents/EBTNG3/cazcomb.csv", row.names = FALSE)


# cazcomb3 <- cazlist %>%
#   group_by(proteinId) %>%
#   mutate(annotation = paste0(annotation, collapse = ", ")) %>%
#   mutate(name = paste0(name, collapse = ", ")
#   ) 
# cazcomb3 <- cazcomb3[!duplicated(cazcomb3$proteinId),]

# # compare mine and Tejas's
# tejcaz <- read.delim("/Users/colleenahern/Downloads/G3_uniqueCAZymes_by_proteinID_mycocosm_.txt", header = T, check.names = FALSE)
# cazcomb3$comp <- 0 
# for (i in 1:nrow(cazcomb3)) {
#   cazcomb3$comp[i] <- ifelse(cazcomb3$proteinId[i]==tejcaz$proteinId[i] && cazcomb3$annotation[i]==tejcaz$annotation[i], TRUE, FALSE)
# }


colnames(cazcomb2) <- c("proteinID", "CAZ_annotation", "CAZ_name")
cazlist2 <- merge(cazcomb2, genelist, "proteinID")
cazlist2 <- cazlist2[,c(1,4,6,2,3,5)]


vary_43 <- rename(vary_43, "geneID" = "Geneid")
vary_43 <- merge(vary_43, cazlist2, "geneID")
vary_43$padjmin <- apply(vary_43[,c("cellobiose_glucose_Padj","FP_glucose_Padj","RCG_glucose_Padj")],1,min, na.rm = TRUE) #Daniel
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]
vary_43

newnames <- lapply(
  colnames(vary_43[,c(2,4,6)]),
  function(x) bquote(.(x)))
head(vary_43)
nrow(vary_43)
#quantile(rowSums(vary_43))
rownum_from <- 1
rownum_to <- 100

vary_44 <- vary_43[c(rownum_from:rownum_to),]
data_matrix_44 <- vary_44[,c(1,2,4,6,12)]
data_matrix_44

ann_cols <- data_matrix_44$CAZ_name
# install.packages("RColorBrewer")
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

#ann_colors = list("GH Family"=c("GH1"="#00CC99","GH10"="#984EA3" ,"GH11"="#FF7F00","GH114"="#FFFF33","GH115"="#4DAF4A","GH13_28"="#DECBE4","GH16"="#996633","GH18"="#999999","GH2"="#E5D8BD","GH24"="#80B1D3","GH26"="#FDB462","GH3"="#B3DE69","GH30_5"="#FCCDE5","GH31"="#E41A1C","GH32"="#CCEBC5","GH39"="#F781BF","GH43"="#FB8072","GH45"="#FFCC00","GH48"="#7FC97F","GH5_1"="#666666","GH5_4"="#999999","GH5_5"="#F0027F","GH5_7"="#FF9933","GH6"="#FB8072","GH8"="#386CB0","GH9"="#FDCDAC","GH95"="#984EA3"),"Multiple Annotations"=c("Y"="#1F78B4","N"="#A6CEE3"))

# install.packages("pheatmap")
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
#GH_regulation_heatmap<-pheatmap(data_matrix_44, annotation_row = my_gene_col_final_S4timecourse, annotation_colors = ann_colors, col=my_palette, cellwidth=55,cellheight=8, cluster_rows=TRUE, cluster_cols=F, 
#  fontsize=14,show_rownames=F, fontsize_row=8, fontsize_col=15, labels_col = as.expression(newnames), border_color=NA,breaks=seq(-6, 6, length.out=300))
#g3_regulation_heatmap<-pheatmap(data_matrix_44, col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
# fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(my_gene_col_g3[,"product"]), border_color=NA,breaks=seq(-4, 6, length.out=300))

g3_regulation_heatmap<-pheatmap(data.matrix(data_matrix_44[,c(2,3,4)]), col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.expression(newnames), labels_row = as.expression(data_matrix_44[,5]), border_color=NA,breaks=seq(-4, 6, length.out=300))



save_pheatmap_png <- function(x, filename, width=1200, height=2000, res = 1500) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

