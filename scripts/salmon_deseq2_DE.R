# load packages ----
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(readr)
library(ggplot2)

# salmon results directory ----
dir("/Users/aless/Desktop/Bioinformatics/Year II/Bioinformatics I/virtual machine/share/RNAseq Analysis/data/salmon_results/")

# create a sapmledata data.frame ----
sampledata <- 
  read.csv("SRR_Acc_List.txt", header = FALSE)
colnames(sampledata) <- "runids"

sampledata$group <- rep(c("WT_no_stress", "WT_stress", "Mut_no_stress", "Mut_stress"), each = 3)

rownames(sampledata) <- sampledata$runids
sampledata[1] <- NULL
sampledata

# quant.sf file paths ----
files <- file.path("data/salmon_results/", row.names(sampledata), "quant.sf")
names(files) <- row.names(colData)

# create tx2gene object ----
txdb <- GenomicFeatures::makeTxDbFromGFF("Saccharomyces_cerevisiae.R64-1-1.104.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

# reorder columns of tx2gene ----
tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
# check tx2gene
head(tx2gene)

# import salmon quantification data ----
txi.salmon <- tximport(files = files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, dropInfReps = TRUE)

# check imported data ----
head(txi.salmon$counts)

# check if sample names match between objects ----
identical(x = rownames(colData), y = colnames(txi.salmon$counts))


# DESeq2 pipeline ----
dds <- DESeqDataSetFromTximport(txi = txi.salmon, colData = sampledata, design = ~group)

# filtering not expressed genes ----
keep <- rowSums(counts(dds)) > 1
table(keep)
dds <- dds[keep, ]
nrow(dds)

dds$group <- relevel(x = dds$group, ref = "WT_no_stress")
dds <- DESeq(dds)

# save the results of DE analysis ----
res <- results(dds)

# look at the results ----
summary(res)
res

# PCA ----
rld <- rlog(dds)
DESeq2::plotPCA(rld, ntop = 500, intgroup = 'group') + 
  theme_bw() + labs(title="PCA Plot")

#comparisons ----
r1 <- results(dds, contrast = c("group", "WT_stress", "WT_no_stress"))
r2 <- results(dds, contrast = c("group", "Mut_stress","Mut_no_stress"))
r3 <- results(dds, contrast = c("group", "Mut_stress", "WT_stress"))
r4 <- results(dds, contrast = c("group", "Mut_no_stress","WT_no_stress"))


# Visualization of results: MA plot ----
DESeq2::plotMA(r1, alpha = 0.05, main = "MA Plot for r1", colSig="firebrick2", colLine="firebrick2")
DESeq2::plotMA(r2, alpha = 0.05, main = "MA Plot for r2", colSig="deepskyblue3", colLine="deepskyblue3")
DESeq2::plotMA(r3, alpha = 0.05, main = "MA Plot for r3", colSig="darkmagenta", colLine="darkmagenta")
DESeq2::plotMA(r4, alpha = 0.05, main = "MA Plot for r4", colSig="chartreuse2", colLine="chartreuse2")

#remove genes with NA values 
upregulated_r1 <- r1[!is.na(r1$padj),]
#select genes with adjusted p-values below 0.05
upregulated_r1 <- upregulated_r1[upregulated_r1$padj < 0.05,]
#select genes with absolute log2 fold change above 1 (two-fold change)
upregulated_r1 <- upregulated_r1[upregulated_r1$log2FoldChange > 1,]

# Gene Ontology - GO Enrichment Analysis ----
library("AnnotationDbi")
library("org.Sc.sgd.db")
upregulated_r1$symbol <- mapIds(org.Sc.sgd.db,
                                keys=rownames(upregulated_r1),
                                column="GENENAME",
                                keytype="ENSEMBL",
                                multiVals="first")


library(clusterProfiler)
GO_BP_1 <- enrichGO(upregulated_r1$symbol, OrgDb = "org.Sc.sgd.db", keyType = "GENENAME", ont = "BP")

dotplot(GO_BP_1, title = "Biological Process for r1")


#remove genes with NA values 
upregulated_r2 <- r2[!is.na(r2$padj),]
#select genes with adjusted p-values below 0.05
upregulated_r2 <- upregulated_r2[upregulated_r2$padj < 0.05,]
#select genes with absolute log2 fold change above 1 (two-fold change)
upregulated_r2 <- upregulated_r2[upregulated_r2$log2FoldChange > 1,]

# Gene Ontology - GO Enrichment Analysis ----
upregulated_r2$symbol <- mapIds(org.Sc.sgd.db,
                                keys=rownames(upregulated_r2),
                                column="GENENAME",
                                keytype="ENSEMBL",
                                multiVals="first")


GO_BP_2 <- enrichGO(upregulated_r2$symbol, OrgDb = "org.Sc.sgd.db", keyType = "GENENAME", ont = "BP")
dotplot(GO_BP_2, title = "Biological Process for r2")


#remove genes with NA values 
upregulated_r3 <- r3[!is.na(r3$padj),]
#select genes with adjusted p-values below 0.05
upregulated_r3 <- upregulated_r3[upregulated_r3$padj < 0.05,]
#select genes with absolute log2 fold change above 1 (two-fold change)
upregulated_r3 <- upregulated_r3[upregulated_r3$log2FoldChange > 1,]

# Gene Ontology - GO Enrichment Analysis ----
upregulated_r3$symbol <- mapIds(org.Sc.sgd.db,
                                keys=rownames(upregulated_r3),
                                column="GENENAME",
                                keytype="ENSEMBL",
                                multiVals="first")


GO_BP_3 <- enrichGO(upregulated_r3$symbol, OrgDb = "org.Sc.sgd.db", keyType = "GENENAME", ont = "BP")
dotplot(GO_BP_3, title = "Biological Process for r3")


#remove genes with NA values 
upregulated_r4 <- r4[!is.na(r4$padj),]
#select genes with adjusted p-values below 0.05
upregulated_r4 <- upregulated_r4[upregulated_r4$padj < 0.05,]
#select genes with absolute log2 fold change above 1 (two-fold change)
upregulated_r4 <- upregulated_r4[upregulated_r4$log2FoldChange > 1,]

# Gene Ontology - GO Enrichment Analysis ----
upregulated_r4$symbol <- mapIds(org.Sc.sgd.db,
                                keys=rownames(upregulated_r4),
                                column="GENENAME",
                                keytype="ENSEMBL",
                                multiVals="first")


GO_BP_4 <- enrichGO(upregulated_r4$symbol, OrgDb = "org.Sc.sgd.db", keyType = "GENENAME", ont = "BP")
dotplot(GO_BP_3, title = "Biological Process for r4")


#Step 7: Differential Expression analysis of raw counts from paper authors

#import data
countData <- read.table('GSE98352_DESeq2_raw_counts.tsv')
head(countData)
#change column names
colnames(countData) <- rownames(colData)
head(countData)

# DESeq2 pipeline ----
new_dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampledata, design = ~group)

# filtering not expressed genes ----
new_keep <- rowSums(counts(new_dds)) > 1
table(new_keep)
new_dds <- new_dds[new_keep, ]
nrow(new_dds)

new_dds$group <- relevel(x = new_dds$group, ref = "WT_no_stress")
new_dds <- DESeq(new_dds)

# save the results of DE analysis ----
new_res <- results(new_dds)

# look at the results ----
summary(new_res)
new_res

# PCA ----
new_rld <- rlog(new_dds)
DESeq2::plotPCA(new_rld, ntop = 500, intgroup = 'group') + 
  theme_bw() + labs(title="PCA Plot for Already Processed Data")

# Group Comparisons ----

# r1.2
r1.2 <- results(new_dds, contrast = c("group", "WT_stress", "WT_no_stress"))
DESeq2::plotMA(r1.2, alpha = 0.05, main = "MA Plot for r1.2", colSig="firebrick2", colLine="firebrick2")

# r2.2
r2.2 <- results(new_dds, contrast = c("group", "Mut_stress","Mut_no_stress"))
DESeq2::plotMA(r2.2, alpha = 0.05, main = "MA Plot for r2.2", colSig="deepskyblue3", colLine="deepskyblue3")

# r3.2
r3.2 <- results(new_dds, contrast = c("group", "Mut_stress", "WT_stress"))
DESeq2::plotMA(r3.2, alpha = 0.05, main = "MA Plot for r3.2", colSig="darkmagenta", colLine="darkmagenta")

# r4.2
r4.2 <- results(new_dds, contrast = c("group", "Mut_no_stress","WT_no_stress"))
DESeq2::plotMA(r4.2, alpha = 0.05, main = "MA Plot for r4.2", colSig="chartreuse2", colLine="chartreuse2")

