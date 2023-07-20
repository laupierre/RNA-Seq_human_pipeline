library(DESeq2)
library(tximport)
library(openxlsx)
library (ggplot2)


system ("mkdir ./output")
system ("cp ./projects/log.out ./output")


tx2gene <- read.delim ("gencode.v44.annotation.txt")
tx2gene <- tx2gene[ ,c("transcript_id", "gene_id")]
colnames (tx2gene) <- c("TXNAME", "GENEID")

dir <- paste (paste (getwd (), "projects", sep="/"), "salmon_results", sep="/")
files <- list.files (path=dir, pattern=".*quant.sf", recursive=TRUE)
files <- paste (dir, files, sep="/")
names(files) <- gsub (".*IIT_", "", gsub ("/quant.sf", "", files))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
res <- txi$counts

anno <- read.delim ("gencode.v43.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

res <- merge (res, anno, by.x="row.names", by.y="gene_id", all.x=TRUE)
colnames (res)[1] <- "Geneid"
res <- res[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", res$gene_type, invert=TRUE), ]


write.xlsx (res, "./output/salmon_gene_lengthScaledTPM_counts.xlsx", rowNames=F)
