library (openxlsx)
library (DESeq2)
library (ggplot2)

system ("mkdir ./output")
system ("cp ./projects/log.out ./output")

anno <- read.delim ("gencode.v44.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)
dim (anno)
# 62703

a <- read.delim ("./projects/star_results/subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
length (intersect (a$Geneid, anno$gene_id))
# 62703

# remove RNAs
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
# remove Y paralogs
a <- a[grep ("PAR_Y", a$Geneid, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub (".*IIT_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"
  
write.xlsx (a, "./output/star_gene_raw_counts.xlsx", rowNames=F)



