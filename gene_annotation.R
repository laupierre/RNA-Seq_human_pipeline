library(rtracklayer)
library(biomaRt)
library (dplyr)

# hg38 based
system ("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz")

my_obj <- import("gencode.v43.annotation.gtf.gz")
my_obj <- as.data.frame (my_obj)
my_obj <- my_obj[my_obj$type == "transcript", ]
my_obj <- my_obj[ ,c("gene_id", "transcript_id", "gene_type", "gene_name", "hgnc_id")]

values <- unique (my_obj$gene_id)
length (values)
# 62703

ensembl <- useEnsembl(biomart = 'genes', 
                       dataset = 'hsapiens_gene_ensembl',
                       version = 109)

res <- getBM(attributes = c('ensembl_gene_id_version', 'external_gene_name', 'chromosome_name', 'start_position' , 'end_position', 'strand', 'description'),
             filters = 'ensembl_gene_id_version',
             values = values, mart = ensembl)      
dim (res)
# 62656

res$chromosome_name <- paste ("chr", res$chromosome_name, sep="")

my_obj <- merge (my_obj, res, by.x="gene_id", by.y="ensembl_gene_id_version", all.x=TRUE)
dim (my_obj)
# 62703
# 252913

transcripts_num <- my_obj %>% group_by (gene_id) %>% summarise (transcripts_number= n ())

my_obj <- merge (my_obj, transcripts_num, by= "gene_id", all.x=TRUE)
my_obj <- my_obj[ ,c(1:10,12,11)]
head (my_obj)

write.table (my_obj, "gencode.v43.annotation.txt", sep="\t", quote=F, row.names=F)



