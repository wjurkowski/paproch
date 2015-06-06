options(digits = 3, quote = F, scientific = F, stringsAsFactors = F, echo=F, header=F)
args <- commandArgs(trailingOnly = TRUE)

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library("biomaRt")

# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

#read data
v <- read.table(args[1], header=F)

# query biomart
results <- getBM(attributes = c("ensembl_transcript_id", "uniprot_swissprot", "uniprot_swissprot_accession", "refseq_mrna"), filters = "ensembl_transcript_id", values = v, mart = mart)
results[,3]

# full results
#  ensembl_transcript_id uniprot_swissprot uniprot_swissprot_accession
#1       ENST00000380152       BRCA2_HUMAN                      P51587