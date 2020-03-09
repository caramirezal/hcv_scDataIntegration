## gene annotation
library(biomaRt)

#writeLines(rownames(rna), '../data/genes_ensg_eltahla.txt')
#saveRDS(rna, '../data/eltahla2016/rna_raw.rds')
rna <- readRDS('../data/eltahla2016/rna_raw.rds')

rownames(rna) %>% head()

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl",
                       host = "www.ensembl.org")
biomart <- getBM(attributes=c('ensembl_gene_id',
                              "hgnc_symbol"), 
                 filters = 'ensembl_gene_id', 
                 values = unique(rownames(rna)), 
                 mart = ensembl)



## dropping duplicated symbols
biomart.s <- biomart.s[!duplicated(biomart$hgnc_symbol), ]
biomart.s <- biomart.s[!duplicated(biomart.s$ensembl_gene_id), ]
length(biomart.s$ensembl_gene_id); length(unique(biomart.s$ensembl_gene_id))
length(biomart.s$hgnc_symbol); length(unique(biomart.s$hgnc_symbol))

## dropping genes in ensg with no hgnc annotation
rna.s <- rna[ rownames(rna) %in% biomart.s$ensembl_gene_id, ]
rna.s <- rna.s[!duplicated(rownames(rna.s)),]
biomart.s <- biomart.s[ biomart.s$ensembl_gene_id %in% rownames(rna.s),]

ord <- match(rownames(rna.s), biomart.s$ensembl_gene_id)
dim(biomart.s); dim(rna.s)

## relabeling gene names
rownames(rna.s) <- biomart.s$hgnc_symbol
RNA
saveRDS(rna.s, '../data/eltahla2016/rna.rds')
