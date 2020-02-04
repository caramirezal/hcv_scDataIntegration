library(Seurat)

rna <- readRDS('../data/rna_seurat_processed_miller2019.rds')
atac <- readRDS('../data/atac_seurat_processed_miller2019.rds')

transfer.anchors <- FindTransferAnchors(reference = rna, 
                                        query = atac, 
                                        features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "rna", 
                                        query.assay = "activity", 
                                        reduction = "cca")