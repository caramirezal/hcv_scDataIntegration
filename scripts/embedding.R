## Embedding 
library(Seurat)

## reading data
rna_ann <- readRDS('../data/rna_annotations_miller2019.rds')
rna <- readRDS('../data/rna_seurat_processed_miller2019.rds')
atac <- readRDS('../data/atac_seurat_processed_miller2019.rds')
transfer.anchors <- readRDS('../data/transfer.anchors.rds')

## Embedding
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = rna$orig.ident, 
                                     weight.reduction = atac[["lsi"]])
atac <- AddMetaData(atac, metadata = celltype.predictions)
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "rna", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, 
                           refdata = refdata, 
                           weight.reduction = atac[["lsi"]])
atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac)
coembed <- NormalizeData(coembed)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:50)