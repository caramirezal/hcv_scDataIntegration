## gene activity matrix
library(Seurat)
library(Matrix)

## barcodes
barcodes <- read.table('../data/satpathy2016/barcodes.tsv', 
                       header = TRUE)  
barcodes.s <- as.character(barcodes$Barcodes)
summaryTCells <- readRDS('../data/satpathy2016/scATAC_TME_TCells_SummarizedExperiment.final.rds')

## peaks
peaks <- read.table('../data/satpathy2016/genes.tsv', header = TRUE)
peaks <- as.character(peaks$Feature)

## counts
counts <- readMM('../data/satpathy2016/matrix.mtx')
rownames(counts) <- as.character(peaks)
colnames(counts) <- as.character(barcodes.s)

## create seurat object
counts_seu <- CreateSeuratObject(counts = counts, project = 'tcells', min.cells = 1, min.features = 1)
remove(counts)
counts_seu <- NormalizeData(counts_seu)
counts_seu <- FindVariableFeatures(counts_seu, 
                                     selection.method = 'vst',
                                     nfeatures = 100000)
var_features <- VariableFeatures(counts_seu)
counts_fs <- counts_seu@assays$RNA@counts[var_features, ]

## generate gene activity matrix
activity.matrix <- CreateGeneActivityMatrix(
        peak.matrix = counts_fs, 
        annotation.file = "../data/Homo_sapiens.GRCh38.98.gtf.gz", 
        seq.levels = c(1:22, "X", "Y"), 
        upstream = 2000, 
        verbose = TRUE
)

saveRDS(activity.matrix, '../data/satpathy2016/activity_matrix.RDS')
gc()

