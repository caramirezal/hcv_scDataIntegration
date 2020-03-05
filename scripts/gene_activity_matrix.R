## gene activity matrix
library(Seurat)
library(Matrix)
library(dplyr)

cat('Reading barcodes\n')
barcodes <- read.table('data/satpathy2016/barcodes.tsv', 
                       header = TRUE)  
barcodes.s <- as.character(barcodes$Barcodes)
summaryTCells <- readRDS('data/satpathy2016/scATAC_TME_TCells_SummarizedExperiment.final.rds')

cat('Reading peaks\n')
peaks <- read.table('data/satpathy2016/genes.tsv', header = TRUE)
peaks <- as.character(peaks$Feature)

cat('Reading peak counts\n')
counts <- readMM('data/satpathy2016/matrix.mtx')
rownames(counts) <- as.character(peaks)
colnames(counts) <- as.character(summaryTCells$Internal_Name)

## create seurat object
<<<<<<< HEAD
#counts_seu <- CreateSeuratObject(counts = counts, project = 'tcells', min.cells = 1, min.features = 1)
#remove(counts)
#counts_seu <- NormalizeData(counts_seu)
#counts_seu <- FindVariableFeatures(counts_seu, 
#                                     selection.method = 'vst',
#                                     nfeatures = 100000)
#var_features <- VariableFeatures(counts_seu)
#counts_fs <- counts_seu@assays$RNA@counts[var_features, ]
=======
counts_seu <- CreateSeuratObject(counts = counts, project = 'tcells', min.cells = 1, min.features = 1)
remove(counts)
counts_seu <- NormalizeData(counts_seu)
counts_seu <- FindVariableFeatures(counts_seu, 
                                     selection.method = 'vst',
                                     nfeatures = 50000)
var_features <- VariableFeatures(counts_seu)
counts_fs <- counts_seu@assays$RNA@counts[var_features, ]
>>>>>>> d9175bde08c5cba9a01edc47f7d6a83deb2e118b

## generate gene activity matrix
#activity.matrix <- CreateGeneActivityMatrix(
#        peak.matrix = counts_fs, 
#        annotation.file = "../data/Homo_sapiens.GRCh38.98.gtf.gz", 
#        seq.levels = c(1:22, "X", "Y"), 
#        upstream = 2000, 
#        verbose = TRUE
#)

#saveRDS(activity.matrix, '../data/satpathy2016/activity_matrix.RDS')
#gc()


#######################################################################
## Processing Satpathy cicero output with seurat 
## to integrate with Mike and Miller data

cat('Getting counts\n')
#summaryTCells <- readRDS('data/satpathy2016/scATAC_TME_TCells_SummarizedExperiment.final.rds')
#peaks <- summaryTCells@assays$data$counts
peaks <- counts
rm(counts)
colnames(peaks) %>% head
rownames(peaks) %>% head()

cat('Setting cicero data\n')
cicero <- readRDS('data/satpathy2016/Log2_Gene_Activity_TME_All_SummarizedExperiment.final.rds')
gene_activity <- cicero@assays$data$logGA
gene_activity <- gene_activity[, colnames(peaks)]
colnames(gene_activity) %>% head()
rownames(gene_activity) %>% head()

cat('dim peaks\n')
dim(peaks)

cat('dim gene activity\n')
dim(gene_activity)

#cat('Dimensions check\n')
#( colnames(peaks) == colnames(gene_activity) ) %>% sum()

#cat('Creating seurat object\n')
#satpathy <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "Satpathy")
## This line is for getting annotations from SummaryTCell file object
#satpathy$'cell_type' <- summaryTCells$T_Cell_Cluster
#satpathy[["ACTIVITY"]] <- CreateAssayObject(counts = gene_activity)
#satpathy <- subset(satpathy, subset = nCount_ATAC > 5000)
#satpathy$tech <- "atac"

#cat('Normalization\n')
#DefaultAssay(satpathy) <- "ACTIVITY"
#satpathy <- FindVariableFeatures(satpathy)
#satpathy <- NormalizeData(satpathy)
#satpathy <- ScaleData(satpathy)

#cat('Running LSI\n')
#DefaultAssay(satpathy) <- "ATAC"
#VariableFeatures(satpathy) <- names(which(Matrix::rowSums(satpathy) > 100))
#satpathy <- RunLSI(satpathy, n = 50, scale.max = NULL)
#satpathy <- RunUMAP(satpathy, reduction = "lsi", dims = 1:50)

#cat('Adding annotations\n')

#cat('Saving data\n')
#saveRDS(satpathy, 'data/satpathy_processed_tcells_seurat.rds')

