## Data preprocessing of HCV project data
library(Seurat)
library(Signac)
library(tidyverse)
library(liger)

## Input: A seurat object
## Performs normalization, scaling, pca and umap 
## Output: A processed seurat object
st_workflow <- function(
        seurat_object,
        n_features = 3000,
        n_pca_dims = 15
){
        cat('Normalizing and finding variable features\n')
        seurat.p <- NormalizeData(seurat_object) %>%
                FindVariableFeatures(selection.method = 'vst',
                                     nfeatures = n_features)
        cat('Scaling and projection\n')
        seurat.p <- ScaleData(seurat.p, 
                              verbose = FALSE) %>% 
                RunPCA(npcs = n_pca_dims, 
                       verbose = FALSE) %>%
                RunUMAP(reduction = "pca", 
                        dims = 1:n_pca_dims) 
        return(seurat.p)
}

## Maike Hoffman data preprocessing
himmer <- read.csv('../data/nina_thimme_raw_counts.csv', 
                   header = TRUE)
r_names <- himmer$X
himmer <- select(himmer, -X)
himmer_mtx <- apply(himmer, 2, as.numeric)
rownames(himmer_mtx) <- gsub('__chr.*', '', r_names) 

## Seurat object construction
himmer_seu <- CreateSeuratObject(counts = himmer_mtx, project = 'himmer', min.cells = 1, assay = 'rna')
himmer_seu <- st_workflow(himmer_seu, n_features = 1000, n_pca_dims = 100)
saveRDS(himmer_seu, '../data/himmer_seu.rds')

## Splitting Hoffman data for liger pipeline
patients <- himmer_seu$orig.ident
himmer_df <- as.data.frame(t(himmer_mtx))
himmer_sp <- split(himmer_df, f = patients) 
himmer_sp <- lapply(himmer_sp, function(x) t(x) )
saveRDS(himmer_sp, '../data/himmer_splitted.rds')

## Standardizing gene names
## Gene names should be in upper case in order to align with human data.
## Run once in order to create ../data/miller2019_upper_Case/
dir.create('../data/miller2019_upper_Case')
file.copy('../data/miller2019/matrix.mtx', '../data/miller2019_upper_Case/')
file.copy('../data/miller2019/barcodes.tsv', '../data/miller2019_upper_Case/')
genes_lc <- read.table('../data/miller2019/genes.tsv', header = FALSE)
genes_lc$V2 <- toupper(genes_lc$V2)
write_tsv(genes_lc, '../data/miller2019_upper_Case/genes.tsv', col_names = FALSE)

### Seurat standard preprocessing
miller_mtx <- Read10X('../data/miller2019_upper_Case/')
miller <- CreateSeuratObject(counts = miller_mtx, project = 'lcmv', assay = 'RNA', min.cells = 1, min.features = 200)
miller <- NormalizeData(miller)
miller <- FindVariableFeatures(miller, selection.method = 'vst', nfeatures = 3000)
miller <- ScaleData(miller, verbose = FALSE)
miller <- RunPCA(miller, npcs = 30, verbose = FALSE)
miller <- RunUMAP(miller, reduction = "pca", dims = 1:30)
miller <- FindNeighbors(miller, dims = 1:30, verbose = FALSE)
miller <- FindClusters(miller, resolution= 0.1, verbose = FALSE)
saveRDS(miller, '../data/miller_seu.rds')

## Miller Cluster annotation
miller$orig.ident <- miller$seurat_clusters
miller_ann <- data.frame('cluster' = miller$orig.ident,
                         'cell_type' = plyr::mapvalues(x = miller$orig.ident,
                                                       from = as.factor(c(3,0,4,2,1)), 
                                                       to = c('Proliferating', 'Effector', 'Naive',
                                                              'Progenitor Ex', 'Terminally Ex')
                         )
)
miller <- AddMetaData(miller, metadata = miller_ann)
saveRDS(miller_ann, '../data/miller_annotations.rds')

## checking no empty intersection of genes
length( intersect( rownames(himmer_sp$DW), rownames(miller_mtx)))

## merging miller and himmer in a list
merged_himmer_miller <- c('Miller_dataset'=miller_mtx, himmer_sp)
saveRDS(merged_himmer_miller, '../data/merged_himmer_miller.rds')

############################################################################################################################
##                                                                                                                        ##
##                          Subsetting Satpathy data                                                                      ##
##                                                                                                                        ##
############################################################################################################################

## barcodes
summaryTCells <- readRDS('../data/satpathy2016/scATAC_TME_TCells_SummarizedExperiment.final.rds')

## peaks
peaks <- read.table('../data/satpathy2016/genes.tsv', header = TRUE)
peaks <- as.character(peaks$Feature)

## counts
counts <- readMM('../data/satpathy2016/matrix.mtx')
rownames(counts) <- as.character(peaks)
colnames(counts) <- summaryTCells$Group_Barcode


## cluster annotation
t_cell_clust <- summaryTCells$T_Cell_Cluster
cluster_labs <- paste0('Cluster', 1:19)        
cell_type <- c('1-Naive CD4 T', '2-Activated CD4', '3-Th1', '4-Memory CD4 T',
               '5-Th17', '6-Th 1', '7-Th 2', '8-Treg 1',
               '9-Treg 2', '10-Treg 3', '11-Treg 4', '12-Effector CD8 T',
               '13-Naive CD8 T', '14-Memory CD8 T', '15-Early TEx', '16-Intermediate TEx',
               '17-Terminal TEx', '18-Other T', '19-Other T')
cluster_cell_types <- plyr::mapvalues(t_cell_clust, 
                                      from = cluster_labs, 
                                      to = cell_type)
TEx <- c('15-Early TEx', 
         '16-Intermediate TEx',
         '17-Terminal TEx')
is_TEx <- sapply(cluster_cell_types, function(x) x %in% TEx)
TEx_names <- summaryTCells$Group_Barcode[is_TEx]

## Dropping non CD8 T Cells and duplicated barcodes
counts <- counts[, is_TEx]
saveRDS(counts, '../data/satpathy_cdtcells_only_counts.rds')




