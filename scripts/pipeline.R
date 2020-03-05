## Run pipeline 
## Must be run in the curry cluster
library(liger)
library(Seurat)
library(Matrix)
library(tidyverse)
library(class)

## standard liger run
run_liger <- function(liger,
                      k=10,
                      liger_path,
                      plot_path,
                      interspecies=NULL){
        cat('liger normalization and scaling\n')
        liger <- normalize(liger)
        liger <- selectGenes(liger, combine = 'union', capitalize=interspecies)
        liger <- scaleNotCenter(liger)
        
        cat('Running liger\n')
        liger <- optimizeALS(liger, k = k) 
        liger <- quantileAlignSNF(liger) 
        
        cat('Run TSNE\n')
        liger = runTSNE(liger)
        
        cat('Plotting visualizations\n')
        pdf(plot_path)
        plotByDatasetAndCluster(liger) #Can also pass in different set of cluster labels to plot
        plotFeature(liger, "nUMI")
        plotWordClouds(liger)
        plotGeneLoadings(liger)
        dev.off()
        
        cat('Saving results\n')
        saveRDS(liger, liger_path)
}

## saving liger results for plotting
save_liger <- function(ligerex,
                       dir_path = '.') {
        cat('Creating list object\n')
        liger_res <- list(
                H=ligerex@H,
                cell_data=ligerex@cell.data,
                H_norm=ligerex@H.norm,
                W=ligerex@W,
                V=ligerex@V,
                tsne_coords=ligerex@tsne.coords,
                alignment_clusters=ligerex@alignment.clusters,
                clusters=ligerex@clusters)
        
        cat('Saving results\n')
        saveRDS(liger_res, dir_path)
}

#############################################################################################
##                                                                                         ##
##               Batch correction and integration of Miller along with                     ##
##                       T CD8 Exhausted cells TEx datasets                                ##
##                                                                                         ##
#############################################################################################
## This section must be run in cluster
cat('Processing liger object\n')
merged_himmer_miller <- readRDS('../data/merged_himmer_miller.rds')
liger <- createLiger(merged_himmer_miller)

cat('Running liger\n')
run_liger(
        liger = liger,
        k = 10,
        liger_path = '../data/liger_himmer_miller.rds',
        plot_path = '../figures/liger_himmer_miller.pdf',  
        interspecies = TRUE
)

liger <- readRDS('../data/liger_himmer_miller.rds')
cat('Saving liger results\n')
save_liger(
        ligerex = liger,
        dir_path = '../data/liger_himmer_miller_tsne.rds'
)

############################################################################################
##                                                                                        ##
##                 Label transfer from Miller to TEx dataset                              ##
##                                                                                        ##
############################################################################################

cat('Loading TEx splitted dataset\n')
TEx_sp <- readRDS('../data/himmer_splitted.rds')

## Creating seurat objects for every sample
TEx_sp_seu <- lapply(1:length(TEx_sp), 
       function(i) CreateSeuratObject(counts = TEx_sp[i][[1]], 
                                      project = 'hcv', 
                                      assay = 'RNA', 
                                      min.cells = 1, 
                                      min.features = 1) 
)
names(TEx_sp_seu) <- names(TEx_sp)

## Loading Miller data
miller <- readRDS('../data/miller_seu.rds')

## Loading Miller annotations
miller_ann <- readRDS('../data/miller_annotations.rds')
miller$'cell_type' <- miller_ann$cell_type

## merging datasets
merged_seu <- TEx_sp_seu
merged_seu$'miller' <- miller

########################################################################################################
## Batch correction using Seurat
for (i in 1:length(merged_seu)) {
        merged_seu[[i]] <- NormalizeData(merged_seu[[i]], verbose = FALSE)
        merged_seu[[i]] <- FindVariableFeatures(merged_seu[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
}


reference.list <- merged_seu[names(TEx_sp)]
anchors <- FindIntegrationAnchors(object.list = reference.list, 
                                  dims = 1:30, 
                                  k.filter = 30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = 30)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = 'pca', dims = 1:30)

## **Plotting integrated data** 
DimPlot(integrated, reduction = 'umap', group.by = 'orig.ident', pt.size = 1.5)


########################################################################################################
## Label transfer using reference based transfer of information
cat('Loading TEx splitted dataset\n')
TEx_sp <- readRDS('../data/himmer_splitted.rds')

## Creating seurat objects for every sample
TEx_sp_seu <- lapply(1:length(TEx_sp), 
                     function(i) CreateSeuratObject(counts = TEx_sp[i][[1]], 
                                                    project = 'hcv', 
                                                    assay = 'RNA', 
                                                    min.cells = 1, 
                                                    min.features = 1) 
)
names(TEx_sp_seu) <- names(TEx_sp)

## Loading Miller data
miller <- readRDS('../data/miller_seu.rds')

## Loading Miller annotations
miller_ann <- readRDS('../data/miller_annotations.rds')
miller$'cell_type' <- miller_ann$cell_type

## merging datasets
merged_seu <- TEx_sp_seu
merged_seu$'miller' <- miller

for (i in 1:length(merged_seu)) {
        merged_seu[[i]] <- NormalizeData(merged_seu[[i]], verbose = FALSE)
        merged_seu[[i]] <- FindVariableFeatures(merged_seu[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
        merged_seu[[i]] <- ScaleData(merged_seu[[i]], verbose = FALSE)
        merged_seu[[i]] <- RunPCA(merged_seu[[i]], npcs = 30, verbose = FALSE)
        merged_seu[[i]] <- RunUMAP(merged_seu[[i]], reduction = "pca", dims = 1:30)
        merged_seu[[i]] <- SCTransform(merged_seu[[i]], verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = merged_seu, nfeatures = 3000)
merged_seu <- PrepSCTIntegration(object.list = merged_seu, anchor.features = features)
reference_dataset <- which(names(merged_seu) == "miller")
anchors <- FindIntegrationAnchors(object.list = merged_seu, normalization.method = "SCT", 
                                       anchor.features = features, reference = reference_dataset)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 50)

integrated <- RunPCA(object = integrated, verbose = FALSE)
integrated <- RunUMAP(object = integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:15, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = 0.1, verbose = FALSE, weights = 100, node.sizes = 10)

integrated$cell_type <- sapply(integrated$cell_type,
                               function(x) ifelse(is.na(x), 'TEx samp', x) )

DimPlot(integrated, group.by = 'orig.ident', reduction = 'umap', pt.size = 1.5)
DimPlot(integrated, group.by = 'cell_type', reduction = 'umap', pt.size = 1.5)

## to be corrected
#predictions <- TransferData(anchorset = anchors, refdata = as.character(miller$cell_type), 
#                            dims = 1:15, weight.reduction = 'umap', k.weight = 50)

cell_type <- integrated$cell_type
preds[1:length(predictions$predicted.id)] <- predictions$predicted.id

#integrated$preds <- preds
integrated$'dataset' <- sapply(integrated$orig.ident, function(x) 
                               ifelse(x=='lcmv', 'Miller', 'Maike'))
#sizes <- rep(1, length(integrated$orig.ident))
#sizes[1:length(predictions$predicted.id)] <- 3

##**Cell type and seurat cluster coincides**
DimPlot(integrated, 
        group.by = 'cell_type', 
        reduction = 'umap', split.by = 'dataset', pt.size = 1.3) + ggtitle('By cell type')
DimPlot(integrated, 
        group.by = 'seurat_clusters', 
        reduction = 'umap', split.by = 'dataset', pt.size = 1.3) + ggtitle('By Seurat cluster')


integrated$'predicted' <- plyr::mapvalues(as.character(integrated$seurat_clusters), 
                                          from = c('0', '1', '2', '3'),
                                          to = c('Effector', 'Terminally Ex',
                                                 'Progenitor Ex', 'Proliferating'))
DimPlot(integrated, 
        group.by = 'predicted', 
        reduction = 'umap', split.by = 'dataset', pt.size = 1.3, label = TRUE, label.size = 8) + NoLegend()

saveRDS(integrated, '../data/integrated_miller_maike.rds')

## Subsetting integrated data to Maike imputed dataset  
integrated_sub <- subset(integrated, dataset != 'Miller')
saveRDS(integrated_sub, '../data/integrated_miller_maike_sub.rds')

##################################################################################################################
##                                                                                                              ##
##                    Imputing cell types from Satpathy to Miller                                               ##
##                                                                                                              ##
##################################################################################################################
## Here 2-D are used
satpathy_tsne <- subset(tsne_ann, dataset == 'Satpathy')
hofmann_tsne <- subset(tsne_ann, dataset == 'Hofmann')
preds <- knn(train = as.matrix(select(satpathy_tsne, tSNE1,tSNE2)), 
             test = as.matrix(select(hofmann_tsne, tSNE1,tSNE2)),
             cl = satpathy_tsne$cell_type, k = 200)

######################################################################################################################
##                                                                                                                  ##
##                   Finding differentially expressed genes in Imputed clusters                                     ##
##                                                                                                                  ##
######################################################################################################################
## Reading hofmann dataset

hofmann_seu <- readRDS('../data/integrated_miller_maike_sub.rds') 

ord <- match(colnames(hofmann_seu), hofmann_tsne$Barcodes)
hofmann_tsne_ord <- hofmann_tsne[ord,]

sum(colnames(hofmann_seu) == as.character(hofmann_tsne_ord$Barcodes ))
hofmann_seu$'seurat_clusters' <- preds

markers <- FindAllMarkers(hofmann_seu, min.pct = 0.25)

head(markers, n=5)





