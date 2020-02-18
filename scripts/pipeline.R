## Run pipeline 
## Must be run in the curry cluster
library(liger)
library(Seurat)

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

###########################################################################################
## Reciprocal PCA
for (i in 1:length(merged_seu)) {
        merged_seu[[i]] <- NormalizeData(merged_seu[[i]], verbose = FALSE)
        merged_seu[[i]] <- FindVariableFeatures(merged_seu[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
        merged_seu[[i]] <- ScaleData(merged_seu[[i]], verbose = FALSE)
        merged_seu[[i]] <- RunPCA(merged_seu[[i]], npcs = 30, verbose = FALSE)
}


############################################################################################

## Standard preprocessing
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
DimPlot(integrated, reduction = 'umap', group.by = 'orig.ident')

query <- merged_seu[["miller"]]
query <- ScaleData(query, verbose = FALSE)
query <- RunPCA(query, npcs = 30, verbose = FALSE)
query <- RunUMAP(query, reduction = "pca", dims = 1:30)
anchors <- FindTransferAnchors(reference = query, query = integrated, 
                                        dims = 1:10, reduction = 'pca')
predictions <- TransferData(anchorset = anchors, refdata = query$orig.ident, 
                            dims = 1:15, weight.reduction = 'pca', k.weight = 20)
query <- AddMetaData(query, metadata = predictions)


########################################################################################################
## Integration based in reference

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
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 20)

integrated <- RunPCA(object = integrated, verbose = FALSE)
integrated <- RunUMAP(object = integrated, dims = 1:30)

integrated$cell_type <- sapply(integrated$cell_type,
                               function(x) ifelse(is.na(x), 'TEx samp', x) )

DimPlot(integrated, group.by = 'orig.ident', reduction = 'umap')
DimPlot(integrated, group.by = 'cell_type', reduction = 'umap')

predictions <- TransferData(anchorset = anchors, refdata = query$cell_type, 
                            dims = 1:15, weight.reduction = 'pca', k.weight = 20)
