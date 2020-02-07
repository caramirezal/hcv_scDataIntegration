## Liger implementation
## liger test
library(liger)
library(dplyr)

## read Miller rna data 
rna <- readRDS('~/sc/singleCellDeconvolution/data/cell_lines_rna.RDS')

## Read Satpathy atac data
atac <- readRDS('~/sc/singleCellDeconvolution/data/cell_lines_activity_matrix.RDS')
colnames(atac) <- paste0(colnames(atac), '.atac')

## preprocessing data
liger <- createLiger(list(rna=rna, atac=atac))
liger <- normalize(liger)
liger <- selectGenes(liger, combine = 'union')
liger <- scaleNotCenter(liger)


## running liger
liger <- optimizeALS(liger, k = 4) 
liger <- quantileAlignSNF(liger) #SNF clustering and quantile alignment

## saving results
saveRDS(liger, 'data/cell_lines_liger.RDS')


ligerex = runTSNE(ligerex)
# for larger datasets, may want to use UMAP instead
ligerex = runUMAP(ligerex)

pdf("liger_results_lcmv.pdf")
plotByDatasetAndCluster(ligerex) #Can also pass in different set of cluster labels to plot
plotFeature(ligerex, "nUMI")
plotWordClouds(ligerex)
plotGeneLoadings(ligerex)
dev.off()