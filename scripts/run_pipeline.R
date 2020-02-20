## Run pipeline 
## Must be run in the curry cluster
library(liger)
#library(Signac)
library(Seurat)
#library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v75)
set.seed(1234)

## gene expression + atac integration
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

######################################################################################
##     										    ##
##           Integration of Himmer and Miller data                                  ##
##								                    ##
######################################################################################

#cat('Processing liger object\n')
#merged_himmer_miller <- readRDS('../data/merged_himmer_miller.rds')
#liger <- createLiger(merged_himmer_miller)

#cat('Running liger\n')
#run_liger(
#        liger = liger,
#        k = 10,
#        liger_path = '../data/liger_himmer_miller.rds',
#        plot_path = '../figures/liger_himmer_miller.pdf',  
#        interspecies = TRUE
#)

#liger <- readRDS('../data/liger_himmer_miller.rds')
## Saving liger results
#save_liger(
#        ligerex = liger,
#        dir_path = '../data/liger_himmer_miller_tsne.rds'
#)

#######################################################################
## Integration of Mike to Satpathy

merge <- readRDS('data/merged_himmer_miller.rds')
merge <- merge[names(merge)!='Miller_dataset']
satpathy <- readRDS('data/score_activity_matrix_cd8tcells.rds')
merge$'satpathy' <- satpathy
sapply(merge, class)
names(merge)
sapply(merge, function(x) head(colnames(x)))
sapply(merge, function(x) head(rownames(x)))

liger <- createLiger(merge)


cat('Running liger\n')
run_liger(
        liger = liger,
        k = 20,
        liger_path = 'data/liger_maike_satpathy.rds',
        plot_path = 'figures/liger_maike_satpathy.pdf',  
        interspecies = TRUE
)

#liger <- readRDS('../data/liger_himmer_miller.rds')
## Saving liger results
save_liger(
        ligerex = liger,
        dir_path = '../data/liger_maike_satpathy_tsne.rds'
)

