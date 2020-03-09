## Run pipeline 
## Must be run in the curry cluster
library(liger)
#library(Signac)
library(Seurat)
#library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v75)
set.seed(333)
library(SeuratWrappers)

## gene expression + atac integration
run_liger <- function(liger,
                      k=10,
                      liger_path,
                      plot_path,
                      interspecies=NULL,
                      lambda=5){
        cat('liger normalization and scaling\n')
        liger <- normalize(liger)
        liger <- selectGenes(liger, combine = 'union', capitalize=interspecies)
        liger <- scaleNotCenter(liger)
        
        cat('Running liger\n')
        liger <- optimizeALS(liger, k = k, lambda = lambda) 
        liger <- quantileAlignSNF(liger) 
        
        cat('Run TSNE\n')
        liger = runTSNE(liger)
        
        cat('Plotting visualizations\n')
        pdf(plot_path)
        plotByDatasetAndCluster(liger) #Can also pass in different set of cluster labels to plot
        plotFeature(liger, "nUMI")
        #plotWordClouds(liger)
        #plotGeneLoadings(liger)
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
## Integration of Maike to Satpathy

#merge <- readRDS('data/merged_himmer_miller.rds')
#merge <- merge[names(merge)!='Miller_dataset']
#satpathy <- readRDS('data/score_activity_matrix_TEx.rds')
#cat('Satpathy dim\n')
#dim(satpathy)
#merge$'satpathy' <- satpathy
#sapply(merge, class)
#names(merge)
#sapply(merge, function(x) head(colnames(x)))
#sapply(merge, function(x) head(rownames(x)))

#liger <- createLiger(merge)


#cat('Running liger\n')
#run_liger(
#        liger = liger,
#        k = 5, lambda = 0.1,
#        liger_path = 'data/liger_maike_satpathy_lambda=0.1.rds',
#        plot_path = 'figures/liger_maike_satpathy_lambda=0.1.pdf',  
#        interspecies = TRUE
#)


#liger <- readRDS('data/liger_maike_satpathy_lambda=0.1.rds')
## Saving liger results
#save_liger(
#        ligerex = liger,
#        dir_path = 'data/liger_maike_satpathy_tsne_lambda=0.1.rds'
#)


#calcAlignment(liger)
#calcAgreement(liger)
#calAlignment(liger)

#cat('Suggest K\n')
#pdf('figures/suggest_k.pdf')
#suggestK(liger, k.test=seq(5,10,20), plot.log2=T)
#dev.off()
#suggestLambda(liger, k=5)


##############################################################################
##								            ##
##             Alignment of Hofmann and Satpathy using Seurat               ##
##                                                                          ##
##############################################################################
scores <- readRDS('../data/score_activity_matrix_TEx.rds')
sam_ann <- readRDS('../data/score_activity_matrix_annotations_TEx.rds')

satpathy_seu <- CreateSeuratObject(
        counts = scores, 
        project = 'HCV', 
        assay = 'SCORES', 
        min.cells = 1, 
        min.features = 1
)
#satpathy_seu <- NormalizeData(satpathy_seu)
#satpathy_seu <- FindVariableFeatures(satpathy_seu)

hofmann_seu <- readRDS('../data/integrated_miller_maike_sub.rds')

seu_list <- list(hofmann_seu, satpathy_seu)
cat('Normalizing data')
for (i in 1:length(seu_list)){
    seu_list[[i]] <- NormalizeData(seu_list[[i]])  
    seu_list[[i]] <- FindVariableFeatures(seu_list[[i]])
}
#merge_seu <- merge(seu_list[[1]], seu_list[[2]])
#merge_seu <- ScaleData(merge_seu, split.by = "orig.ident", do.center = FALSE)
#merge_seu <- RunOptimizeALS(merge_seu, k = 5, lambda = 5, split.by = "orig.ident")
#merge_seu <- RunQuantileAlignSNF(merge_seu, split.by = "orig.ident")
#merge_seu <- RunUMAP(merge_seu, dims = 1:ncol(merge_seu[["iNMF"]]), reduction = "iNMF")

#saveRDS(merge_seu, '../data/hofmann_satpathy_liger_seurat.rds')
