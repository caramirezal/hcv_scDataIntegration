## Liger implementation
## liger test
library(liger)
library(dplyr)



## gene expression + atac integration
run_liger <- function(rna_path,
                      gam_path,
                      k=10,
                      liger_path,
                      plot_path,
                      interspecies=NULL){
        cat('Reading rna dataset\n')
        rna <- readRDS(rna_path)
        
        cat('Readin gene activity matrix from ATAC seq\n')
        atac <- readRDS(gam_path)
        colnames(atac) <- paste0(colnames(atac), '.atac')
        
        cat('Preprocessing data\n')
        liger <- createLiger(list(rna=rna, atac=atac))
        liger <- normalize(liger)
        liger <- selectGenes(liger, combine = 'union', capitalize=interspecies)
        liger <- scaleNotCenter(liger)
        
        cat('Running liger\n')
        liger <- optimizeALS(liger, k = k) 
        liger <- quantileAlignSNF(liger) 
        
        cat('Saving results\n')
        saveRDS(liger, liger_path)
        
        cat('Run TSNE\n')
        liger = runTSNE(liger)
        
        cat('Plotting visualizations\n')
        pdf(plot_path)
        plotByDatasetAndCluster(liger) #Can also pass in different set of cluster labels to plot
        plotFeature(liger, "nUMI")
        plotWordClouds(liger)
        plotGeneLoadings(liger)
        dev.off()
}




## setting parameters and paths
#rna_path = '../data/rna_counts_miller.rds'
#gam_path = '../data/gene_activity_matrix_satpathy_cd8.rds'
#k = 10
#liger_path = '../liger_cd8tcells.rds'
#plot_path = '../liger_cd8tcells.pdf'  
#interspecies = TRUE

#cat('reading rna counts')
#rna <- readRDS(rna_path)
        
#cat('Readin gene activity matrix from ATAC seq\n')
#atac <- readRDS(gam_path)
#colnames(atac) <- paste0(colnames(atac), '.atac')
        
#cat('Preprocessing data\n')
#liger <- createLiger(list(rna=rna, atac=atac))
#liger <- normalize(liger)
#liger <- selectGenes(liger, combine = 'union', capitalize=interspecies)
#liger <- scaleNotCenter(liger)
        
#cat('Running liger\n')
#liger <- optimizeALS(liger, k = k) 
#liger <- quantileAlignSNF(liger) 
        
        
#cat('Run TSNE\n')
#liger = runTSNE(liger)

#cat('Saving results to ', liger_path, '\n')
#saveRDS(liger, liger_path)
        
#cat('Plotting visualizations to ', plot_path, '\n')
#pdf(plot_path)
#plotByDatasetAndCluster(liger) #Can also pass in different set of cluster labels to plot
#plotFeature(liger, "nUMI")
#plotWordClouds(liger)
#plotGeneLoadings(liger)
#dev.off()

liger_save <- function(liger,
                       path_dest = '.'
             ){
      res <- list(
                  H_norm = liger@H.norm,
                  H = liger@H,
                  tsne = liger@tsne.coords,
                  W = liger@W,
                  V = liger@V,
                  cluster = liger@clusters
             )
       cat('Saving results \n')
       saveRDS(res, path_dest)
}


#saveRDS(res, '../data/liger_cd8tcells_tsne.rds')
#liger_in <- readRDS('../data/TEx_batch_correction.RDS')
