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
        ligerex = runTSNE(ligerex)
        
        cat('Plotting visualizations\n')
        pdf(plot_path)
        plotByDatasetAndCluster(ligerex) #Can also pass in different set of cluster labels to plot
        plotFeature(ligerex, "nUMI")
        plotWordClouds(ligerex)
        plotGeneLoadings(ligerex)
        dev.off()
}




#run_liger(
#        rna_path = '../data/rna.rds',
#        gam_path = '../data/gene_activity_matrix_cd8tcells.rds',
#        k = 10,
#        liger_path = '../data/liger_cd8tcells.rds',
#        plot_path = '../figures/liger_cd8tcells.pdf',  
#        interspecies = TRUE
#)


liger_in <- readRDS('../data/TEx_liger')



