## Downloading Satpathy AT et al, 2019 data
barcodes_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129785&format=file&file=GSE129785%5FscATAC%2DTME%2DTCells%2Ecell%5Fbarcodes%2Etxt%2Egz'
matrix_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129785&format=file&file=GSE129785%5FscATAC%2DTME%2DTCells%2Emtx%2Egz'
peaks_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129785&format=file&file=GSE129785%5FscATAC%2DTME%2DTCells%2Epeaks%2Etxt%2Egz'

## dowloading data
if ( ! file.exists('../data/satpathy2016/barcodes.tsv')) {
        download.file(barcodes_url, destfile = '../data/hcv/satpathy2016/barcodes.tsv.gz')
        gunzip('../data/hcv/satpathy2016/barcodes.tsv.gz', 
               destname = '../data/hcv/satpathy2016/barcodes.tsv',
               remove = FALSE)  
}
if ( ! file.exists('../data/satpathy2016/matrix.mtx')){
        download.file(matrix_url, destfile = '../data/satpathy2016/matrix.mtx.gz')
        gunzip('../data/hcv/satpathy2016/matrix.mtx.gz', 
               destname = '../data/hcv/satpathy2016/matrix.mtx', 
               remove = FALSE)  
}
if ( ! file.exists('../data/satpathy2016/genes.tsv')){
        download.file(peaks_url, destfile = '../data/satpathy2016/genes.tsv.gz')
        gunzip('../data/hcv/satpathy2016/genes.tsv.gz', 
               destname = '../data/hcv/satpathy2016/genes.tsv',
               remove = FALSE) 
}
if ( ! file.exists('../data/satpathy2016/scATAC_TME_TCells_SummarizedExperiment.final.rds')){
        summarized_url <- 'https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/scATAC_TME_TCells_SummarizedExperiment.final.rds' 
        download.file(summarized_url, destfile = '../data/satpathy2016/scATAC_TME_TCells_SummarizedExperiment.final.rds')
}