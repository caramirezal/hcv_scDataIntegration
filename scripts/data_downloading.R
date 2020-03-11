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


## Downloading Miller Data
barcodes_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122712&format=file&file=GSE122712%5Fbarcodes%2Etsv%2Egz'
genes_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122712&format=file&file=GSE122712%5Fgenes%2Etsv%2Egz'
counts_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122712&format=file&file=GSE122712%5Fmatrix%2Emtx%2Egz' 
if (! file.exists('../data/miller2019/barcodes.tsv')){
        download.file(barcodes_url, destfile = '../data/miller2019/GSE122712_barcodes.tsv.gz')
        gunzip('../data/miller2019/GSE122712_barcodes.tsv.gz',
               destname = '../data/miller2019/barcodes.tsv', 
               remove = FALSE)
}
if (! file.exists('../data/miller2019/genes.tsv')) {
        download.file(genes_url, destfile = '../data/miller2019/GSE122712_genes.tsv.gz')
        gunzip('../data/miller2019/GSE122712_genes.tsv.gz', 
               destname = '../data/miller2019/genes.tsv', 
               remove = FALSE)
}
if (! file.exists('../data/miller2019/matrix.mtx.gz')) {
        download.file(counts_url, destfile = '../data/miller2019/matrix.mtx.gz')
        gunzip('../data/miller2019/matrix.mtx.gz', 
               destname = '../data/miller2019/matrix.mtx', 
               remove = FALSE)
}

## Satpathy bam files
trace_url <- 'https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP192525'
download.file(trace_url, destfile = '../data/satpathy2016/srp_trace.html')
bam_html <- readLines('../data/satpathy2016/srp_trace.html')
bam_lines <- bam_html[grep('bam.bam', bam_html)] 
bam_names <- gsub('\\.bam.*', '.bam.1', bam_lines)
bam_urls <- gsub('.*https', 'https', bam_names)
bam_comm <- sapply(bam_urls, function(s) 
        paste('wget --tries=30 -P data/satpathy2016/bam/',
              s)
)

writeLines(bam_comm, '../scripts/download_satpathy_bam.sh')


