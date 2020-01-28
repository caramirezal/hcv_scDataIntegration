## Cistopic 
library(Matrix)
library(cisTopic)



cistopic_preprocessing <- function(counts) {
        cat('Preprocessing data... \n Selects number of peaks bigger than 5 and less than 10 percent of the number of cells... \n')
        peak.sum=rowSums(counts)
        max.ratio=0.01
        min.peak=5
        peak.selected=rownames(counts)[peak.sum<ceiling(ncol(counts)*max.ratio) & peak.sum>min.peak]
}
