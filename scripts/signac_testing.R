## Signac testing
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
set.seed(1234)

## Raw data
if ( ! file.exists('data/signac_testing_data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5') ) {
      cat('Downloading Raw data\n')
      url <- 'http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5'
      download.file(url, 'data/signac_testing_data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')
}
## Metadata
if ( ! file.exists('data/signac_testing_data/atac_v1_pbmc_10k_singlecell.csv') ) {
      cat('Downloading Metadata\n')
      url <- 'http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv'
      download.file(url, 'data/signac_testing_data/atac_v1_pbmc_10k_singlecell.csv')
}
## Fragments
if ( ! file.exists('data/signac_testing_data/atac_v1_pbmc_10k_fragments.tsv.gz') ) {
      cat('Downloading Fragments\n')
      url <- 'http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz'
      download.file(url, 'data/signac_testing_data/atac_v1_pbmc_10k_fragments.tsv.gz')
}
## Fragments index
if ( ! file.exists('data/signac_testing_data/atac_v1_pbmc_10k_fragments.tsv.gz.tbi') ) {
      cat('Downloading Fragments index\n')
      url <- 'http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi'
      download.file(url, 'data/signac_testing_data/atac_v1_pbmc_10k_fragments.tsv.gz.tbi')
}

cat('Reading peaks counts from  h5\n')
counts <- Read10X_h5(filename = "data/signac_testing_data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")

cat('Reading metadata\n')
metadata <- read.csv(
  file = "data/signac_testing_data/atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

cat('Creating Seurat object\n')
pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

#cat('Setting fragments to Seurat object\n')
#fragment.path <- 'data/signac_testing_data/atac_v1_pbmc_10k_fragments.tsv.gz'
#pbmc <- SetFragments(
#  object = pbmc,
#  file = fragment.path
#)

#cat('Calculation of nucleosome signals\n')
#pbmc <- NucleosomeSignal(object = pbmc)
#pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
#pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

#head(pbmc$nucleosome_signal) 
#pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 10, 'NS > 10', 'NS < 10')

# create granges object with TSS positions
#gene.ranges <- genes(EnsDb.Hsapiens.v75)
#gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]

#tss.ranges <- GRanges(
#  seqnames = seqnames(gene.ranges),
#  ranges = IRanges(start = start(gene.ranges), width = 2),
#  strand = strand(gene.ranges)
#)

#seqlevelsStyle(tss.ranges) <- 'UCSC'
#tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

#cat('to save time use the first 2000 TSSs')
#pbmc <- TSSEnrichment(object = pbmc, tss.positions = tss.ranges[1:2000])
#pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')

#pbmc <- subset(pbmc, subset = peak_region_fragments > 1000 & peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05 & nucleosome_signal < 10 & TSS.enrichment > 2)
#pbmc

#cat('plotting results\n')
#pdf('figures/signac_testing.pdf')
#VlnPlot(
#  object = pbmc,
#  features = c('pct_reads_in_peaks', 'blacklist_ratio', 'nucleosome_signal', 'peak_region_fragments'),
#  pt.size = 0.1,
#  ncol = 4) + NoLegend()
#FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
#TSSPlot(pbmc, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()
#dev.off()


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(
  object = pbmc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 1:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 1:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE)

pdf('figures/signac_testing.pdf')
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()

