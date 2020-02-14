
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

## Read in pvalues for TFs
## pvals <- read.delim('pvalues_file.txt')

## Find which Genes have promoters with enriched TF's and open chromatin regions
chromosome_locs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
promoters <- getPromoterSeq(chromosome_locs, Hsapiens, upstream=2000, downstream=2000)

## Functional analysis to find (KEGG) pathways that are enriched
