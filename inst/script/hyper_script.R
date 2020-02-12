
library(BiocParallel)
library(MotifDb)
library(TFBSTools)
library(JASPAR2018)
library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(KEGGREST)
library(GenomicFeatures)

## Step 1: Match JASPAR TFs to chromosomes
opts <- list()
opts['species'] <- '9606'
pfm <- getMatrixSet(JASPAR2018, opts)
pwm <- toPWM(pfm)

bs <- BSgenome.Hsapiens.UCSC.hg19
chrs <- bplapply(list(17), function(i) {
    bs[[paste0('chr', i)]]
})

siteseqs <- bplapply(seq_along(chrs), function(i) {
    seqs <- searchSeq(pwm, chrs[[i]], seqname = paste0('chr', i), min.score="80%", strand = "*")
    as(seqs, "GRanges")
})

## Step 2: Match motifs to peaks mapping and motifs to promoters

pro <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)
ipro <- ranges(pro)
peaks <- import.bed('/Users/da42327_ca/Monocytes/01_results/NA_peaks.gappedPeak')
ipeaks <- ranges(peaks)
imotifs <- ranges(siteseqs)

motif_to_peak_mappings <- bplapply(siteseqs, function(site) {
    subsetByOverlaps(site, peaks)
})

site_to_pro_counts <- bplapply(motif_to_peak_mappings, function(mot) {
    table(countOverlaps(mot, pro) > 0)
    
})

## Step 3: perform phyper and aquire pvalues for each TFBS

## Step 4: Download KEGG pathways and aquire gene transcripts

chromosome_locs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

## Find enriched peaks / TFBS gene overlaps in Pathways

## functional enrichment test
