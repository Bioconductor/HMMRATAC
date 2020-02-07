
library(MotifDb)
library(TFBSTools)
library(JASPAR2018)
library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

opts <- list()
opts['species'] <- '9606'
pfm <- getMatrixSet(JASPAR2018, opts)
pwm <- toPWM(pfm)

bs <- BSgenome.Hsapiens.UCSC.hg19
chrs <- lapply(bs, function(i) {
    bs[[paste0('chr', i)]]
})

siteseqs <- lapply(seq_along(chrs), function(i) {
    seqs <- searchSeq(pwm, chrs[[i]], seqname = paste0('chr', i), min.score="80%", strand = "*")
    as(seqs, "GRanges")
})

pro <- promoters(TXDb.Hsapiens.UCSC.hg19.knownGene)
peaks <- import.bed('gappedPeak')

motif_to_peak_mappings <- lapply(siteseqs, function(site) {
    subsetByOverlaps(site, peaks)
})

site_to_pro_counts <- lapply(motif_to_peak_mappings, function(mot) {
    table(countOverlaps(mot, pro) > 0)
})
