#!/mnt/lustre/usr/global/centos7/R/3.6.1/bin/R

library(BiocParallel)
library(MotifDb)
library(TFBSTools)
library(JASPAR2018)
library("BSgenome.Hsapiens.UCSC.hg19")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(KEGGREST)
library(GenomicFeatures)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0)
    stop("Provide gappedPeak file as an argument")

## Step 1: Match JASPAR TFs to chromosomes
opts <- list()
opts['species'] <- '9606'
pfm <- getMatrixSet(JASPAR2018, opts)
pwm <- toPWM(pfm)
pwm <- pwm[1]

bs <- BSgenome.Hsapiens.UCSC.hg19
chrs <- bplapply(list(17), function(i) {
    bs[[paste0('chr', i)]]
})

peaks <- import.bed(paste0(args, "_peaks.gappedPeak")) #import.bed('/Users/da42327_ca/Monocytes/01_results/NA_peaks.gappedPeak')
pro <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)

siteseqs <- bplapply(seq_along(chrs), function(i) {
    seqs <- searchSeq(pwm, chrs[[i]], seqname = paste0('chr', i), min.score="80%", strand = "*")
    as(seqs, "GRanges")
})

## Step 2: Perform ORA to find signficant TFBS's
run_seq3 <- function(seqname, motif_id, pwms, bsgenome, peaks, promoters) {
    message("Motif ID: ", motif_id, " -- Seqname: ", seqname)
    ## Find binding sites on chromosome
    siteseq <- searchSeq(pwms[motif_id], bsgenome[[seqname]], seqname = seqname, min.score="80%", strand = "*")
    tfbs <- as(siteseq, "GRanges")
    promoters_chr <- promoters[seqnames(promoters) == seqname]
    peaks_chr <- peaks[seqnames(peaks) == seqname]

    ## Get promoters with TF
    idx <- promoters_chr %over% tfbs
    promoters_with_tf <- promoters_chr[idx]
    promoters_without_tf <- promoters_chr[!idx]

    ## Summary statistics for phyper
    npromoters_with_tf <- sum(countOverlaps(promoters_chr, tfbs) > 0)
    npromoters_without_tf <- length(promoters_chr) - npromoters_with_tf
    
    npromoters_with_tf_in_open_chromatin <- sum(countOverlaps(peaks_chr, promoters_with_tf) > 0)
    npromoters_without_tf_in_open_chromatin <- sum(countOverlaps(peaks_chr, promoters_without_tf) > 0)
    
    ## Return values
    c(q = npromoters_with_tf_in_open_chromatin,
      m = npromoters_with_tf,
      n = npromoters_without_tf,
      k = npromoters_with_tf_in_open_chromatin + npromoters_without_tf_in_open_chromatin)
}

perform_run_by_motif <- function(motif_id, pwm, bs, peaks, pro) {
    seqnames <- paste0('chr', 1:22)
    motif_result <- lapply(seqnames, run_seq3, motif_id = motif_id, pwms = pwm, bsgenome = bs, peaks = peaks, promoters = pro)
    red <- Reduce(`+`, motif_result)
    do.call(phyper, as.list(red))
}

motif_names <- names(pwm)
motif_symbols <- vapply(pwm, function(motif) {
    motif@tags$symbol
}, character(1))
total_results <- bplapply(motif_names, perform_run_by_motif, pwm = pwm, bs = bs, peaks = peaks, pro = pro)
total_results <- matrix(total_results)
rownames(total_results) <- motif_symbols
write.table(total_results, file = paste0(args, "_total_results.txt"), colnames = TRUE)

##    c(q = npromoters_with_tf_in_open_chromatin, ## drawn white balls = promoters with TF in open chromatin
##      m = npromoters_with_tf, ## total white = promoters with TF
##      n = npromoters_without_tf, ## total black = promoters without TF
##      k = npromoters_with_tf_in_open_chromatin + npromoters_without_tf_in_open_chromatin) ## all promoters within open chromatin
## Step 2: Match motifs to peaks mapping and motifs to promoters

