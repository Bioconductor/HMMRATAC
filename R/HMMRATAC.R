#' Perform HMMRATAC
#'
#' @description Perform HMMRATAC
#'
#' @usage
#'  MMRATAC <- function(bam,
#'                      index,
#'                      genome,
#'                      means = c(50, 200, 400, 600),
#'                      stddev = c(20, 20, 20, 20),
#'                      fragem = TRUE,
#'                      minmapq = 30,
#'                      upper = 20,
#'                      lower = 10,
#'                      zscore = 100,
#'                      zscoreprescan = 0,
#'                      output_dir = tempdir(),
#'                      blacklist = tempfile(fileext=".bed"),
#'                      peaks = TRUE,
#'                      kmeans = 3,
#'                      training,
#'                      bedgraph = FALSE,
#'                      minlen = 200,
#'                      score = c("max", "ave", "med", "fc", "zscore", "all"),
#'                      bgscore = FALSE,
#'                      trim = 0,
#'                      window = 25000000,
#'                      model,
#'                      modelonly = FALSE)
#'
#'
#'
#' @param bam Sorted BAM file containing the ATAC-seq reads.
#' @param index Index file for the sorted BAM File.
#' @param genome Two column, tab delimited file containing genome size
#'     stats.
#' @param means numeric vector of initial mean values for the fragment
#'     distribution. Default = c(50,200,400,600).
#' @param stddev numeric vector of initial standard deviation values
#'     for fragment distribution. Default = c(20,20,20,20).
#' @param fragem logical(1) of whether to perform EM training on the
#'     fragment distribution. Default = TRUE.
#' @param minmapq numeric(1) of ,inimum mapping quality of reads to
#'     keep.  Default = 30.
#' @param upper numeric(1) of Upper limit on fold change range for
#'     choosing training sites. Default = 20.
#' @param lower numeric(1) of lower limit on fold change range for
#'     choosing training sites. Default = 10.
#' @param zscore numeric(1) of Zscored read depth to mask during
#'     Viterbi decoding. Default = 100.
#' @param zscoreprescan numeric(1) of Minimum zscored read depth to be
#'     included in Viterbi decoding. Default = 0.
#' @param output character(1) base name (including directory path) of
#'     output files. Default: file.path(tempfile(), basename(bam)).
#' @param blacklist character(1) of name of bed file of blacklisted
#'     regions to exclude.
#' @param peaks logical(1) of whether to report peaks in bed format.
#'     Default = TRUE.
#' @param kmeans numeric(1) of number of states in model. Default =
#'     3. If not k=3, recommend NOT calling peaks, use bedgraph.
#' @param training character(1) of name of BED file of training
#'     regions to use for training model, instead of foldchange
#'     settings.
#' @param bedgraph logical(1) Whether to report whole genome bedgraph
#'     of all state anntations. Default = FALSE.
#' @param minlen numeric(1) of minimum length of open region to call
#'     peak.  Note: peaks must be set. Default = 200.
#' @param score character(1) either "max", "ave", "med", "fc",
#'     "zscore", "all".  What type of score system to use for
#'     peaks. Can be used for ranking peaks. Default = max.
#' @param bgscore logical(1) of whether to add the HMMR score to each
#'     state annotation in bedgraph. Note: this adds considerable
#'     time.  Default = FALSE.
#' @param trim numeric(1) of how many signals from the end to trim off
#'     (ie starting with tri signal then di etc). This may be useful
#'     if your data doesn't contain many large fragments. Default = 0.
#' @param window numeric(1) of size of the bins to split the genome
#'     into for Viterbi decoding. To save memory, the genome is split
#'     into <int> long bins and viterbi decoding occurs across each
#'     bin. Default = 25000000.  Note: For machines with limited
#'     memory, it is recomended to reduce the size of the bins.
#' @param model character(1) of the name of a binary model file
#'     (generated from previous HMMR run) to use instead of creating
#'     new one.
#' @param modelonly logical(1) of whether or not to stop the program
#'     after generating model. Default = false.
#' @param overwrite logical(1) overwrite existing `output` files?
#' @param verbose logical(1) report progress?
#'
#' @return matrix A matrix(?)
#'
#' @examples
#' example_data <- system.file(package="HMMRATAC", "extdata")
#' example_files <- dir(example_data, recursive = TRUE, full = TRUE)
#' basename(example_files)
#'
#' output_directory <- tempfile()
#' dir.create(output_directory)
#' output <- file.path(
#'     output_directory,
#'     sub(".bam$", "", basename(example_files)[[1]])
#' )
#' output
#'
#' outputs <- HMMRATAC(
#'    bam = example_files[[1]],
#'    index = example_files[[2]],
#'    genome = example_files[[3]],
#'    output = output
#' )
#' outputs
#'
#' @import rJava
#' @importFrom Rsamtools BamFile
#' @importFrom rtracklayer BEDFile
#' @export
HMMRATAC <- function(bam,
                     index,
                     genome,
                     means = c(50, 200, 400, 600),
                     stddev = c(20, 20, 20, 20),
                     fragem = TRUE,
                     minmapq = 30,
                     upper = 20,
                     lower = 10,
                     zscore = 100,
                     zscoreprescan = 0,
                     output = sub(".bam$", "", basename(bam)),
                     blacklist = tempfile(fileext=".bed"),
                     peaks = TRUE,
                     kmeans = 3,
                     training,
                     bedgraph = FALSE,
                     minlen = 200,
                     score = c("max", "ave", "med", "fc", "zscore", "all"),
                     bgscore = FALSE,
                     trim = 0,
                     window = 25000000,
                     model,
                     modelonly = FALSE,
                     overwrite = FALSE,
                     verbose = FALSE)
{
    output_exts <- c(".model", "_peaks.gappedPeak", "_summits.bed")

    stopifnot(
        .is_scalar_character(bam),
        .is_scalar_character(index),
        .is_scalar_character(genome),
        .is_scalar_character(output),
        dir.exists(dirname(output)),
        .is_scalar_logical(overwrite),
        overwrite || !any(file.exists(paste0(output, output_exts))),
        .is_scalar_logical(verbose)
    )

    if(!missing(means)) {
        ## '-m 1,2,3,4'
        means <- paste(means, collapse = ",")
    }

    if(!missing(stddev)) {
        ## FIXME: parse stddev
        stddev <- paste(stddev, collapse = ",")
    }

    if(!missing(fragem)) {
        if (fragem) fragem <- 'True'
        else fragem <- 'False'
    }

    if (!missing(score))
        score <- match.arg(score)

    args <- c(
        '-b', bam,
        '-i', index,
        '-g', genome,
        '-o', output,
        if (!missing(means)) c('-m', means),
        if (!missing(stddev)) c('-s', stddev),
        if (!missing(fragem)) c('-f', fragem),
        if (!missing(minmapq)) c("-q", minmapq),
        if (!missing(upper)) c('-u', upper),
        if (!missing(lower)) c('-l', lower),
        if (!missing(zscore)) c('-z', zscore),
        if (!missing(blacklist)) c('-e', blacklist),
        if (!missing(peaks)) c('-p', peaks),
        if (!missing(kmeans)) c('-k', kmeans),
        if (!missing(training)) c('-t', training),
        if (!missing(bedgraph)) c('--bedgraph', bedgraph),
        if (!missing(minlen)) c('--minlen', minlen),
        if (!missing(score)) c('--score', score),
        if (!missing(bgscore)) c('--bgscore', bgscore),
        if (!missing(trim)) c('--trim', trim),
        if (!missing(window)) c('--window', window),
        if (!missing(model)) c('--model', model),
        if (!missing(modelonly)) c('--modelonly', modelonly)
    )

    if (verbose) {
        arglines <- paste(args, collapse = " ")
        arglines <- unlist(regmatches(
            arglines,
            gregexpr("--?[^-]+", arglines)
        ))
        message("args:\n  ", paste(arglines, collapse = "\n  "))
    }

    hmmr <- rJava::new(J("HMMR_ATAC.Main_HMMR_Driver"))
    status <- hmmr$main(args)

    outputs <- dir(dirname(output), full.names = TRUE)
    keep <- startsWith(outputs, paste0(output, "_"))
    outputs[keep]
}
