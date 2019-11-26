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
#' @param genome Two column, tab delimited file containing genome size stats.
#' @param means numeric vector of initial mean values for the fragment
#'      distribution. Default = c(50,200,400,600).
#' @param stddev numeric vector of initial standard deviation values for
#'      fragment distribution. Default = c(20,20,20,20).
#' @param fragem logical(1) of whether to perform EM training on the fragment
#'      distribution. Default = TRUE.
#' @param minmapq numeric(1) of ,inimum mapping quality of reads to keep.
#'      Default = 30.
#' @param upper numeric(1) of Upper limit on fold change range for choosing
#'      training sites. Default = 20.
#' @param lower numeric(1) of lower limit on fold change range for choosing
#'      training sites. Default = 10.
#' @param zscore numeric(1) of Zscored read depth to mask during Viterbi
#'      decoding. Default = 100.
#' @param zscoreprescan numeric(1) of Minimum zscored read depth to be included
#'      in Viterbi decoding. Default = 0.
#' @param output character(1) of directory of output files. Default = tempdir().
#' @param blacklist character(1) of name of bed file of blacklisted regions to
#'      exclude.
#' @param peaks logical(1) of whether to report peaks in bed format.
#'      Default = TRUE.
#' @param kmeans numeric(1) of number of states in model. Default = 3. If not
#'      k=3, recommend NOT calling peaks, use bedgraph.
#' @param training character(1) of name of BED file of training regions to use
#'      for training model, instead of foldchange settings.
#' @param bedgraph logical(1) Whether to report whole genome bedgraph of all
#'      state anntations. Default = FALSE.
#' @param minlen numeric(1) of minimum length of open region to call peak.
#'      Note: peaks must be set. Default = 200.
#' @param score character(1) either "max", "ave", "med", "fc", "zscore", "all".
#'      What type of score system to use for peaks. Can be used for ranking
#'      peaks. Default = max.
#' @param bgscore logical(1) of whether to add the HMMR score to each state
#'      annotation in bedgraph. Note: this adds considerable time.
#'      Default = FALSE.
#' @param trim numeric(1) of how many signals from the end to trim off
#'      (ie starting with tri signal then di etc). This may be useful if your
#'      data doesn't contain many large fragments. Default = 0.
#' @param window numeric(1) of size of the bins to split the genome into for
#'      Viterbi decoding. To save memory, the genome is split into <int> long
#'      bins and viterbi decoding occurs across each bin. Default = 25000000.
#'      Note: For machines with limited memory, it is recomended to reduce the
#'      size of the bins.
#' @param model character(1) of the name of a binary model file
#'      (generated from previous HMMR run) to use instead of creating new one.
#' @param modelonly logical(1) of whether or not to stop the program after
#'      generating model. Default = false.
#'
#' @return matrix A matrix(?)
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
                     output_dir = tempdir(),
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
                     modelonly = FALSE)
{

args <- c()

if(missing(bam)) {
    stop("The mandatory argument 'bam' is missing.")
}
args <- c(args, '-b', bam)

if(missing(index)) {
    stop("The mandatory argument 'index' is missing.")
}
args <- c(args, '-i', index)

if(missing(genome)) {
    stop("The mandatory argument 'genome' is missing.")
}
args <- c(args, '-g', genome)

if(!missing(means)) {
    ## FIXME: parse means
    args <- c(args, '-m', means)
}

if(!missing(stddev)) {
    ## FIXME: parse stddev
    args <- c(args, '-s', stddev)
}

if(!missing(fragem)) {
    if (fragem) fragem <- 'True'
    else fragem <- 'False'
    args <- c(args, '-f', fragem)
}

if(!missing(minmapq)) {
    args <- c(args, '-q', minmapq)
}

if(!missing(upper)) {
    args <- c(args, '-u', upper)
}

if(!missing(lower)) {
    args <- c(args, '-l', lower)
}

if(!missing(zscore)) {
    args <- c(args, '-z', zscore)
}

if(!missing(output_dir)) {
    args <- c(args, '-o', output_dir)
}

if(!missing(blacklist)) {
    args <- c(args, '-e', blacklist)
}

if(!missing(peaks)) {
    args <- c(args, '-p', peaks)
}

if(!missing(kmeans)) {
    args <- c(args, '-k', kmeans)
}

if(!missing(training)) {
    args <- c(args, '-t', training)
}

if(!missing(bedgraph)) {
    args <- c(args, '--bedgraph', bedgraph)
}

if(!missing(minlen)) {
    args <- c(args, '--minlen', minlen)
}

if(!missing(score)) {
    args <- c(args, '--score', score)
}

if(!missing(bgscore)) {
    args <- c(args, '--bgscore', bgscore)
}

if(!missing(trim)) {
    args <- c(args, '--trim', trim)
}

if(!missing(window)) {
    args <- c(args, '--window', window)
}

if(!missing(model)) {
    args <- c(args, '--model', model)
}

if(!missing(zscore)) {
    args <- c(args, '--modelonly', modelonly)
}



# .jcall("HMMR_ATAC.Main_HMMR_Driver", "V", "main", .jarray(list(), "java/lang/String"))
hmmr <- rJava::new(J("HMMR_ATAC.Main_HMMR_Driver"))

browser()

#args <- as.list(args)

#args <- lapply(args, function(x) {.jnew}

hmmr$main(args)
#hmmr$main(.jarray(args, "java/lang/String"))
}
