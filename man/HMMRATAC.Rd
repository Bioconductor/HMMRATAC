% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMMRATAC.R
\name{HMMRATAC}
\alias{HMMRATAC}
\title{Perform HMMRATAC}
\usage{
MMRATAC <- function(bam,
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
}
\arguments{
\item{bam}{Sorted BAM file containing the ATAC-seq reads.}

\item{index}{Index file for the sorted BAM File.}

\item{genome}{Two column, tab delimited file containing genome size
stats.}

\item{means}{numeric vector of initial mean values for the fragment
distribution. Default = c(50,200,400,600).}

\item{stddev}{numeric vector of initial standard deviation values
for fragment distribution. Default = c(20,20,20,20).}

\item{fragem}{logical(1) of whether to perform EM training on the
fragment distribution. Default = TRUE.}

\item{minmapq}{numeric(1) of ,inimum mapping quality of reads to
keep.  Default = 30.}

\item{upper}{numeric(1) of Upper limit on fold change range for
choosing training sites. Default = 20.}

\item{lower}{numeric(1) of lower limit on fold change range for
choosing training sites. Default = 10.}

\item{zscore}{numeric(1) of Zscored read depth to mask during
Viterbi decoding. Default = 100.}

\item{zscoreprescan}{numeric(1) of Minimum zscored read depth to be
included in Viterbi decoding. Default = 0.}

\item{output}{character(1) base name (including directory path) of
output files. Default: file.path(tempfile(), basename(bam)).}

\item{blacklist}{character(1) of name of bed file of blacklisted
regions to exclude.}

\item{peaks}{logical(1) of whether to report peaks in bed format.
Default = TRUE.}

\item{kmeans}{numeric(1) of number of states in model. Default =
3. If not k=3, recommend NOT calling peaks, use bedgraph.}

\item{training}{character(1) of name of BED file of training
regions to use for training model, instead of foldchange
settings.}

\item{bedgraph}{logical(1) Whether to report whole genome bedgraph
of all state anntations. Default = FALSE.}

\item{minlen}{numeric(1) of minimum length of open region to call
peak.  Note: peaks must be set. Default = 200.}

\item{score}{character(1) either "max", "ave", "med", "fc",
"zscore", "all".  What type of score system to use for
peaks. Can be used for ranking peaks. Default = max.}

\item{bgscore}{logical(1) of whether to add the HMMR score to each
state annotation in bedgraph. Note: this adds considerable
time.  Default = FALSE.}

\item{trim}{numeric(1) of how many signals from the end to trim off
(ie starting with tri signal then di etc). This may be useful
if your data doesn't contain many large fragments. Default = 0.}

\item{window}{numeric(1) of size of the bins to split the genome
into for Viterbi decoding. To save memory, the genome is split
into <int> long bins and viterbi decoding occurs across each
bin. Default = 25000000.  Note: For machines with limited
memory, it is recomended to reduce the size of the bins.}

\item{model}{character(1) of the name of a binary model file
(generated from previous HMMR run) to use instead of creating
new one.}

\item{modelonly}{logical(1) of whether or not to stop the program
after generating model. Default = false.}

\item{overwrite}{logical(1) overwrite existing `output` files?}

\item{verbose}{logical(1) report progress?}
}
\value{
matrix A matrix(?)
}
\description{
Perform HMMRATAC
}
\examples{
example_data <- system.file(package="HMMRATACData", "data")
example_files <- dir(example_data, recursive = TRUE, full = TRUE)
basename(example_files)

output_directory <- tempfile()
dir.create(output_directory)
output <- file.path(
    output_directory,
    sub(".bam$", "", basename(example_files)[[1]])
)
output

outputs <- HMMRATAC(
   bam = example_files[[1]],
   index = example_files[[2]],
   genome = example_files[[3]],
   output = output,
   window = 2500000
)
outputs

}
