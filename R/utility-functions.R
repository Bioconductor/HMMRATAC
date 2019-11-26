#' @importFrom rtracklayer import.bed
gappedPeaks_to_bed <- function(output_dir)
{
    peak_file <- list.files(pattern = 'gappedPeak')
    import.bed(peak_file[1], extraCols = c(peakScore = 'character', V2 = 'character', V3 = 'character'))
}
