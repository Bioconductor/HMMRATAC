#' @importFrom rtracklayer import.bed
gappedPeaks_to_bed <- function(output_dir)
{
    peak_file <- list.files(pattern = 'gappedPeak')
    import.bed(peak_file[1], extraCols = c(peakScore = 'character', V2 = 'character', V3 = 'character'))
}

.is_scalar_character <-
    function(x, allow.na = FALSE, allow.zchar = FALSE)
{
    is.character(x) && length(x) == 1L && (allow.na || !is.na(x)) &&
        (allow.zchar || nzchar(x))
}

.is_scalar_logical <-
    function(x, allow.na = FALSE)
{
    is.logical(x) && length(x) == 1L && (allow.na || !is.na(x))
}
