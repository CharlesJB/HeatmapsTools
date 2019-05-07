#' Import bedgraph files
#'
#' @param filename Paths to a valid bedGraph file.
#' @param filter_negative_coverages Convert negative coverage values to 0?
#'     (Default: TRUE)
#' @param keep_standard_chromosomes Remove alternative chromosomes? 
#'     (Default: TRUE)
#' @param genome The name of the genome (i.e.: "hg38", "mm10", etc...)
#'
#' @return A list of GRanges (one element per file).
#'
#' @examples
#' filenames <- get_demo_bdg()
#' cov <- import_bedgraph(filenames)
#'
#' @import rtracklayer
#' @import purrr
#' @importFrom stringr str_detect
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom GenomeInfoDb seqlevels<-
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqinfo
#'
#' @export
import_bedgraph <- function(filename,
                            filter_negative_coverages = TRUE,
                            keep_standard_chromosomes = TRUE,
                            genome = NULL) {
    stopifnot(file.exists(filename))
    stopifnot(stringr::str_detect(filename, "bdg$|bdg\\.gz$"))
    stopifnot(is.logical(keep_standard_chromosomes))
    stopifnot(is.logical(filter_negative_coverages))

    genome <- get_si(genome, keep_standard_chromosomes)
    gr <- import(filename, format = "bedGraph")

    if (filter_negative_coverages) {
        gr$score[gr$score < 0] <- 0
    }
    if (keep_standard_chromosomes) {
        gr <- GenomeInfoDb::keepStandardChromosomes(gr,
                                pruning.mode = "coarse")
    }
    if (!is.null(genome)) {
        seqlevels(gr) <- seqlevels(genome)
        seqinfo(gr) <- genome
    }
    gr
}

get_si <- function(genome, keep_standard_chromosomes) {
    if (!is.null(genome)) {
        stopifnot(is.character(genome))
        stopifnot(length(genome) == 1)
        stopifnot(any(names(si) == genome))
        genome <- si[[genome]]
        if (keep_standard_chromosomes) {
            genome <- GenomeInfoDb::keepStandardChromosomes(genome)
        }
    }
    genome
}
