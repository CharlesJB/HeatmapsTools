#' Get a bedGraph file for demo
#'
#' @return The path to a valid bdg file
#'
#' @examples
#' bdg <- get_demo_bdg()
#'
#' @export
get_demo_bdg <- function() {
    system.file("extdata/demo.bdg", package="HeatmapsTools", mustWork = TRUE)
}

#' Get an indexed bedGraph file for demo
#'
#' @return The path to a valid indexed bdg file
#'
#' @examples
#' bdg <- get_demo_bdg_idx()
#'
#' @export
get_demo_bdg_idx <- function() {
    system.file("extdata/demo.bdg.gz", package="HeatmapsTools", mustWork = TRUE)
}

#' Get a bed file for demo
#'
#' @return The path to a valid bed file
#'
#' @examples
#' bed <- get_demo_bed_files()
#'
#' @export
get_demo_bed_files <- function() {
    c(system.file("extdata/demo.bed", package="HeatmapsTools", mustWork = TRUE),
      system.file("extdata/demo2.bed", package="HeatmapsTools", mustWork = TRUE))

}

#' Return a list of heatmaps for demo
#'
#' @return A \code{list} of \code{Heatmap} objects
#'
#' @param partitions Should the heatmaps be splitted? \code{TRUE} or
#' \code{FALSE}. Default: \code{FALSE}.
#'
#' @examples
#' bdg <- get_demo_heatmap_list()
#'
#' @import purrr
#' @importFrom magrittr %>%
#'
#' @export
get_demo_heatmap_list <- function(partitions = FALSE) {
    cov <- import_bedgraph(get_demo_bdg())
    peaks <- map(get_demo_bed_files(), rtracklayer::import)

    if (!partitions) {
        heatmap1 <- produce_heatmap(cov = cov,
                                    peaks = peaks[[1]],
                                    name = "demo1",
                                    force_seqlevels = TRUE)
        heatmap2 <- produce_heatmap(cov = cov,
                                    peaks = peaks[[2]],
                                    name = "demo2",
                                    force_seqlevels = TRUE)
    } else {
        partitions <- c(rep(1, 10), rep(2, 10)) %>%
		          as.character
        heatmap1 <- produce_heatmap(cov = cov,
                                    peaks = peaks[[1]],
                                    name = "demo1",
                                    partitions = partitions,
                                    force_seqlevels = TRUE)
        heatmap2 <- produce_heatmap(cov = cov,
                                    peaks = peaks[[2]],
                                    name = "demo2",
                                    partitions = partitions,
                                    force_seqlevels = TRUE)

    }
    list(heatmap1 = heatmap1, heatmap2 = heatmap2)
}
