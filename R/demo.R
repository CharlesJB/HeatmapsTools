#' Get a bedGraph file for demo
#'
#' @return The path to a valid bdg file
#'
#' @examples
#' bdg <- get_demo_bdg()
#'
#' @export
get_demo_bdg <- function() {
    system.file("extdata/demo.bdg", package="heatmaps", mustWork = TRUE)
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
    system.file("extdata/demo.bdg.gz", package="heatmaps", mustWork = TRUE)
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
    c(system.file("extdata/demo.bed", package="heatmaps", mustWork = TRUE),
      system.file("extdata/demo2.bed", package="heatmaps", mustWork = TRUE))

}
