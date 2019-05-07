#' Combine the heatmaps objects
#'
#' This function prepare the ready for drawing \code{HeatmapList-class} object.
#' 
#' @return This function invisibly returns a \code{HeatmapList-class} object
#' for which the layout has been created.
#'
#' @param heatmap_list The \code{list} of \code{Heatmap} objects produced with
#' the \code{produce_heatmap} function.
#' @param partitions If you want to split the heatmap based on different group,
#' you must specify the groups (or partitions). The \code{partition} must be a
#' \code{vector} of the same length as the \code{peaks} param used when
#' creating the \code{Heatmap} objects with the \code{produce_heatmap}
#' function. All heatmaps must be produced using the same \code{peaks} object.
#' If the \code{vector} is empty (i.e.: \code{character(0)}, the heatmap won't
#' be splitted. Default: \code{character(0)}.
#' 
#' @examples \dontrun{
#' htl <- get_demo_heatmap_list()
#' htl <- prepare_heatmap_list(htl)
#' }
#'
#' @import purrr
#' @import EnrichedHeatmap
#' @importFrom ggthemes colorblind_pal
#'
#' @export
prepare_heatmap_list <- function(heatmap_list, partitions = character(0),
                                 seed = 99841) {
    stopifnot(is(heatmap_list, "list"))
    stopifnot(all(map_lgl(heatmap_list, is, "Heatmap")))
    stopifnot(is.character(partitions))
    if (length(partitions) > 0) {
        stopifnot(length(partitions) == nrow(heatmap_list[[1]]@matrix))
        stopifnot(length(unique(partitions)) > 1)
    }
    stopifnot(is.numeric(seed))

    if (length(partitions) > 0) {
        g <- unique(partitions)
        set.seed(seed)
        col_blind <- ggthemes::colorblind_pal()(length(g))
        heatmaps_partitions <- ComplexHeatmap::Heatmap(partitions,
                                       col = structure(col_blind),
                                       name = "",
                                       show_row_names = FALSE,
                                       width = grid::unit(3, "mm"))
        heatmap_list <- c(heatmap_list, heatmaps_partitions)
    }
    purrr::reduce(heatmap_list, ComplexHeatmap::`+.AdditiveUnit`)
}
