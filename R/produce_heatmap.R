#' Convert coverage to heatmap
#'
#' @return An heatmap object
#'
#' @param cov A coverage generated with the \code{import_bedgraph} function.
#' @param peaks A \code{GRanges} object representing the peaks to import. Peaks
#' will be resized to 5000 base pairs to make sure all regions are of the same
#' length.
#' @param name The name of the heatmap that will be shown at the top.
#' @param partitions If you want to split the heatmap based on different group,
#' you must specify the groups (or partitions). The \code{partition} must be a
#' \code{vector} of \code{character} of the same length as the \code{peaks}
#' param. If the \code{vector} is empty (i.e. \code{character(0)}, the heatmap
#' won't be splitted. Default: \code{character(0)}.
#' @param col A \code{circlize::colorRamp2} object. If \code{NULL}, default
#' values will be used.
#' @param col_order A \code{numeric} object corresponding to the order of the
#' colors in the right annotation. To use when the colors are not matchning
#' with the top_annotation. Default: \{NULL}.
#' @param force_seqlevels If \code{TRUE}, remove regions that are not found
#'                in the coverage. Corresponds to \code{pruning.mode =
#'                "coarse"} in \code{?seqinfo}. Default: \code{FALSE}.
#' @param seed Set a seed value for the color selection. Only used when
#' \code{partition} is not \code{NULL}. The value must be the same as the one
#' used for the \code{prepare_heatmap_list} function. Default: 99841.
#' @param top_anno_ylim Manually set the ylim value for the top_annotation.
#' If \code{NULL}, value is inferred using the data. Default: \code{NULL}.
#' @param show_row_names Show matrix rownames on left side.
#' If \code{FALSE}, row names are hidden Default: \code{FALSE}.
#' @param row_labels Manually set row label on the left side.
#' If \code{NULL}, no row label. Default: \code{NULL}.
#'
#' @examples
#' bdg <- get_demo_bdg()
#' bed_file <- get_demo_bed_files()[[1]]
#' peaks <- import(bed_file)
#' cov <- import_bedgraph(bdg)
#' heatmap <- produce_heatmap(cov, peaks, "demo")
#'
#' @import EnrichedHeatmap
#' @import GenomeInfoDb
#' @import circlize
#' @importFrom magrittr %>%
#'
#' @export


produce_heatmap <- function(cov, 
                            peaks, 
                            name, 
                            partitions = character(0), 
                            col = NULL, 
                            col_order = NULL, 
                            force_seqlevels = FALSE, 
                            seed = 99841, 
                            top_anno_ylim = NULL, 
                            show_row_names = FALSE, 
                            row_labels = NULL){
  stopifnot(is(cov, "GRanges"))
  stopifnot(length(cov) > 0)
  stopifnot(is(peaks, "GRanges"))
  stopifnot(length(peaks) > 0)
  stopifnot(is.character(name))
  stopifnot(nchar(name) > 0)
  stopifnot(is.character(partitions) | is.factor(partitions))
  if (length(partitions) > 0) {
    stopifnot(length(partitions) == length(peaks))
    stopifnot(length(unique(partitions)) > 1)
  }
  stopifnot(is(force_seqlevels, "logical"))
  stopifnot(is.numeric(seed))
  stopifnot(is.null(top_anno_ylim) | is.numeric(top_anno_ylim))
  stopifnot(is.logical(show_row_names))
  stopifnot(is.null(row_labels) | is.character(row_labels) | is.list(row_labels))
  
  if (!force_seqlevels) {
    seqnames_cov <- GenomeInfoDb::seqnames(cov) %>% 
      unique %>% 
      as.character
    seqnames_peaks <- GenomeInfoDb::seqnames(peaks) %>% 
      unique %>% 
      as.character
    stopifnot(all(seqnames_peaks %in% seqnames_cov))
  } else {
    GenomeInfoDb::seqlevels(peaks, pruning.mode = "coarse") <- 
      GenomeInfoDb::seqlevels(cov)
  }
  
  peaks <- GenomicRanges::resize(peaks, 1, fix = "center")
  
  m <- EnrichedHeatmap::normalizeToMatrix(cov, peaks,
                                          value_column = "score", # TODO: Add param
                                          extend = 5000, # TODO: Add param
                                          mean_mode = "w0", # TODO: Add param
                                          w = 100) # TODO: Add param
  
  if (is.null(col)) {
    col <- quantile(m, c(0.0, 0.7, 0.95), na.rm = TRUE) %>% 
      circlize::colorRamp2(c("white", "white", "red"))
  }
  
  if (length(partitions) > 0) {
    g <- unique(partitions)
    set.seed(seed)
    col_blind <- (ggthemes::colorblind_pal())(length(g))
    if (!is.null(col_order)) {
      col_blind <- col_blind[col_order]
    }
    
    axis_params <- list(facing = "left")
    ha <- ComplexHeatmap::HeatmapAnnotation(
      lines = EnrichedHeatmap::anno_enriched(
        gp = grid::gpar(col = col_blind),
        axis_param = list(facing = "inside"), 
        ylim = top_anno_ylim))
  } else {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      lines = EnrichedHeatmap::anno_enriched(
        axis_param = list(facing = "inside"),
        ylim = top_anno_ylim))
  }
  
  if (is.null(row_labels)) {
    EnrichedHeatmap::EnrichedHeatmap(m, col = col, name = name, 
                                     column_title = name, top_annotation = ha, 
                                     show_row_names = show_row_names, 
                                     row_names_side = "left")
  } else {
    EnrichedHeatmap::EnrichedHeatmap(m, col = col, name = name, 
                                     column_title = name, top_annotation = ha, 
                                     show_row_names = show_row_names, 
                                     row_names_side = "left",
                                     row_labels = row_labels[rownames(m)], 
                                     row_names_gp = gpar(fontsize = 7))
  }
}
