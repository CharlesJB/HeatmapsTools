#' Convert coverage to heatmap
#'
#' @return An heatmap object
#'
#' @param cov A coverage generated with the \code{import_bedgraphs} function
#' @param peaks A \code{GRanges} object representing the peaks to import
#' @param column_title The name of the heatmap that will be shown at the top
#' @param force_seqlevels If \code{TRUE}, remove regions that are not found
#'                in the coverage. Corresponds to \code{pruning.mode =
#'                "coarse"} in \code{?seqinfo}. Default: \code{FALSE}.
#'
#' @examples
#' bdg <- get_demo_bdg()
#' narrow_peak_file <- peakimport::get_demo_file("narrowPeak")
#' peaks <- peakimport::import_peaks(narrow_peak_file)
#' cov <- import_bedgraphs(bdg)
#' heatmap <- produce_heatmap(cov, peaks, "demo")
#'
#' @import EnrichedHeatmap
#' @import GenomeInfoDb
#' @import circlize
#' @importFrom magrittr %>%
#'
#' @export
produce_heatmap <- function(cov, peaks, column_title, force_seqlevels = FALSE) {
	stopifnot(is(cov, "GRanges"))
	stopifnot(length(cov) > 0)
	stopifnot(is(peaks, "GRanges"))
	stopifnot(length(peaks) > 0)

	stopifnot(is(force_seqlevels, "logical"))
	if (!force_seqlevels) {
		seqnames_cov <- GenomeInfoDb::seqnames(cov) %>% unique  %>% as.character
		seqnames_peaks <- GenomeInfoDb::seqnames(peaks) %>% unique %>% as.character
		stopifnot(all(seqnames_peaks %in% seqnames_cov))
	} else {
		GenomeInfoDb::seqlevels(peaks, pruning.mode = "coarse") <-
			GenomeInfoDb::seqlevels(cov)
	}

	m <- normalizeToMatrix(cov, peaks,
					  value_column = "score", # TODO: Add param
					  extend = 5000, # TODO: Add param
					  mean_mode = "w0", # TODO: Add param
					  w = 100) # TODO: Add param
	col <- quantile(m, c(0.0, 0.7, 0.95), na.rm = TRUE) %>%
		circlize::colorRamp2(c("white", "white", "red"))
	EnrichedHeatmap(m, col = col, column_title = column_title)
}
