#' Extract ordered matrix from prepared heatmaps list object
#'
#' @return A list of dataframes
#'
#' @param prepared_hm_list A \code{HeatmapList-class} object for which the 
#' layout has been created and which has been drawn up in advance. 
#' 
#' @examples \dontrun{
#' htl <- get_demo_heatmap_list()
#' htl <- prepare_heatmap_list(htl)
#' htl_df <- extract_count_matrix(htl)
#' }
#'
#' @importFrom magrittr %>%
#'
#' @export
extract_count_matrix <- function(prepared_hm_list){
  is.partitions <- FALSE
  
  row_order <- unlist(row_order(prepared_hm_list))
  hm_list <- prepared_hm_list@ht_list
  
  if("partitions" %in% names(hm_list)){
    is.partitions <- TRUE
    partitions <- prepared_hm_list@ht_list$partitions@matrix[row_order]
  }
  
  track_names <- names(hm_list)[names(hm_list) != "partitions"]
  
  hm_list_df <- list()
  
  for (track_name in track_names){
    hm_list_df[[track_name]] <- hm_list[[track_name]]@matrix %>%
      as.data.frame
    hm_list_df[[track_name]] <- hm_list_df[[track_name]][row_order,]
    
    if (is.partitions){
      hm_list_df[[track_name]]$partitions <- partitions
    }
  }
  
  hm_list_df
}