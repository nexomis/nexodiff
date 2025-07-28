#' Helper function for ExprData$filter_and_set_selected_ids
#'
#' @param current_selected_ids private$selected_ids (see ExprData)
#' @param annotation private$annotation (see ExprData)
#' @param main_etag see ExprData$filter_and_set_selected_ids
#' @param values see ExprData$filter_and_set_selected_ids
#' @param filtered_var see ExprData$filter_and_set_selected_ids
#' @param filter_type see ExprData$filter_and_set_selected_ids
#' @return Vector of filtered expression tag IDs
filter_etags <- function(
  current_selected_ids, annotation, main_etag,
  values, filtered_var = "type", filter_type = "keep"
) {

  if (! filtered_var %in% c("tax_id", "tax_name", "type")) {
    logging::logerror("Wrong value for filtered_var argument")
    stop()
  }

  filter_dict <- annotation$generate_translate_dict(
    main_etag, filtered_var)

  if (filter_type == "keep") {
    keep_ids <- names(filter_dict)[filter_dict %in% values]
  } else if (filter_type == "excl") {
    keep_ids <- names(filter_dict)[!filter_dict %in% values]
  } else {
    logging::logerror("Wrong value for filter_type argument")
    stop()
  }

  intersect(keep_ids, current_selected_ids)
}
