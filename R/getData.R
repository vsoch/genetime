#' getData
#'
#' Returns all RNA probes defined in developmental enrichment package
#' @keywords developmental gene lookup
#' @param justgenes Return only gene list (and not full data). Default is FALSE
#' @export
#' @examples
#' getData()

getData = function(justgenes=FALSE){
  
  data("interpolated_timeseries")
  
  if (justgenes){
    return(unique(interpolated$genes)) 
  } else {
   return(interpolated)    
  }
}