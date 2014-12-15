#' probeLookup
#'
#' Returns overlap of user gene list with RNA associated genes defined in developmental enrichment package
#' @param genelist A list of gene symbol identifiers
#' @keywords developmental RNA gene lookup
#' @export
#' @examples
#' probeLookup()

probeLookup = function(genelist){
 genes = getData(justgenes=TRUE) 
 overlap = intersect(genes,genelist)
 cat("Found overlap of",length(overlap),"genes associated with RNA probes in developmental enrichment set of",length(genes),"\n")
 return(intersect(genes,genelist))
}