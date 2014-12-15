#' getColors
#'
#' Returns 16 colors for brain regions in developmental enrichment package
#' @keywords brain region colors
#' @export
#' @examples
#' colors = getColors()

getColors = function(){
  # 16 unique (hand selected) colors for brain regions
  colors = c("cadetblue3","slateblue4","tomato","chartreuse2","cornflowerblue","chocolate1","darkred","firebrick1","forestgreen","darkturquoise","purple","gold2","gray35","deeppink2","navajowhite3","red")
  return(colors)
}

