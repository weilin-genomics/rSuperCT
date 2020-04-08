#' plot histogram of predicted cell types
#' @param object CellESet object.
#' @importFrom ggplot2 ggplot geom_bar geom_text ylab xlab theme element_text
#' theme_classic labs aes_string guides
#' @export
#' @return
#' A ggplot object
#' @examples
#' \dontrun{
#' table(myces@meta.data$pred_types)
#' plotHist(myces)
#' }
plotHist <- function(object){
  classOfObj <- class(x = object)
  if(classOfObj == 'CellESet'){
    pred_types <- object@meta.data$pred_types
  } else if(classOfObj == "data.frame"){
    pred_types <- object$pred_types
  } else{
    stop('A CellESet object expected.')
  }
  types_cnt <- as.data.frame(x = sort(x = table(pred_types), decreasing = TRUE))
  types_cnt$prop  <- paste0(round(100 * types_cnt$Freq/sum(types_cnt$Freq), 2), '%')
  g <- ggplot(data = types_cnt, mapping = aes_string(x = 'pred_types', y = 'Freq', fill = 'pred_types')) +
    geom_bar(stat = 'identity') +
    guides(fill = FALSE) +
    geom_text(mapping = aes_string(label = 'prop'), data = types_cnt, stat = 'identity', vjust = -0.5)  +
    theme_classic() +
    ylab(label = 'Frequency') +
    xlab(label = NULL) +
    labs(title = 'Distribution of predicted cell types') +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.7))
  return(g)
}
