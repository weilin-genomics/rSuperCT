#' Import expression profiles from other object as object of class CellESet
#' @param object An object of class \code{seurat}, \code{Seurat}, \code{matrix}(can be sparse) or \code{data.frame}
#' in which rows represent for features and columns for cells.
#' @importFrom methods new validObject
#' @import Matrix
#' @return
#' An object of class CellESet(Cell Expression Set) with raw data stored in slot \code{raw.data}.
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' myces <- ImportData(pbmc_small)
#' myces
#' }
ImportData <- function(object){
  classOfObj <- class(object)[1]
  if(classOfObj == 'matrix' | classOfObj == 'dgCMatrix' | classOfObj == 'dgTMatrix')
    data <- object
  if(classOfObj == 'data.frame')
    data <- as.matrix(x = object)
  if(classOfObj == 'seurat')
    data <- object@raw.data
  if(classOfObj == 'Seurat')
    data <- object@assays$RNA@counts
  if(is.null(x = rownames(x = data)))
    stop('Feature names not found.')
  data <- as(object = data, Class = 'dgCMatrix')
  if(is.null(x = colnames(x = data)))
    colnames(x = data) <- paste0('cell_', 1:ncol(x = data))
  meta.data <- data.frame(
    nGene = apply(X = data, MARGIN = 2, FUN = function(x) sum(x > 0)),
    nUMI = apply(X = data, MARGIN = 2, FUN = function(x) sum(x)),
    id = 'cell',
    check.rows = FALSE,
    row.names = colnames(x = data),
    stringsAsFactors = FALSE
  )
  ces <- new(Class = 'CellESet',
             raw.data = data,
             meta.data = meta.data,
             use.cells = colnames(x = data),
             version = '0.1.0')
  validObject(object = ces)
  return(ces)
}
