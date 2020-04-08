#' class definition for class CellESet
#' @slot raw.data Raw matrix for expression profile.
#' @slot meta.data Data frame for storing the meta information of cells.
#' @slot use.cells Cells used for prediction. Default is all cells.
#' @slot version package version
#' @importFrom methods setClass
CellESet <- setClass(Class = 'CellESet',
                     slots = c(raw.data = 'ANY',
                               meta.data = 'data.frame',
                               use.cells = 'vector',
                               version = 'ANY')
)
#' show method
#' @param object A CellESet object
#' @importFrom methods setMethod
#' @docType methods
setMethod(
  f = 'show',
  signature = 'CellESet',
  definition = function(object) {
    cat('An object of class ',
        class(object), ': ',
        nrow(object@raw.data),
        ' genes across ',
        ncol(object@raw.data),
        ' cells.\n',
        sep = '')
    invisible(NULL)
  }
)
