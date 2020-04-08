#' prepare the expression profile and complete it if there's genes not found
#' @param mat features x barcodes matrix
#' @param homogenes homogenes
#' @param species Either 'human' or 'mouse' for now.
#' @export
#' @return a (0,1)-matrix
trans_dge <- function(mat, homogenes, species = 'human'){
  homogenes <- switch(EXPR = species, human = homogenes[,1], mouse = homogenes[,2])
  loc <- match(x = homogenes, table = rownames(x = mat))
  mat <- mat[loc,]
  mat[is.na(x=mat)] <- 0
  rownames(x = mat) <- homogenes
  mat[mat > 0] <- 1
  return(mat)
}

