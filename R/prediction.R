#' Make predictions of cell types based on the single-cell RNA-seq digital expression
#' profiles using a supervised classifier, SuperCT
#' @param object CellESet object or features x barcodes expression matrix
#' @param species Either 'human' or 'mouse' for now
#' @param model Choose a supported model. i.e 'CellTypes', 'CellStates'
#' Please refer to https://github.com/weilin-genomics/SuperCT for more details.
#' @param use.cells Character vector specify cells to make prediction for
#' @param results.dir Specify directory to save downloaded required files for the prediction
#' @importFrom utils download.file read.csv untar
#' @return
#' Predicted cell identities saved in \code{object@meta.data[['pred_types']]}
#' or a data frame with cell identities
#' @export
#' @references
#' Xie Peng and Gao Mingxuan (2019) \emph{SuperCT: a supervised-learning framework for
#' enhanced characterization of single-cell transcriptomic profiles},
#' \url{https://doi.org/10.1093/nar/gkz116} \emph{Nucleic Acids Research}
#' @examples
#' \dontrun{
#' myces <- PredCellTypes(myces, species = 'human', model = '38CellTypes', results.dir = '.')
#' }
PredCellTypes <- function(object, species = 'human', model = '38CellTypes', use.cells = NULL, results.dir = '.'){
  classOfObj <- class(x = object)
  if(classOfObj == 'CellESet'){
    mat <- as.matrix(x = object@raw.data)
    if(is.null(x = use.cells)){
      use.cells <- object@use.cells
    } else{
      use.cells <- intersect(x = use.cells, y = colnames(x = object@raw.data))
    }
  } else if(classOfObj == 'matrix'){
    mat <- object
    if(is.null(x = use.cells)){
      use.cells <- colnames(x = object)
    } else{
      use.cells <- intersect(x = use.cells, y = colnames(x = object))
    }
  } else{
      stop('A CellESet or matrix object expected.')
  }
  if(length(x = use.cells) == 0)
    stop('No available cell names.')
  if(!species %in% c('human', 'mouse'))
    stop('Unrecognized species. Specify human or mouse.')
  mat <- mat[, use.cells, drop = FALSE]
  # prepare required files used in prediction
  dir.create(path = results.dir, showWarnings = FALSE)
  if(startsWith(x = results.dir, prefix = '~'))
    results.dir <- gsub(pattern = '~', replacement = Sys.getenv('HOME'), results.dir, fixed = TRUE)
  dirname <- paste0(results.dir, '/', model)
  filename <- paste0(dirname, '.tar.gz')
  if(!dir.exists(paths = dirname)){
    if(!file.exists(filename)){
      download.file(url = paste0('https://raw.githubusercontent.com/weilin-genomics/rSuperCT_models/blob/master/',
                                 model, '/', model, '.tar.gz'),
                    destfile = filename)
    } else{
      untar(tarfile = filename, exdir = results.dir)
    }
  }
  all_files <- list.files(path = dirname, full.names = TRUE, recursive = TRUE)
  genes_file <- all_files[grep(pattern = 'genes.csv', x = all_files, perl = TRUE)][1]
  if(!file.exists(genes_file))
    stop('genes file missing.\n')
  homogenes <- read.csv(file = genes_file, header = FALSE, stringsAsFactors = FALSE)

  types_file <- all_files[grep(pattern = 'types.csv', x = all_files, perl = TRUE)][1]
  if(!file.exists(types_file))
    stop('types file missing.\n')
  celltypes <- read.csv(file = types_file, header = TRUE, stringsAsFactors = FALSE)

  bias_files <- all_files[grep(pattern = 'b(\\d)+.csv', x = all_files, perl = TRUE)]
  weights_files <- all_files[grep(pattern = 'w(\\d)+.csv', x = all_files, perl = TRUE)]
  if(length(x = bias_files) != length(x = weights_files))
    stop(paste0("\nbias files: ", length(x = bias_files), '\n',
                "weights_files: ", length(x = weights_files), '\n',
                "Equal length expected.\n"))

  mat <- trans_dge(mat = mat, species = species, homogenes = homogenes)
  # use trained model to predict cell types
  for(wf in weights_files){
    suffix <- unlist(x = regmatches(x = wf, m = gregexpr(pattern = "(\\d)+.csv", text = wf, perl = TRUE)))
    bf <- grep(pattern = suffix, bias_files, value = T)
    w_tmp <- read.csv(file = wf, header = FALSE)
    b_tmp <- read.csv(file = bf, header = FALSE)

    w_tmp <- as.matrix(x = w_tmp)
    b_tmp <- as.matrix(x = b_tmp)

    mat <- t(x = w_tmp) %*% mat
    b_mat <- matrix(data = rep(x = b_tmp, ncol(x = mat)),
                    nrow = nrow(x = mat),
                    ncol = ncol(x = mat))
    mat <- mat + b_mat
    mat[mat < 0] <- 0 #ReLU transformation
  }
  # get identity of each cell
  pred_types <- celltypes[max.col(m = t(x = mat)), 'celltype']
  if(classOfObj == 'CellESet'){
    object@meta.data$pred_types <- pred_types
    return(object)
  } else{
    return(data.frame(cell = colnames(x = mat), pred_types = pred_types))
  }
}
