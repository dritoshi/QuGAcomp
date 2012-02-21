##
## set functions
##
setGeneric("intersectNum", function(object) { standardGeneric("intersectNum") })
setMethod("intersectNum", "QuGAcomp",
  function(object) {
    sum(object@data[[1]] == object@data[[2]])
  }    
)
setGeneric("setdiffNum", function(object) { standardGeneric("setdiffNum") })
setMethod("setdiffNum", "QuGAcomp",
  function(object) {
    right.num <- sum(object@data[[1]] < object@data[[2]])
    left.num  <- sum(object@data[[1]] > object@data[[2]])
  
    return(list(left=left.num, right=right.num))
  }    
)
setGeneric("eachsetNum", function(object) { standardGeneric("eachsetNum") })
setMethod("eachsetNum", "QuGAcomp",
  function(object) {
    right.num <- length(object@data[[2]])
    left.num  <- length(object@data[[1]])
    return(list(left=left.num, right=right.num))
  }    
)
setGeneric("unionNum", function(object) { standardGeneric("unionNum") })
setMethod("unionNum", "QuGAcomp",
  function(object) {
    sum(unlist(eachsetNum(object))) - intersectNum(object)
  }    
)

#' Caluculate cross table
#'
#' This function caluculate cross table (contingency table)
#'
#' @usage crossTable(quga)
#' @param quga QuGAcomp object
#' @export
#' @importMethodsFrom IRanges table
#' @examples
#' data(quga)
#'
#' quga.mat <- crossTable(quga)
setGeneric("crossTable", function(object) { standardGeneric("crossTable") })
setMethod("crossTable", "QuGAcomp",
  function(object) {
    # match.vector <- IRanges::table( object@data[[1]] + object@data[[2]] )
    match.vector <- table( object@data[[1]] + object@data[[2]] )
    diff.vector  <- as.vector(unlist(setdiffNum(object)))
    mat <- matrix(
      c(match.vector[1], diff.vector[1], diff.vector[2], match.vector[3]),
      nrow=2
    )
    return(mat)
  }
)
