#' Caluculate Mutual information
#'
#' This function caluculate Mutual information
#'
#' @usage mutualInformation(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.mi <- mutualInformation(quga)
setGeneric("mutualInformation", function(object) { standardGeneric("mutualInformation") })
setMethod("mutualInformation", "QuGAcomp",
  function(object) {
    intersectoin.num <- intersectNum(object)
    setdiff.num      <- setdiffNum(object)
    eachset.num      <- eachsetNum(object)
    union.num        <- unionNum(object)
    return( - log2(intersectoin.num / prod(unlist(eachset.num)) * union.num) )
  }    
)

#' Caluculate dice index
#'
#' This function caluculate dice index
#'
#' @usage dice(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.dice <- dice(quga)
setGeneric("dice", function(object) { standardGeneric("dice") })
setMethod("dice", "QuGAcomp",
  function(object) {
    intersectoin.num <- intersectNum(object)
    setdiff.num      <- setdiffNum(object)
    eachset.num      <- eachsetNum(object)    
    union.num        <- unionNum(object)
    return( 2 * intersectoin.num / sum(unlist(eachset.num)) )
  }    
)

#' Caluculate Jaccard index
#'
#' This function caluculate Jaccard index
#'
#' @usage jaccard(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.jaccard <- jaccard(quga)
setGeneric("jaccard", function(object) { standardGeneric("jaccard") })
setMethod("jaccard", "QuGAcomp",
  function(object) {
    intersectoin.num <- intersectNum(object)
    setdiff.num      <- setdiffNum(object)
    eachset.num      <- eachsetNum(object)        
    union.num        <- unionNum(object)
    return( intersectoin.num / union.num )
  }    
)

#' Caluculate Simpson index
#'
#' This function caluculate Simpson index
#'
#' @usage simpson(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.simpson <- simpson(quga)
setGeneric("simpson", function(object) { standardGeneric("simpson") })
setMethod("simpson", "QuGAcomp",
  function(object) {
    intersectoin.num <- intersectNum(object)
    setdiff.num      <- setdiffNum(object)
    eachset.num      <- eachsetNum(object)    
    union.num        <- unionNum(object)
    return( intersectoin.num / min(unlist(eachset.num)) )
  }    
)

#' Caluculate Cosine coefficent
#'
#' This function caluculate Cosine coefficent
#'
#' @usage cosineCoef(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.cosine <- cosineCoef(quga)
setGeneric("cosineCoef", function(object) { standardGeneric("cosineCoef") })
setMethod("cosineCoef", "QuGAcomp",
  function(object) {
    intersectoin.num <- intersectNum(object)
    setdiff.num      <- setdiffNum(object)
    eachset.num      <- eachsetNum(object)    
    union.num        <- unionNum(object)
    return( intersectoin.num / sqrt(prod(unlist(eachset.num))) )
  }    
)

#' Caluculate Pearson product-moment correlation coefficient
#'
#' This function caluculate Pearson product-moment correlation coefficient
#'
#' @usage pearsonCoef(quga.ctable)
#' @param quga.ctable QuGAcomp CrossTable object
#' @export
#' @examples
#' data(quga)
#'
#' quga.pearson <- pearsonCoef(quga)
setGeneric("pearsonCoef", function(object) { standardGeneric("pearsonCoef") })
setMethod("pearsonCoef", "QuGAcomp",
  function(object) {
    cor(object@data[[1]], object@data[[2]])
  }
)

##
## Coefficient of association (contingency)
##

#' Caluculate chi squre of cross table
#'
#' This function caluculate chi squre of cross table (contingency table)
#' 
#' @usage chisq(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.chisq <- chisq(quga)
setGeneric("chisq", function(object) { standardGeneric("chisq") })
setMethod("chisq", "QuGAcomp",
  function(object) {
    Oij <- crossTable(object)
    Eij <- outer(rowSums(Oij), colSums(Oij)) / sum(Oij)
    sum( (Oij-Eij)^2 / Eij )
    #chisq.test(Oij)$p.value
  }    
)

#' Caluculate Phi coefficient of cross table
#'
#' This function caluculate Phi coefficient of cross table (contingency table).
#' (= mean square contingency coefficient)
#' (= Pearson correlation coefficient for two binary variables)
#'
#' @usage phiCoef(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.phi <- phiCoef(quga)
setGeneric("phiCoef", function(object) { standardGeneric("phiCoef") })
setMethod("phiCoef", "QuGAcomp",
  function(object) {
    Oij <- crossTable(object)    
    sqrt( chisq(object) / sum(Oij) )    
  }    
)

#' Caluculate Contingency coefficient C of cross table
#'
#' This function caluculate Contingency coefficient C of cross table (contingency table).
#'
#' @usage contingencyCcoef(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.contingencyC <- contingencyCcoef(quga)
setGeneric("contingencyCcoef", function(object) { standardGeneric("contingencyCcoef") })
setMethod("contingencyCcoef", "QuGAcomp",
  function(object) {
    Oij     <- crossTable(object)
    mychisq <- chisq(object)
    sqrt( mychisq / (sum(Oij) + mychisq) )
  }
)

#' Caluculate cramer's coefficient V of cross table
#'
#' This function caluculate cramer's coefficient V of cross table (contingency table).
#'
#' @usage cramerCoef(quga)
#' @param quga QuGAcomp object
#' @export
#' @examples
#' data(quga)
#'
#' quga.cramer <- cramerCoef(quga)
setGeneric("cramerCoef", function(object) { standardGeneric("cramerCoef") })
setMethod("cramerCoef", "QuGAcomp",
  function(object) {
    Oij <- crossTable(object)    
    phiCoef(object) / sqrt( min(dim(Oij)) -1)
  }
)
