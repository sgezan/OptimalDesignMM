#' Update a Criterion-Value after Swap based on a design from a RCB design with random blocks and fixed treatments.
#'
#' @param matdf an xperimental design
#' @param criteria either ``A" or ``D"
#' @param K a matrix calculated from the original matdf
#' @param Rinv is an inverse of the spatial correlation matrix
#'
#' @return trace or log(determinant) value
#'
#' @examples
#' # Example 1
#' trt = length(1:9);criteria="A"
#' blocks = 2; rb=3;cb=3;Tr=3;Tc=6;rhox=0.6;rhoy=0.6;h2=0.3;s20=0
#'
#' des1 = rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#' ans1 = VarCov.rcbd.blocksRandom(matdf=des1, Tr, Tc, criteria="A")
#' attributes(ans1)
#'
#' # Example 2
#' criteria="D"
#' ans2 = VarCov.rcbd.blocksRandom(matdf=des1, Tr, Tc, criteria="D")
#' attributes(ans2)
#'
#' # Try swapping a pair of treatments and recalculate a criterion value
#' newmat <- SwapPair(des1)
#' which(des1[,"Treatments"] != newmat[,"Treatments"])
#' DesLayout(des1, trt, cb, rb, blocks)
#' DesLayout(newmat, trt, cb, rb, blocks)
#'
#' K = as.matrix(ans1$K)
#' Rinv = as.matrix(ans1$Rinv)
#'
#' mats2 <- NewValue.rcbd.blocksRandom(newmat, criteria, K=K,Rinv = Rinv)
#' mats2
#'
#' @export
#'
#' @seealso  \code{\link{Optimize.rcbd}} to improve the designs
#'
NewValue.rcbd.blocksRandom <- function(matdf, criteria, K, Rinv) {
  if(nrow(matdf)==length(unique(matdf[,"Treatments"]))){
    X <- as.matrix(matdf[, "Treatments"])
    colnames(X)<-NULL
  }
  if(nrow(matdf) > length(unique(matdf[,"Treatments"]))){
    X <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatments"])-1)
    colnames(X)<-NULL
  }
  X <- as.matrix(X)
  Rinv = as.matrix(Rinv)

  C11 <- Matrix::crossprod(X, Rinv) %*% X

  temp0 <- C11 - t(X) %*% K %*% X
  M <- solve(temp0)
  M  <- as(M, "sparseMatrix")
  M = round(M,7)
  X = Matrix::drop0(X)

  if (criteria == "A") {
    traceI <- sum(Matrix::diag(M))
    return(traceI)
  }
  if (criteria == "D") {
    doptimI <- log(Matrix::det(M))
    return(doptimI)
  }
}
