#' Update a Criterion-Value after Swap based on a design from a RCB design
#'
#' \code{NewValue.rcbd} calculates an updated variance-covariance matrix and its trace or log of determinant after each swap of pairs of treatments from a RCB design
#'
#' @param matdf an xperimental design
#' @param criteria either ``A" or ``D"
#' @param Rinv a matrix calculated from the design
#' @param C11 a matrix calculated from the design
#' @param K a matrix calculated from matdf
#' @return trace or log(determinant) value
#'
#' @examples
#' trt = length(1:4);criteria="A"
#' blocks = 2; rb=2;cb=2;Tr=2;Tc=4;rhox=0.9;rhoy=0.9;h2=0.1;s20=0
#'
#' des1 = rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#' ans1 = VarCov.rcbd.fixedAll(matdf=des1,Tr,Tc, criteria="A")
#'
#' newmat <- SwapPair(des1)
#' which(des1[,"Treatments"] != newmat[,"Treatments"])
#' DesLayout(des1, trt, cb, rb, blocks)
#' DesLayout(matdf=newmat, trt, cb, rb, blocks)
#'
#' mats2 <- NewValue.rcbd.fixedAll(matdf = des1, criteria="A", Rinv=ans1$Rinv,
#'          K=ans1$K, C11=ans1$C11)
#' mats2
#'
#' @export
#'
#' @seealso  \code{\link{VarCov.rcbd.fixedAll}}

NewValue.rcbd.fixedAll <- function(matdf, criteria, Rinv, C11, K) {
  W <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatments"]) - 1)
  colnames(W) <- NULL
  W <- as.matrix(W)
  W = Matrix::drop0(W)

  M = solve(Matrix::crossprod(W,Rinv) %*% W - Matrix::crossprod(W,K) %*% W)
  M <- as(M, "sparseMatrix")
  M = round(M,7)

  if (criteria == "A") {
    return(traceI = sum(Matrix::diag(M)))
  }
  if (criteria == "D") {
    deTm = Matrix::det(M)
    return(doptimI = log(deTm))
  }
}
