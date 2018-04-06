#' Update a Criterion-Value after Swap based on a design from a RCB design
#'
#' \code{NewValue.rcbd} calculates an updated variance-covariance matrix and its trace or log of determinant after each swap of pairs of treatments from a RCB design
#'
#' @param matdf an xperimental design
#' @param criteria either ``A" or ``D"
#' @param Rinv a matrix calculated from the original matdf
#' @param Ginv a matrix calculated from the original matdf
#' @param K a matrix calculated from the original matdf
#'
#' @return trace or log(determinant) value
#'
#' @examples
#' trt = length(1:30);criteria="A"
#' blocks = 6; rb=5;cb=6;Tr=15;Tc=12;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#'
#' matdf = rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#'
#' # If pedigree is available
#' # data(ped30fs)
#' # Amat <- GenA(male=ped30fs[,"male"],female=ped30fs[,"female"])
#' # Amat <- as.matrix(Amat[-c(1:5),-c(1:5)])
#' # G <- VarG*as.matrix(Amat)
#' # Ginv <- round(chol2inv(chol(as.matrix(G))),7)
#'
#' # else,  Treatments are independent
#' m = length(unique(matdf[,"Treatments"]))
#' Ginv <- round((1/VarG) * Matrix::Diagonal(m),7)
#' # Ginv <- as.matrix(Ginv)
#' # Ginv <- as(Ginv, "sparseMatrix")
#'
#' Rinv <- Rmatrix(matdf,VarG=0.1)
#' Rinv[1:5,1:5]
#'
#' res <- VarCov.rcbd(matdf,Ginv,Rinv,criteria="A")
#' NewValue.rcbd(matdf, Rinv, Ginv, K=as.matrix(res$K), criteria="A")
#'
#' res2 <- VarCov.rcbd(matdf,Ginv,Rinv,criteria="A",K = as.matrix(res$K), Update=TRUE)
#' res2
#'
#' @export
#'
#' @seealso  \code{\link{Optimize.rcbd}} to improve the designs

NewValue.rcbd <- function(matdf, Rinv, Ginv, K,criteria="A") {
  Z <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatments"]) - 1)
  Z <- as.matrix(Z)
  temp0 <- t(Z) %*% Rinv %*% Z + Ginv - t(Z) %*% K %*% Z
  C22 <- solve(temp0)
  if (criteria == "A") {
    traceI <- sum(Matrix::diag(C22))
    return(traceI)
  }
  if (criteria == "D") {
    doptimI <- log(Matrix::det(C22))
    return(doptimI)
  }
}
