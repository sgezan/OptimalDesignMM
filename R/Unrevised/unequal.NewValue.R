#' Update the trace or determinant after a design alteration
#'
#' \code{unequal.NewValue} calculates an updated variance-covariance matrix and its trace or log of determinant after a swapping pairs of treatment within or accross blocks or following a  replacement of some treatments as given by a list of contraints for the minimum or maximum number of allowable treatments.
#'
#' @param matdf an xperimental design
#' @param criteria either ``A" or ``D"
#' @param Rinv a matrix calculated from the original matdf
#' @param Ginv a matrix calculated from the original matdf
#' @param K a matrix calculated from the original matdf
#'
#' @return trace or log(determinant) value
#'
#' @references
#' Mramba, Lazarus. K. and Gezan, Salvador. A. (2016), Improving Unequally Replicated,
#' Incomplete Block and Augmented Experimental Designs with Spatially
#' and Genetically Correlated Observations, Submitted to the Journal of
#' Theoretical and Applied Genetics
#'
#'
#' @examples
#' trt = length(c(1:30));criteria="A"
#' blocks = 3; rb=5;cb=6;Tr=5;Tc=18;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#' min.u = sample(1:3,trt ,replace=TRUE)
#' max.u = sample(3:5,trt ,replace=TRUE)
#' des1 = unequal.RBD(trt,blocks,max.u,min.u,rb,cb,Tr,Tc,plot=TRUE)
#'
#' data(ped30hs)
#' Amat <- GenA(male=ped30hs[,"male"], female = ped30hs[,"female"])
#' Amat <- as.matrix(Amat[-c(1:5), -c(1:5)])
#'
#' ans1 <- unequal.VarCov(des1$matdf,Tr,Tc, rhox=0.6,rhoy=0.6,
#'            VarG=0.3,nugget=0, criteria="A",Amat=TRUE)
#'
#' attributes(ans1)
#' ans1$traceI
#' ans1$ef
#'
#' newmat <- SwapPair(des1$matdf)
#' unequal.NewValue(newmat, criteria="A", Rinv=ans1$Rinv,
#'  Ginv = ans1$Ginv, K=ans1$K)
#'
#' @export
#'
#' @seealso \code{\link{unequal.RBD}}

unequal.NewValue <- function(matdf, criteria="A", Rinv, Ginv, K) {
  Zg <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatments"]) - 1)
  Zb <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Reps"]) - 1)
  Z <- Matrix::cBind(Zb,Zg)
  Z <- as.matrix(Z)
  temp0 <- t(Z) %*% Rinv %*% Z + Ginv - t(Z) %*% K %*% Z
  C22 <- solve(temp0)
  bb <- length(unique(matdf[,"Reps"]))
  C22 <- C22[-(1:bb), -(1:bb)]
  C22 <- as(C22, "sparseMatrix")

  if (criteria == "A") {
    traceI <- sum(Matrix::diag(C22))
    return(traceI)
  }
  if (criteria == "D") {
    doptimI <- log(Matrix::det(C22))
    return(doptimI)
  }
}

