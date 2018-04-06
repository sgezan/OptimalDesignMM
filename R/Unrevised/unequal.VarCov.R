#' Generate the traces or log of determinants and other related matrices
#'
#' \code{unequal.VarCov} calculates traces or log of determinants from information matrix
#'
#'
#' @param matdf an experimental layout
#' @param rhox spatial correlation along the rows
#' @param rhoy spatial correlation along the columns
#' @param VarG narrow-sense heritability
#' @param nugget nugget error
#' @param Tr total number of rows
#' @param Tc total number of columns
#' @param criteria either "A" or "D"
#' @param Amat a numerator relationship matrix. If FALSE, no pedigree is needed
#' @param sigBl variance of the blocks effect. If FALSE, it is calculated internally
#' @param regular a logical statement, if FALSE, a different algorithm of AR1 is used
#'
#' @return either a trace value or log of determinant from calculated based on Henderson mixed models solutions and other matrices such as an inverse of G, inverse of R and K.
#'
#' @references
#' Mramba, Lazarus. K. and Gezan, Salvador. A. (2016), Improving Unequally Replicated,
#' Incomplete Block and Augmented Experimental Designs with Spatially
#' and Genetically Correlated Observations, Submitted to the Journal of
#' Theoretical and Applied Genetics
#'
#' @examples
#' trt = length(1:30);criteria="A"
#' blocks = 3; rb=5;cb=6;Tr=5;Tc=18;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#' min.u = sample(1:3,trt,replace=TRUE)
#' max.u = sample(3:5,trt,replace=TRUE)
#' des1 = unequal.RBD(trt,blocks,max.u,min.u,rb,cb,Tr,Tc,plot=TRUE)
#'
#' data(ped30hs)
#' Amat <- GenA(male=ped30hs[,"male"], female = ped30hs[,"female"])
#' Amat <- as.matrix(Amat[-c(1:5), -c(1:5)])
#'
#' ans1 <- unequal.VarCov(des1$matdf,Tr,Tc, rhox=0.9,rhoy=0.9,
#' VarG=0.1,nugget=0, criteria="A",Amat=TRUE)
#'
#' attributes(ans1)
#' ans1$traceI
#' ans1$ef
#'
#' @export
#' @seealso \code{\link{DesLayout}} for a layout with treatments in string format, \code{\link{rcbd}}
#'
#'
unequal.VarCov <- function(matdf,Tr,Tc,rhox=0,rhoy=0,VarG=0.3,nugget=0,criteria="A",Amat=FALSE,sigBl=FALSE,regular=TRUE)
{
  if(regular ==TRUE & Tr*Tc != nrow(matdf)) stop("check Tr by Tc dimensions")
  X <- matrix(1,nrow = nrow(matdf))
  # determine number of blocks
  bb <- length(unique(matdf[,"Reps"]))
  if(is.numeric(sigBl))
  {
    Binv <- (1/sigBl)*Matrix::Diagonal(bb)
  }else{
    sigBl <- 0.2*(1 - VarG)
    Binv <- (1/sigBl)*Matrix::Diagonal(bb)
    Binv <- as(Binv, "sparseMatrix")
  }

  s2e <- (1 - nugget) * (1 - VarG - sigBl)
  stopifnot(s2e > 0)

  m = length(unique(matdf[,"Treatments"]))
  if(is.matrix(Amat)) {
    Gg <- VarG * as.matrix(Amat)
    Gg <- round(solve(Gg),7)
    Gg <- as(Gg, "sparseMatrix")
  }else{
    Gg <- (1/VarG) * Matrix::Diagonal(m)
    Gg <- as(Gg, "sparseMatrix")
  }
  Ginv <- Matrix::bdiag(Binv,Gg)

  Zg <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatments"]) - 1)
  Zb <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Reps"]) - 1)
  Z <- Matrix::cBind(Zb,Zg)

  if(regular==FALSE){
    R <- Matrix::Diagonal(nrow(matdf))
    for(i in 1:(nrow(matdf)-1)) {
      x1  <- matdf[,"Col"][i]
      y1  <- matdf[,"Row"][i]
      for (j in (i+1):nrow(matdf)){
        x2 <- matdf[,"Col"][j]
        y2  <- matdf[,"Row"][j]
        R[i,j]<-(rhox^abs(x2 -x1))*(rhoy^abs(y2 -y1))
      }
    }
    R = as.matrix(round(s2e*R,7))
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
    R <- as(R, "sparseMatrix")
    Rinv <- round(chol2inv(chol(R)),7)
    Rinv <- as(Rinv, "sparseMatrix")
  }

  if(regular==TRUE){
    sigx <- Matrix::Diagonal(Tc)
    sigx <- rhox^abs(row(sigx) - col(sigx))
    sigy <- Matrix::Diagonal(Tr)
    sigy <- rhoy^abs(row(sigy) - col(sigy))
    R <- round(s2e * kronecker(sigy, sigx),7)
    R <- as(R, "sparseMatrix")
    Rinv <- round(chol2inv(chol(R)),7)
    Rinv <- as(Rinv, "sparseMatrix")
  }

  C11 <- t(X) %*% Rinv %*% X
  C11inv <- 1/C11
  K <- round(Rinv %*% X %*% C11inv %*% t(X) %*% Rinv ,7)
  K <- as(K, "sparseMatrix")
  Z <- as.matrix(Z)
  temp0 <- t(Z) %*% Rinv %*% Z + Ginv - t(Z) %*% K %*% Z
  C22 <- solve(temp0)
  C22 <- round(C22[-(1:bb), -(1:bb)],7)
  C22 <- as(C22, "sparseMatrix")

  r <- max(matdf[,"Reps"])
  ef <- (2*(s2e+nugget)/r)/mean(lower.tri(C22))

  if (criteria == "A") {
    return(c(traceI = sum(Matrix::diag(C22)), Ginv = Ginv, Rinv = Rinv, K = K,ef=ef))
  }
  if (criteria == "D") {
    deTm = Matrix::det(C22)
    return(c(doptimI = log(deTm), Ginv = Ginv, Rinv = Rinv, K=K,ef=ef))
  }
}
