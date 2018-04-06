#' Generate the traces or log of determinants and other related matrices
#'
#' \code{unequal.Augmented.VarCov} calculates relevant matrices and traces or determinants from a variance-covariance information matrix
#'
#'
#' @param matdf an experimental layout
#' @param rhox spatial correlation along the rows
#' @param rhoy spatial correlation along the columns
#' @param VarG treatments variability
#' @param nugget unstructured errors (measurement errors)
#' @param rb number of rows per block
#' @param cb number of columns per block
#' @param criteria either "A" or "D"
#' @param Amat a numerator relationship matrix. If FALSE, no pedigree is needed
#' @param sigBl variance of the blocks effect. If FALSE, it is calculated internally
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
#' CheckPlots = c(paste0(c("C"),1:2))
#' trt.list = c(paste0(c("g"),1:80))
#' Reps.Per.Block = 5
#' rb = 10
#' cb = 5
#' blocks = 2
#' matdf = rcbd.Augmented(blocks,trt.list,CheckPlots,
#'                       Reps.Per.Block, rb,cb, plot=TRUE)
#'
#' newblocks <- matrix(trt.list,ncol=blocks,byrow = FALSE)
#' chk= matrix(rep(rep(CheckPlots,Reps.Per.Block),blocks),ncol = blocks)
#' newmat = rbind(chk,newblocks)
#' Trts = apply(newmat, 2, sample)
#' Treatments  <- matrix(Trts,ncol=1,byrow = FALSE)
#' Treatments <- as.numeric(as.factor(Treatments))
#' trt <- length(Treatments)
#' DesLayout(matdf,trt,cb,rb,blocks) # Displays experimental layout
#'
#' ans1 = unequal.Augmented.VarCov(matdf,rb,cb,rhox=0.6,rhoy=0.6,
#' criteria="A")
#' attributes(ans1)
#' ans1$traceI
#'
#' @export
#' @seealso \code{\link{VarCov.rcbd}},\code{\link{unequal.VarCov}}
#'
#'
unequal.Augmented.VarCov <- function(matdf,rb,cb,rhox=0,rhoy=0,VarG=0.3,nugget=0,criteria="A",
                                     Amat=FALSE,sigBl=FALSE){
  X <- matrix(1,nrow = nrow(matdf))
  X<-Matrix(X)
  bb <- length(unique(matdf[,"Reps"]))
  if(is.numeric(sigBl))
  {
    Binv <- (1/sigBl)*Matrix::Diagonal(bb)
    s2e <- (1 - nugget) * (1 - VarG - sigBl)
  }else{
    sigBl <- 0.2*(1 - VarG)
    s2e <- (1 - nugget) * (1 - VarG - sigBl)
    Binv <- (1/sigBl)*Matrix::Diagonal(bb)
  }
  stopifnot(s2e > 0)
  m = length(unique(matdf[,"Treatments"]))
  if(is.matrix(Amat)){
    Gg <- VarG * as.matrix(Amat)
    Gg <- Matrix::drop0(Gg)
    Gg <- round(chol2inv(chol(Gg)),7)
  }else{
    Gg <- (1/VarG) * Matrix::Diagonal(m)
    Gg <- as(Gg, "sparseMatrix")
  }
  Ginv <- Matrix::bdiag(Binv,Gg)
  Zg <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatments"]) - 1)
  Zb <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Reps"]) - 1)
  Z <- cbind(Zb,Zg)
  if(rhox==0 & rhoy==0){
    h = nrow(matdf)/bb
    rinv = (1/s2e) * Matrix::Diagonal(h)
    Rinv <- do.call(bdiag, replicate(bb, rinv, simplify=FALSE))
  }
  else{
    matX <- subset(matdf,matdf[,"Reps"]==1)
    sigx <- Matrix::Diagonal(cb)
    sigx <- rhox^abs(row(sigx) - col(sigx))
    sigy <- Matrix::Diagonal(rb)
    sigy <- rhoy^abs(row(sigy) - col(sigy))
    R <- round(s2e * kronecker(sigy, sigx),7)
    R <- as(R, "sparseMatrix")
    Rinv <- round(chol2inv(chol(R)),7)
    Rinv <- as(Rinv, "sparseMatrix")
    Rinv <- do.call(Matrix::bdiag, replicate(bb, Rinv, simplify=FALSE))
  }

  X <- as.matrix(X)
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
