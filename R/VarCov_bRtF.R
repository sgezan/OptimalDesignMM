#' Generate variance-covariance matrix of treatments, traces and/or log(determinants) and other
#' related matrices for block designs (RCBD or IBD) with random block effects and fixed
#' treatment effects.
#'
#' \code{VarCov_bRtF} generates the variance-covariance matrix (information matrix) of treatment
#' effects and calculates traces and/or log(determinants) for a linear mixed model with
#' random block effects and fixed treatment effects.
#'
#' @param matdf an experimental design (layout) based on a randomized complete block designs (RCBD),
#' or incomplete block designs (IBD) where 'Treatment' is the column of effects of interest.
#' @param criteria indicates the optimization criteria to report. It can be 'A' for A-optimal or
#' 'D' for D-optimal criteria. Default is 'A'.
#' @param s2Bl variance of the blocks effect. Default is 0.1
#' @param Rinv an inverse of the error variance-covariance matrix.
#' @param K an intermediate matrix calculated from the original layout that was previously obtained.
#'
#' @return either a trace value or log of determinant of the variance-covariance matrix
#' (information matrix) for treatments, together with the K matrix to use in future runs.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2016) Generating experimental designs for spatially and
#' genetically correlated data using mixed models, Submitted to Australian and New Zealand
#' Journal of Statistics.
#'
#' @examples
#' # Example 1: A-optimality
#' matdf <- rcbd(nblock=2, ntrt=9, rb=3, cb=3)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#' Rinv <- Rmatrix(matdf=matdf, VarE=0.7, rhox=0.6, rhoy=0.6, regular=TRUE)
#' resD <- VarCov_bRtF(matdf=matdf, criteria="A", s2Bl=0.15, Rinv=Rinv)
#' attributes(resD)
#' resD$OptimScore
#' resD <- VarCov_bRtF(matdf=matdf, criteria="A", s2Bl=0.15, Rinv=Rinv, K=resD$K) # K is provided
#'
#' # Example 2: D-optimality
#' VarCov_bRtF(matdf=matdf, criteria="D", s2Bl=0.15, Rinv=Rinv, K=resD$K)$OptimScore
#'
#' @export
#' @seealso \code{\link{rcbd}}

VarCov_bRtF <- function(matdf, criteria="A", s2Bl=0.1, Rinv=NULL, K=NULL) {

  if(nrow(matdf)==length(unique(matdf[,"Treatment"]))){
    X <- as.matrix(matdf[, "Treatment"])
    colnames(X)<-NULL
  }
  if(nrow(matdf) > length(unique(matdf[,"Treatment"]))){
    X <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatment"]) - 1)
    colnames(X)<-NULL
  }

  # Obtaining Ginverse for random blocks
  b <- length(unique(matdf[,"Rep"]))  # Number of blocks
  Binv <- Gmatrix(ng=b, VarG=s2Bl)   # Independent random effects with variance s2Bl

  Z <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Rep"]) - 1)
  colnames(Z) <- NULL
  Z <-  as.matrix(Z)
  X <-  as.matrix(X)
  Rinv <- as.matrix(Rinv)

  # A different approach
  C11 <- Matrix::crossprod(X, Rinv) %*% X
  C11  <- as(C11, "sparseMatrix")
  #C22 <- t(Z) %*% Rinv  %*% Z + Binv
  #C22inv <- solve(C22)
  #C12 <- t(X) %*% Rinv  %*% Z
  #temp0 <- C11 - C12 %*% C22inv %*% t(C12)
  #M <- solve(temp0)

  if(is.null(K)){
    C22 <- t(Z) %*% Rinv  %*% Z + Binv
    C22inv <- solve(C22)
    K <- Rinv %*% Z %*% C22inv %*% Matrix::crossprod(Z, Rinv)
  }

  #K <- round(K,7)
  #K <- Matrix::drop0(K)  # Eliminationg zeros

  #temp0 <- C11 - t(X) %*% K %*% X
  temp0 <- C11 - Matrix::crossprod(X, K) %*% X

  M <- solve(temp0)
  #M <- round(M,7)
  #M <- as(M, "sparseMatrix")

  # Calculating Optimum Criteria over matrix M (of ALL fixed treatment effects)
  if(criteria == "A"){
    OptimScore <- sum(Matrix::diag(M))
    return(list(OptimScore=OptimScore,K=K))
  }
  if(criteria == "D"){
    OptimScore <- log(Matrix::det(M))
    return(list(OptimScore=OptimScore,K=K))
  }
}

# Need to check that ther important matrix is M not C22 (in reality C11)
