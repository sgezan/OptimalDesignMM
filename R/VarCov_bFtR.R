#' Generate variance-covariance matrix of treatments, traces and/or log(determinants) and other
#' related matrices for a RCBD with fixed blocks effects and random treatment effects.
#'
#' \code{VarCov_bFtR} generates the variance-covariance matrix (information matrix) of treatment
#' effects and calculates traces and/or log(determinants) for a linear mixed model with
#' fixed block effects and random treatment effects.
#'
#' @param matdf an experimental design (layout) based on a randomized complete block designs (RCBD),
#' where 'Treatment' is the column of effects of interest
#' @param criteria indicates the optimization criteria to report. It can be 'A' for A-optimal or
#' 'D' for D-optimal criteria. Default is 'A'.
#' @param Ginv a variance-covariance matrix from pedigree or molecular data previously generated
#' @param Rinv an inverse of the error variance-covariance matrix.
#' @param K an intermediate matrix calculated from the original layout that was previously obtained.
#'
#' @return either a trace value or log of determinant of the variance-covariance matrix
#' (information matrix) for treatments, together with the K matrix to use in future runs.
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2016) Generating experimental designs for spatially and
#' genetically correlated data using mixed models, Submitted to Australian and New Zealand
#' Journal of Statistics.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @examples
#' # Example 1: Regular-grid experiment with independent random effects
#' matdf <- rcbd(nblock=2, ntrt=30, rb=5, cb=6)
#' Rinv <- Rmatrix(matdf=matdf, VarE=0.7, rhox=0.6, rhoy=0.6, regular=TRUE)
#' Ginv <- Gmatrix(ng=30, VarG=0.3)   # Independent random effects for a heritability of 0.3
#' resD <- VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv)  # K is not provided but calculated
#' resD$OptimScore
#' VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv, K=resD$K)$OptimScore # K is provided
#'
#' # Example 2: Regular-grid experiment with 2 blocks and 30 related genotypes (treatments)
#' matdf <- rcbd(nblock=2, ntrt=30, rb=5, cb=6)
#' data(ped30fs)
#' Amat <- makeA(ped30fs)
#' Amat <- as.matrix(Amat[-c(1:5),-c(1:5)])   # Only offspring, eliminating the 5 parents.
#' Ginv <- Gmatrix(VarG=1, G=Amat)
#' Rinv <- Rmatrix(matdf=matdf, VarE=0.7, rhox=0.6, rhoy=0.6, regular=TRUE)
#' resPed <- VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv)
#' resPed$OptimScore
#'
#' @export
#' @seealso \code{\code{\link{rcbd} \link{Rinv}}, \code{\link{Ginv}}

VarCov_bFtR <- function(matdf, criteria="A", Ginv=NULL, Rinv=NULL, K=NULL) {

  # Checking if Rinv and Ginv are not provided
  if(is.null(Rinv)){
    stop('Rinv was not provided')
  }
  if(is.null(Ginv)){
    stop('Ginv was not provided')
  }

  # Obtaining Z matrix
  Z <- Matrix::sparse.model.matrix(~as.factor(matdf[,"Treatment"]) - 1)
  Z <- as.matrix(Z)

  # Obtaining Rinv and Ginv matrix (and its inverse) for spatial analysis
  Rinv <- Matrix::drop0(round(Rinv,7))
  Ginv <- Matrix::drop0(round(Ginv,7))

  if(is.null(K)){

    # Obtaining X matrix
    if(nrow(matdf) == length(unique(matdf[,"Treatment"]))){
      X <- as.matrix(matdf[, "Rep"])
    }
    if(nrow(matdf) > length(unique(matdf[,"Treatment"]))){
      if (length(unique(matdf[,"Rep"])) == 1) {
        X <- matrix(data=1,nrow=nrow(matdf),ncol=1)
      } else {
        X <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Rep"])-1)
      }
    }

    # Obtaining C22 (only for Random Treatments)
    C11 <- Matrix::crossprod(as.matrix(X), as.matrix(Rinv)) %*% as.matrix(X)
    C11inv <- solve(C11)
    k1 <- Rinv %*% as.matrix(X)
    k2 <- Matrix::tcrossprod(as.matrix(C11inv), as.matrix(X))
    k3 <- k2 %*% Rinv
    K <- k1 %*% k3
    K <- as(K, "sparseMatrix")
    temp0 <- Matrix::crossprod(Z, Rinv) %*% Z + Ginv - Matrix::crossprod(Z, K) %*% Z

    # Rounding K matrix
    K <- round(K,7)

  } else {
    temp0 <- t(Z) %*% Rinv %*% Z + Ginv - t(Z) %*% K %*% Z
  }

  C22 <- solve(temp0)
  C22 <- as(C22, "sparseMatrix")  # This is the M(lambda) matrix

  # Calculating Optimum Criteria (over matrix C22 of ALL random effects)
  if(criteria == "A"){
    OptimScore <- sum(Matrix::diag(C22))
    return(list(OptimScore=OptimScore,K=K))
  }
  if(criteria == "D"){
    OptimScore <- log(Matrix::det(C22))
    return(list(OptimScore=OptimScore,K=K))
  }

}
