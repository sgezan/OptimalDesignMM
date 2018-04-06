#' Generate variance-covariance matrix of treatments, traces and/or log(determinants) and other
#' related matrices for an RCBD with both block and treatment fixed effects.
#'
#' \code{VarCov_bRtF} generates the variance-covariance matrix (information matrix) of treatment
#' effects and calculates traces and/or log(determinants) for a linear model with both
#' block and treatment fixed effects, from a provided RCDB layout.
#'
#' @param matdf an experimental design (layout) based on a randomized complete block designs (RCBD),
#' where 'Treatment' is the column of effects of interest.
#' @param criteria indicates the optimization criteria to report. It can be 'A' for A-optimal or
#' 'D' for D-optimal criteria. Default is 'A'.
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
#' Example: Regular-grid experiment with independent random effects
#' matdf <- rcbd(nblock=2, ntrt=30, rb=5, cb=6)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#' Rinv <- Rmatrix(matdf=matdf, VarE=1, rhox=0.6, rhoy=0.6, regular=TRUE)
#' resD <- VarCov_bFtF(matdf=matdf, criteria="A")  # K is not provided but calculated
#' resD$OptimScore
#' VarCov_bFtF(matdf=matdf, criteria="A", K=resD$K)$OptimScore # K is provided
#'
#' @export
#' @seealso \code{\link{rcbd}}

VarCov_bFtF <- function(matdf, criteria="A", Rinv=NULL) {

  # X matrix of Reps, W matrix of treatments
  b <- length(unique(matdf[,"Rep"]))  # Number of blocks
  t <- length(unique(matdf[,"Treatment"]))  # Number of blocks

  bb <- as.factor(matdf[, "Rep"])
  gg <- as.factor(matdf[, "Treatment"])
  if(nrow(matdf)==length(unique(matdf[,"Treatment"]))){
    X <- as.matrix(matdf[, "Rep"])
    colnames(X)<-NULL
  }
  if(nrow(matdf) > length(unique(matdf[,"Treatment"]))){
    X <- Matrix::sparse.model.matrix(~bb-1)
    X[,2] <- X[,2]-X[,1]
    X <- X[,2:max(matdf[, "Rep"])]
    colnames(X)<-NULL
  }
  W <- Matrix::sparse.model.matrix(~gg-1)
  colnames(W)<-NULL

  # Defining some of the matrices
  X <- as.matrix(X)
  W <- as.matrix(W)
  Rinv <- as.matrix(Rinv)
  Xf <- cbind(X,W)

  C11 <- Matrix::crossprod(Xf, Rinv) %*% Xf
  C11inv <- solve(C11)
  M <- C11inv[-c(1:(b-1)),-c(1:(b-1))]
  #M <- as(M, "sparseMatrix")  # This is the M(lambda) matrix for all random effects

  # Calculating Optimum Criteria over matrix M (of ALL fixed treatment effects)
  if(criteria == "A"){
    OptimScore <- sum(Matrix::diag(M))
    return(list(OptimScore=OptimScore))
  }
  if(criteria == "D"){
    OptimScore <- log(Matrix::det(M))
    return(list(OptimScore=OptimScore))
  }
}
