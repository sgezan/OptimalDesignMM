#' Generates error surfaces
#'
#' \code{Error.Surfaces} Generates error surfaces to be incorporated in the simulation of the phenotype
#' @import Matrix
#' @param Tr total number of rows
#' @param Tc total number of columns
#' @param rhox spatial correlation along the rows
#' @param rhoy spatial correlation along the columns
#' @param VarG variance of the treatment effects.
#' @param nugget  un-structured residual error (nugget effect)
#'
#' @details This function is needed for internal purposes only to generate error
#' surfaces for incorporation into the simulated phenotype
#'
#' @return Returns a vector of error surfaces
#'
#' @references
#' Mramba, L. K. and Gezan, S. A. (2016), Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
#'
#' @examples
#' blocks = 4; Tr=10;Tc=12;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#' surf = Error.Surfaces(Tr,Tc)
#'
#' @export
#' @seealso  \code{\link{rcbd}}, \code{\link{VarCov.rcbd}}
#'

Error.Surfaces<-function(Tr,Tc,rhox=0,rhoy=0,VarG=0.3,nugget=0){
  # calculation of base V(N*N) matrix with AR(1)*AR(1) covariance model
  e1=rnorm(Tr*Tc);e2=rnorm(Tr*Tc)   # e1,e2: independent vectors of errors
  N = Tr*Tc
  V = diag(N)
  sigx <- Matrix::Diagonal(Tc)
  sigx <- rhox^ abs(row(sigx) - col(sigx))
  sigy <- Matrix::Diagonal(Tr)
  sigy <- rhoy ^ abs(row(sigy) - col(sigy))
  s2e <- 1-VarG-nugget
  V<- s2e * kronecker(sigy, sigx)  # takes 0.01 second
  L = chol(V)
  PAT = t(L)%*%e1
  PAT = (PAT - mean(PAT))/sd(PAT) # standardized patches
  ZST = sqrt(1-VarG)*(PAT*sqrt(1-nugget)+e2*sqrt(nugget))
  ZST
}
