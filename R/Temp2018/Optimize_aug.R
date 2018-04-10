#' Optimize an Augmented/Unreplicated/P-replicated experimental design
#'
#' \code{Optimize_aug} optimizes an initial experimental design using a simple pairwise
#' swap algorithm.
#'
#' @param matdf an experimental layout (from function augmented)
#' @param niter number of iterations desired for optimization. Default is 1.
#' @param trt.type indicates the assumptions for the factor treatment: "F" for fixed effect and
#' "R" for random effect.
#' @param block.type indicates the assumptions for the factor block: "F" for fixed effect and
#' "R" for random effect.
#' @param OptimScore optimization score criteria (A or D) from the current design matdf.
#' @param criteria for optimization, either "A" or "D". Default is "A"
#' @param Rinv matrix calculated for matdf (from function Rmatrix)
#' @param Ginv matrix calculated for matdf (from function Gmatrix)
#' @param K matrix calculated for matdf (from respective function VarCov_bFtR, VarCov_bFtF,
#' VarCov_bRtR or VarCov_bRtF)
#'
#' @return a data frame with the summary of all results (Step, OptimScore and ODE for each of
#' the succesfull swaps. Also the final improved experiment with 4 columns: Row, Col, Rep, Treatment.
#'
#' @examples
#' # Example 1: Optimizes an Augmented-RCB design with 4 blocks and 36 treatments + 3 checks
#' matdf <- augmented(nblock=4, test.trt=c(36,1), check.trt=c(3,1), rb=4, cb=3)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
#' Rinv <- Rmatrix(matdf=matdf, VarE=0.7, rhox=0.9, rhoy=0.9, regular=TRUE)
#' Ginv <- Gmatrix(ng=39, VarG=0.3)   # Independent random effects for a heritability of 0.3
#' resD <- VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv)
#' resp <- Optimize_aug(matdf=matdf, niter=1000, trt.type='F', block.type='R',
#'                      OptimScore=resD$OptimScore, criteria="A", Ginv=Ginv, Rinv=Rinv, K=resD$K)
#' head(resp$results)
#' newmatdf <- resp$matdf
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
#'
#'
#'
#' # Example 2: Optimizes a CRB design, 136 treatments + 3 checks (with 8 replications each)
#' matdf <- augmented(nblock=1, test.trt=c(136,1), check.trt=c(3,8), rb=16, cb=10)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
#' Rinv <- Rmatrix(matdf=matdf, VarE=0.7, rhox=0.9, rhoy=0.9, regular=TRUE)
#' Ginv <- Gmatrix(ng=139, VarG=0.3)   # Independent random effects for a heritability of 0.3
#' ##########
#' resD <- VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv)  # Wrong becuase it is a single block
#' resp <- Optimize_aug(matdf=matdf, iter=5000, s=1, OptimScore=resD$OptimScore,
#'        criteria="A", Ginv=Ginv, Rinv=Rinv, K=resD$K)
#' head(resp$results)
#' newmatdf <- resp$matdf
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
#'
#'
#'#' # Test for Varcov_bFtR for CRD
#' matdf <- crd(nrep=3, ntrt=12, rb=6, cb=6)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#' Rinv <- Rmatrix(matdf=matdf, VarE=0.7, rhox=0.6, rhoy=0.6, regular=TRUE)
#' Ginv <- Gmatrix(ng=12, VarG=0.3)   # Independent random effects for a heritability of 0.3
#' resD <- VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv)  # K is not provided but calculated
#' resD$OptimScore
#'
#' @export
#'
#' @seealso \code{\link{augmented}} for design of initial experiments

Optimize_aug<- function(matdf, niter=1, trt.type='R', block.type='F', OptimScore,
                        criteria="A", s2Bl=NA, Ginv=NULL, Rinv=NULL, K=NULL) {

  # Create a table to store results
  results <- data.frame(Step=0, Score=OptimScore, ODE=0)
  print(results)
  OldScore <- OptimScore
  Oldmatdf <- matdf

  for (i in 1:niter) {

     newmatdf <- SwapMethods(matdf=Oldmatdf, pairs=1, swapmethod="within")

     if (block.type == "F" & trt.type == "F") {
         resD <- VarCov_bFtF(matdf=newmatdf, criteria=criteria)
     } else if (block.type == "R" & trt.type == "F") {
         if(is.null(s2Bl)){ stop('Variance of Blocks, s2Bl, was not provided') }
         resD <- VarCov_bRtF(matdf=newmatdf, criteria=criteria, s2Bl=s2Bl, Rinv=Rinv, K=K)
     } else if (block.type == "F" & trt.type == "R") {
         resD <- VarCov_bFtR(matdf=newmatdf, criteria=criteria, Ginv=Ginv, Rinv=Rinv, K=K)
     } else if (block.type == "R" & trt.type == "R") {
         if(is.null(s2Bl)){ stop('Variance of Blocks, s2Bl, was not provided') }
         resD <- VarCov_bRtR(matdf=newmatdf, criteria=criteria, s2Bl=s2Bl, Ginv=Ginv, Rinv=Rinv, K=K)
     }

     NewScore <- resD$OptimScore
     if(NewScore < OldScore){
         ODEtemp <- 100*((OptimScore - NewScore)/OptimScore)
         OldScore <- NewScore
         Oldmatdf <- newmatdf
         temp <- c(i, OldScore, ODEtemp)
         print(temp)
         results <- rbind(results, temp)
     }
  }
  results <- data.frame(results)
  return(list(results=results,matdf=newmatdf))
}


