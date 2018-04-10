#' Optimize a RCB Design using a pairwise swap algorithm
#'
#' \code{Optimize.ibd} optimizes an experimental design using a simple pairwise swap algorithm.
#'
#' @import Matrix
#'
#' @param matdf an  experimental layout (from function rcbd)
#' @param n the number of iterations desired
#' @param traceI a trace or determinant value calculated from matdf  ???
#' @param criteria either ``A" or ``D"
#' @param Rinv a matrix calculated for matdf (from function Rmatrix)
#' @param Ginv a matrix calculated for matdf (from function Gmatrix)
#' @param K a matrix calculated for matdf (from function VarCov.rcbd)
#' @param plot a logical argument for obtaining plots
#'
#' @return a vector of all traces from the n iterations,
#'  a  matrix with two rows: one with the accepted (successful) swaps
#'  and the other column with the iteration index,  and an optimal
#'  ``improved" design
#'
#' @examples
#' ## Example 1: optimize a single generated design
#' VarG = 0.3; rhox = 0.6; rhoy = 0.6; nugget = 0; criteria="A"
#' blocks = 3; rb = 3; cb = 5; Tr = 9; Tc = 5; trt = length(c(1:15))
#' matdf<- rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#' Rinv = Rmatrix(matdf,VarG,rhox=0.6,rhoy=0.6,nugget=0)
#' m = length(unique(matdf[,"Treatments"]))
#' Ginv <- round((1/VarG) * Matrix::Diagonal(m),7)
#'
#' res <- VarCov.rcbd(matdf,Ginv,Rinv,rhox=0.6,rhoy=0.6,VarG=0.3,criteria="A")
#' attributes(res)
#' traceI=res$OptimScore
#' K = res$K
#'
#' ans <- Optimize.rcbd(matdf=matdf,n=10,traceI=res$OptimScore,
#'        criteria="A",Rinv=Rinv, Ginv=Ginv,K=K,plot=TRUE)
#' attributes(ans)
#'
#'
#' ## TEST SAG EXAMPLE
#'
#'matdf <- rcbd(nblock=4, ntrt=30, rb=5, cb=6)
#' head(matdf)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#'
#'
#' ntrt <- 30
#' nblock <- 2; rb <- 5; cb <- 6
#' VarG <- 0.3; VarE <- 0.7; rhox <- 0.6; rhoy <- 0.6;
#' matdf <- rcbd(nblock=nblock,ntrt=ntrt,rb=rb,cb=cb)
#' head(matdf)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)

#' Rinv <- Rmatrix(matdf,VarE,rhox,rhoy)
#' Ginv <- Gmatrix(trt,VarG)
#' resD <- VarCov.rcbd(matdf,criteria="A",Ginv,Rinv)  # K is not provided but calculated
#' OptimScore <- resD$OptimScore
#' K <- resD$K
#'
#' ans <- Optimize.rcbd(matdf=matdf,n=500,OptimScore=OptimScore,
#'        criteria="A",Ginv=Ginv,Rinv=Rinv,K=K,plot=TRUE)
#'
#'
#' resD <- VarCov.rcbd(matdf,criteria="A",Ginv,Rinv,K) # K is provided
#'
#' (InitScore<-resD$OptimScore)
#' (OldScore<-resD$OptimScore)
#' Oldmatdf<-matdf
#' print(OldScore)
#' iter<-800
#' for (i in 1:iter) {
#'   #newmatdf <- SwapMethods(matdf,pairs=1,swapmethod="within")
#'   #newmatdf <- SwapMethods(matdf,pairs=4,swapmethod="across")
#'   newmatdf <- SwapMethods(Oldmatdf,pairs=2,swapmethod="any")
#'   NewScore <- VarCov.rcbd(newmatdf,criteria="A",Ginv,Rinv,K)$OptimScore # K is provided
#'   if(NewScore < OldScore){
#'     OldScore <- NewScore
#'     Oldmatdf <- newmatdf
#'     print(OldScore)
#'   }
#' }
#'
#'
#' @export
#'
#' @seealso \code{\link{rcbd}} for design of initial experiments

Optimize.rcbd<- function(matdf,n=100,OptimScore,criteria="A",Ginv,Rinv,K,plot=FALSE) {
  newmatdf <- matdf
  Score <- OptimScore
  mat <- NULL
  mat <- rbind(mat, c(value = Score, iterations = 0))
  Design_best <- newmatdf
  Des <- list()
  TRACE <- c()
  newmatdf <- SwapPair(matdf = matdf)

  for (i in 2:n) {
    newmatdf <- SwapPair(matdf = newmatdf)   ### This changes....
    TRACE[i] <- VarCov.rcbd(newmatdf,criteria,Ginv,Rinv,K)[[1]]
    Des[[i]] <- newmatdf
    if (VarCov.rcbd(newmatdf,criteria,Ginv,Rinv,K)[[1]] < Score) {
      print(sprintf("Swapping within blocks: %d", i, "complete\n",
                    sep = ""))
      Design_best <- Des[[i]] <- newmatdf
      Design_best <- newmatdf
      Score <- VarCov.rcbd(newmatdf,criteria,Ginv,Rinv,K)[[1]]
      mat <- rbind(mat, c(value = Score, iterations = i))
    }
    if (VarCov.rcbd(newmatdf,criteria,Ginv,Rinv,K)[[1]] > Score & nrow(mat) <= 1) {
      newmatdf <- matdf
      Des[[i]] <- matdf
      Design_best <- matdf
    }
    if (VarCov.rcbd(newmatdf,criteria,Ginv,Rinv,K)[[1]] > Score & nrow(mat) > 1) {
      newmatdf <- Des[[length(Des) - 1]]
      Des[[i]] <- newmatdf
      Design_best <- newmatdf
    }
  }


  if(plot==TRUE){
    Design_best = as.data.frame(Design_best)
    P = ggplot2::ggplot(data=Design_best)+geom_text(aes(x=Design_best$Col,
                                                  y=as.factor(Design_best$Row),
                                                  label=Design_best$Treatments,
                                                  col=as.factor(Design_best$Reps)), size=6) +
      scale_y_discrete("Row coordinates") +
      scale_x_discrete("Column coordinates") +
      scale_color_discrete("Block") +
      ggtitle(substitute(paste("Improved design after n = ", n,~ " iterations", sep="")))+
      theme(text = element_text(size=15,face="bold"),
            plot.title = element_text(size=20, face="bold",vjust=2),
            axis.text=element_text(size=17,face="bold"),
            axis.title=element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      coord_fixed()
    print(P)
  }

  ODE = (((mat[1,"value"]) - (mat[nrow(mat),"value"]))/(mat[1,"value"]))*100
  print(sprintf("ODE due to swapping pairs of treatments within blocks is: %f", ODE, "complete\n",
                sep = ""))
  list(TRACE = c(as.vector(mat[1, "value"]), TRACE[!is.na(TRACE)]), mat = mat,
       Design_best = Design_best)
}
