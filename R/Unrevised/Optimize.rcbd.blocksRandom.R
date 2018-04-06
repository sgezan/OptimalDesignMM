#' Optimize a RCB Design with random  blocks and fixed treatments using a pairwise swap algorithm
#'
#' @import Matrix
#'
#' @param matdf  an  experimental layout
#' @param n the number of iterations required
#' @param traceI a trace or determinant value calculated from matdf;
#' @param criteria either ``A" or ``D"
#' @param K a matrix calculated from matdf
#' @param Rinv is an inverse of the spatial correlation matrix
#' @param plot a logical argument for obtaining plots
#'
#' @return a vector of all traces from the n iterations,
#'  a  matrix with two rows: one with the accepted (successful) swaps
#'  and the other column with the iteration index,  and an optimal
#'  ``improved" design
#'
#' @examples
#' # Example
#' h2 = 0.3; rhox = 0.6; rhoy = 0.6; s20 = 0; criteria="A"
#' blocks = 3; rb = 3; cb = 3; Tr = 9; Tc = 3; trt = length(c(1:9))
#'
#' matdf<- rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#' res <- VarCov.rcbd.blocksRandom(matdf,Tr,Tc,criteria="A",rhox = 0.6, rhoy = 0.6,h2 = 0.3)
#'
#' attributes(res)
#' traceI=res$traceI; criteria="A"
#' Rinv = as.matrix(res$Rinv); K = as.matrix(res$K)
#'
#' ans <- Optimize.rcbd.blocksRandom(matdf,n=10,traceI=traceI,criteria,K=K,Rinv=Rinv,plot=TRUE)
#' attributes(ans)
#'
#' @export
#'
#' @seealso \code{\link{rcbd}} for design of initial experiments
#'
Optimize.rcbd.blocksRandom <- function(matdf,n=10,traceI,criteria="A",K,Rinv,plot=FALSE) {
  newmatdf <- matdf
  trace <- traceI
  mat <- NULL
  mat <- rbind(mat, c(value = trace, iterations = 0))
  Design_best <- newmatdf
  Des <- list()
  TRACE <- c()
  newmatdf <- SwapPair(matdf = matdf)
  for (i in 2:n) {
    newmatdf <- SwapPair(matdf = newmatdf)
    TRACE[i] <- NewValue.rcbd.blocksRandom(matdf=newmatdf, criteria, K, Rinv)
    Des[[i]] <- newmatdf
    if (NewValue.rcbd.blocksRandom(matdf=newmatdf, criteria, K, Rinv) < trace) {
      print(sprintf("Swapping within blocks: %d", i, "complete\n",
                    sep = ""))
      Design_best <- Des[[i]] <- newmatdf
      Design_best <- newmatdf
      trace <- NewValue.rcbd.blocksRandom(matdf=newmatdf, criteria,  K, Rinv)
      mat <- rbind(mat, c(trace = trace, iterations = i))
    }
    if (NewValue.rcbd.blocksRandom(matdf=newmatdf, criteria, K, Rinv) > trace & nrow(mat) <= 1) {
      newmatdf <- matdf
      Des[[i]] <- matdf
      Design_best <- matdf
    }
    if (NewValue.rcbd.blocksRandom(matdf=newmatdf, criteria,  K, Rinv) > trace & nrow(mat) > 1) {
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
