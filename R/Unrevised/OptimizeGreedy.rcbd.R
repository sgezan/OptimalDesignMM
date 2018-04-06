#' Optimize a RCB Design using a greedy algorithm
#'
#' \code{OptimizeGreedy.rcbd} optimizes an RCB experimental design using a greedy swap algorithm.
#'
#' @import Matrix
#'
#' @param matdf  an  experimental layout
#' @param n the number of iterations required
#' @param traceI a trace or determinant value calculated from matdf;
#' @param criteria either ``A" or ``D"
#' @param gsize number of treatments to be swapped per iteration. This should be an even number
#' @param Rinv a matrix calculated from matdf
#' @param Ginv a matrix calculated from matdf
#' @param K a matrix calculated from matdf
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
#' ans <- OptimizeGreedy.rcbd(matdf,n=10,traceI,criteria,gsize=4,Rinv,Ginv,K,plot=TRUE)
#' attributes(ans)
#'
#' @export
#'
#' @seealso \code{\link{rcbd}} for design of initial experiments

OptimizeGreedy.rcbd<- function(matdf,n=10,traceI,criteria="A",gsize=4,Rinv,Ginv,K,plot=FALSE) {
  newmatdf <- matdf
  trace <- traceI
  mat <- NULL
  mat <- rbind(mat, c(value = trace, iterations = 0))
  Design_best <- newmatdf
  Des <- list()
  TRACE <- c()
  newmatdf <- SwapGreedy(matdf = matdf,gsize = gsize)
  for (i in 2:n) {
    newmatdf <- SwapGreedy(matdf = newmatdf,gsize = gsize)
    TRACE[i] <- VarCov.rcbd(newmatdf,Ginv,Rinv,criteria,K = as.matrix(K), Update=TRUE)[[1]]
    Des[[i]] <- newmatdf
    if (VarCov.rcbd(newmatdf,Ginv,Rinv,criteria,K = as.matrix(K), Update=TRUE)[[1]] < trace) {
      print(sprintf("Swapping greedly within blocks: %d", i, "complete\n",
                    sep = ""))
      Design_best <- Des[[i]] <- newmatdf
      Design_best <- newmatdf
      trace <- VarCov.rcbd(newmatdf,Ginv,Rinv,criteria,K = as.matrix(K), Update=TRUE)[[1]]
      mat <- rbind(mat, c(trace = trace, iterations = i))
    }
    if (VarCov.rcbd(newmatdf,Ginv,Rinv,criteria,K = as.matrix(K), Update=TRUE)[[1]] > trace & nrow(mat) <= 1) {
      newmatdf <- matdf
      Des[[i]] <- matdf
      Design_best <- matdf
    }
    if (VarCov.rcbd(newmatdf,Ginv,Rinv,criteria,K = as.matrix(K), Update=TRUE)[[1]] > trace & nrow(mat) > 1) {
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
  print(sprintf("ODE due to greedly swapping pairs of treatments within blocks is: %f", ODE, "complete\n",
                sep = ""))
  list(TRACE = c(as.vector(mat[1, "value"]), TRACE[!is.na(TRACE)]), mat = mat,
       Design_best = Design_best)
}
