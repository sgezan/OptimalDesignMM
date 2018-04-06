#' Improves an unequally replicated design by replacing treatments using a list
#'
#' \code{unequal.Optimize.Rpl} improves an experimental design by replacing treatments within blocks
#'
#'
#' @param matdf  an  experimental layout
#' @param n the number of iterations required
#' @param traceI a trace or determinant value calculated from matdf;
#' @param criteria either ``A" or ``D"
#' @param Rinv a matrix calculated from matdf
#' @param Ginv a matrix calculated from matdf
#' @param K a matrix calculated from matdf
#' @param trt a numeric vector or string vector for  treatments
#' @param max.u the maximum number treatments available
#' @param min.u the minimum number of treatments allowed
#' @param plot a logical argument for obtaining plots
#'
#' @return a vector of all traces from the n iterations,
#'  a  matrix with two rows: one with the accepted (successful) swaps
#'  and the other column with the iteration index,
#'  and an improved design
#'
#' @export
#'
#' @examples
#' trt = length(1:30);criteria="A"
#' blocks = 3; rb=5;cb=6;Tr=5;Tc=18;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#' min.u = sample(1:3,trt,replace=TRUE)
#' max.u = sample(3:5,trt,replace=TRUE)
#'
#' des1 = unequal.RBD(trt,blocks,max.u,min.u,rb,cb,Tr,Tc,plot=TRUE)
#'
#' data(ped30hs)
#' Amat <- GenA(male=ped30hs[,"male"], female = ped30hs[,"female"])
#' Amat <- as.matrix(Amat[-c(1:5), -c(1:5)])
#' ans1 <- unequal.VarCov(des1$matdf,Tr,Tc, rhox=0.6,rhoy=0.6,
#'            VarG=0.3,nugget=0, criteria="A",Amat=TRUE)
#' attributes(ans1)

#' traceI <- ans1$traceI
#' Rinv = as.matrix(ans1$Rinv); Ginv = as.matrix(ans1$Ginv); K = as.matrix(ans1$K)
#'
#' Results <- unequal.Optimize.Rpl(matdf=des1$matdf,n=10,
#' traceI,criteria="A",Rinv,Ginv,K,trt,max.u,min.u,plot=TRUE)
#'
#' @seealso \code{\link{SwapPair}}, \code{\link{unequal.NewValue}}
#'
#'
unequal.Optimize.Rpl<- function(matdf,n,traceI,criteria="A",Rinv,Ginv,K,trt,max.u,min.u,plot=FALSE)
{
  newmatdf <- matdf
  trace <- traceI
  mat <- NULL
  mat <- rbind(mat, c(value = trace, iterations = 0))
  Design_best <- newmatdf
  Des <- list()
  newmatdf <- ReplaceTrts(matdf=matdf,trt,max.u,min.u)$matdf
  TRACE <- c()
  for (i in 2:n) {
    newmatdf <- ReplaceTrts(matdf=newmatdf,trt,max.u,min.u)$matdf
    TRACE[i] <- unequal.NewValue(matdf=newmatdf,  criteria, Rinv, Ginv, K)
    Des[[i]] <- newmatdf
    if (unequal.NewValue(matdf=newmatdf, criteria, Rinv, Ginv, K) < trace) {
      print(sprintf("Replacing treatments: %d", i, "complete\n",
                    sep = ""))
      Design_best <- Des[[i]] <- newmatdf
      Design_best <- newmatdf
      trace <- unequal.NewValue(matdf=newmatdf, criteria, Rinv, Ginv, K)
      mat <- rbind(mat, c(trace = trace, iterations = i))
    }
    if (unequal.NewValue(matdf=newmatdf, criteria, Rinv, Ginv, K) > trace &
        nrow(mat) <= 1) {
      newmatdf <- matdf
      Des[[i]] <- matdf
      Design_best <- matdf
    }
    if (unequal.NewValue(matdf=newmatdf, criteria, Rinv, Ginv, K) > trace &
        nrow(mat) > 1) {
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
  ODE = ((mat[1,"value"] - mat[nrow(mat),"value"])/mat[1,"value"])*100
  print(sprintf("ODE due to replacing treatments is: %f", ODE, "complete\n",
                sep = ""))
  list(ODE = ODE, TRACE = c(as.vector(mat[1, "value"]), TRACE[!is.na(TRACE)]), mat = mat,
       Design_best = Design_best)
}
