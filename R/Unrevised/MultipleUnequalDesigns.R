#' Randomly generate multiple unequally-replicated designs
#'
#' @import Matrix
#'
#' @param DesN number of initial designs to be generated
#' @param trt  number of treatments
#' @param blocks number of blocks (replicate)
#' @param max.u maximum available treatments
#' @param min.u minimum available treatments
#' @param rb number of rows per block
#' @param cb number of columns per block
#' @param Tr total number of rows
#' @param Tc total number of columns
#' @param criteria either ``A" or ``D"
#' @param VarG  heritability of the trt
#' @param rhox  spatial correlations along the x-coordinates (rows)
#' @param rhoy  spatial correlations along the y-coordinates (columns)
#' @param nugget  un-structured residual error (nugget effect)
#' @param Amat is logical, if TRUE, a numerator relationship matrix is to be provided
#' @param regular is logical, if FALSE, a different procedure is used to calculate AR1
#' @param sigBl variance of blocks. It can be specified or left blank in which case, it will be calculated internally
#' @param plot a logical argument for obtaining plots
#'
#' @return A vector of all traces, with computed mean and the best design with the smallest criterion value.
#'
#' @examples
#' trt = length(1:30);criteria="A"
#' blocks = 6; rb=5;cb=6;Tr=15;Tc=12
#' rhox=0;rhoy=0;VarG=0.1
#' nugget=0
#' min.u = rep(4,trt)
#' max.u = rep(8,trt)
#' des1 = MultipleUnequalDesigns(DesN=5, trt,blocks,max.u,
#' min.u,rb,cb,Tr,Tc, rhox,rhoy,VarG, criteria="A", Amat=TRUE,plot=TRUE)
#' attributes(des1)
#'
#' @export


MultipleUnequalDesigns <- function(DesN, trt,blocks,max.u,min.u,rb,cb,Tr,Tc, rhox=0,rhoy=0,VarG=0.3,                                   criteria="A", Amat=FALSE, nugget=0,regular=TRUE,sigBl=FALSE,plot=FALSE) {
  matrix0 <- list()
  initialValues1 <- c()
  initialValues2 <- c()
  for (i in 1:DesN) {
    print(sprintf("generating initial design: %d", i, "complete\n",
                  sep = ""))
    flush.console()
    matrix0[[i]] <- unequal.RBD(trt,blocks,max.u,min.u,rb,cb,Tr,Tc,regular)$matdf
    initialValues1[i] <- unequal.VarCov(matdf = matrix0[[i]],Tr,Tc,rhox,rhoy,VarG,nugget,criteria,Amat,sigBl,regular)[[1]]
    a <- which.min(initialValues1)
    newmatdfA <- matrix0[a][[1]]
    min_initialValues1 <- initialValues1[a][[1]]
  }
  if(plot==TRUE){
    newmatdfA = as.data.frame(newmatdfA)
    P = ggplot2::ggplot(data=newmatdfA)+geom_text(aes(x=newmatdfA$Col,
                                                      y=as.factor(newmatdfA$Row),
                                                      label=newmatdfA$Treatments,
                                                      col=as.factor(newmatdfA$Reps)), size=6) +
      scale_y_discrete("Row coordinates") +
      scale_x_discrete("Column coordinates") +
      scale_color_discrete("Block") +
      ggtitle("Initially improved design") +
      theme(text = element_text(size=15,face="bold"),
            plot.title = element_text(size=20, face="bold",vjust=2),
            axis.text=element_text(size=17,face="bold"),
            axis.title=element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      coord_fixed()
    print(P)
  }
  return(list(newmatdf = newmatdfA, initialValues = initialValues1,
              min_initialValues = min_initialValues1,
              mean_initialValues = mean(initialValues1)))
}
