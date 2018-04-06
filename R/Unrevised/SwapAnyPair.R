#' Randomly Swap any Pair of Treatments either within or across blocks
#'
#' \code{unequal.SwapAcross} Randomly swaps any pairs of treatments either within or across blocks
#'
#' @param matdf an xperimental design
#' @param plot a logical statement, if true, a plot is printed
#'
#' @return an experimental design with swapped pairs of treatments
#'
#' @references
#' Mramba, Lazarus. K. and Gezan, Salvador. A. (2016), Improving Unequally Replicated,
#' Incomplete Block and Augmented Experimental Designs with Spatially
#' and Genetically Correlated Observations, Submitted to the Journal of
#' Theoretical and Applied Genetics
#'
#' @examples
#' blocks = 4; trt = length(1:16); Tr = 8; Tc = 8; rb = 4; cb = 4
#' matdf <- rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#' newdes <- SwapAnyPair(matdf,plot=TRUE)
#' which(matdf[,"Treatments"] != newdes[,"Treatments"])
#'
#' @export
#' @seealso \code{\link{SwapPair}}, \code{\link{ReplaceTrts}}, \code{\link{unequal.SwapAcross}}
#'
SwapAnyPair <- function(matdf,plot=FALSE) {
  mat0 <- matdf[,"Treatments"]
  sam <- sample(mat0,2, replace = TRUE)
  mat0[c(sam[1],sam[2])] <- mat0[c(sam[2],sam[1])]
  matdf[,"Treatments"] <- mat0
  matdf[order(matdf[,"Row"],matdf[,"Col"]),]
  if(plot==TRUE){
    matdf = as.data.frame(matdf)
    P = ggplot2::ggplot(data=matdf)+geom_text(aes(x=matdf$Col,
                                                  y=as.factor(matdf$Row),
                                                  label=matdf$Treatments,
                                                  col=as.factor(matdf$Reps)), size=6) +
      scale_y_discrete("Row coordinates") +
      scale_x_discrete("Column coordinates") +
      scale_color_discrete("Block") +
      ggtitle("Swapping within blocks") +
      theme(text = element_text(size=15,face="bold"),
            plot.title = element_text(size=20, face="bold",vjust=2),
            axis.text=element_text(size=17,face="bold"),
            axis.title=element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      coord_fixed()
    print(P)
  }
  return( matdf[order(matdf[,"Row"],matdf[,"Col"]),])
}

