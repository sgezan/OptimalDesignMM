#' Generate an Augmented/Unreplicated/P-replicated experimental design
#'
#' \code{rcbd.Augmented} Generates an initial RCB design with replicated control treatments and unreplicated/partially replicated test treatments (Note: it only works on regular designs)
#'
#' @param blocks number of blocks (full replicate)
#' @param trt.list a numeric vector or string vector for  treatments
#' @param CheckPlots control treatments replicated in every block
#' @param Reps.Per.Block number of replications of CheckPlots per block
#' @param rb number of rows per block
#' @param cb number of columns per block
#' @param plot a logical argument for obtaining plots
#'
#' @return an experimental design in matrix form with 4 columns: Row, Col, Reps, Treatments
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2016), Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @examples
#' ## Example 1: Generates a RCB design with 4 blocks and 160 treatments + 2 checks
#' blocks <- 4;
#' test <- 160; nrep.test <- 1   # 160 treatments replicated once
#' control <- 2; nrep.control <- 5    # 2 checks replicated 5 times
#' rb <- 10; cb <- 5;
#' matdf <- rcbd.Augmented(blocks,test,nrep.test,control,nrep.control,rb,cb,plot=TRUE)
#'
#' ## Example 2: Generates a RCB design with 2 blocks and 20 treatments (double replicated) + 5 checks
#' blocks <- 2;
#' test <- 20; nrep.test <- 2   # 160 treatments replicated once
#' control <- 5; nrep.control <- 1    # 2 checks replicated 5 times
#' rb <- 5; cb <- 5;
#' matdf <- rcbd.Augmented(blocks,test,nrep.test,control,nrep.control,rb,cb,plot=TRUE)
#'
#'
#' @export
#' @seealso \code{\link{rcbd}}, \code{\link{unequal.RBD}}
#'

rcbd.Augmented <- function(blocks,test,nrep.test=1,control,nrep.control=1,rb,cb,plot=FALSE) {

  test.names <- c(1:test)
  control.names <- 1000+c(1:control)
  test.plot <- matrix(rep(test.names,nrep.test),ncol=blocks,byrow=FALSE)   # This can be improved for random asignation
  chk.plot <- matrix(rep(rep(control.names,nrep.control),blocks),ncol = blocks)
  newmat <- rbind(chk.plot,test.plot)

  Trts <- apply(newmat, 2, sample)
  Treatments <- matrix(Trts,ncol=1,byrow = FALSE)
  N <- length(Treatments)
  Reps <- rep(1:blocks,each=length(Treatments)/blocks)
  Row <- c(rep(c(rep(1:rb,each=cb)),blocks))
  Col <- c(rep(c(rep(1:cb,rb)),blocks))                     # Problems here
  cord <-cbind(Row=Row,Col=Col)
  matdf <- cbind(cord,Reps=Reps,Treatments)
  #row.names(matdf)<-NULL

  # Generating layout plot if requested
  if(plot==TRUE){
    matdf <- as.data.frame(matdf)
    P <- ggplot2::ggplot(data=matdf)+geom_text(aes(x=Col,
                                                   y=as.factor(Row),
                                                   label=Treatments,
                                                   col=as.factor(Reps)), size=6) +
      scale_y_discrete("Row coordinates") +
      scale_x_discrete("Column coordinates") +
      scale_color_discrete("Block") +
      ggtitle("Initial design") +
      theme(text = element_text(size=15,face="bold"),
            plot.title = element_text(size=20, face="bold",vjust=2),
            axis.text = element_text(size=17,face="bold"),
            axis.title = element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      coord_fixed()
    print(P)
  }
  matdf[order(matdf[,"Row"],matdf[,"Col"]),]
  row.names(matdf)<-NULL
  matdf
}


# Comments
# - Need to check for other designs and also make flexible when Tc and Tr are not provided.
# - Is it coord=list() good for initial values
# - Make an example that is read for an irregular that is in a file
# - Add Figure
# - Also, if design is regular but coordinates are provided?
