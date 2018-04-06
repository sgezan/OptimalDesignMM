#' Generate a Randomized Complete Block Design
#'
#' \code{rcbd} Generates an initial RCB design.
#'
#' @import ggplot2
#' @param blocks number of blocks (replicate)
#' @param trt number of  treatments per block
#' @param rb number of rows per block
#' @param cb number of columns per block
#' @param Tr total number of rows
#' @param Tc total number of columns
#' @param coord matrix with coordinates (x and y) and blocks for each experimental unit, in the form: Row,Col,Reps
#' @param regular logical statement, if FALSE user has to provide the coordinates and replicate for all experimental units'
#' @param plot a logical argument for obtaining plots
#'
#' @return an experimental design in matrix form with 4 variables: Rows, Cols, Blocks and Treatments
#'
#' @references
#' Mramba, L. K. and Gezan, S. A. (2016), Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
#'
#' @examples
#' ## Example 1: Generates a regular-grid RCB design with 4 blocks and 30 treatments
#' blocks <- 4; trt <- 30
#' rb <- 5; cb <- 6; Tr <- 10; Tc <- 12
#' matdf <- rcbd(blocks,trt,rb,cb,Tr,Tc,plot=TRUE)
#'
#' ## Example 2: Generates an irregular-grid RCB design with 3 blocks and 9 treatments and plots it.
#' blocks <- 3; trt <- 9;
#' Tr <- 6; Tc <- 6; rb <- 3; cb <- 3
#' Row <- c(rep(1:3,each=3),rep(4:6,each=3),rep(4:6,each=3))  # Row position of exp. units
#' Col <- c(rep(1:3,3),rep(1:3,3),rep(4:6,3)) # Column position of exp. units
#' Reps <- c(rep(1:3,each=9))
#' coord <- cbind(Row,Col,Reps)
#'
#' matdf <- rcbd(blocks,trt,rb,cb,Tr,Tc,coord,regular=FALSE,plot=TRUE)
#' DesLayout(matdf,trt,cb,rb,blocks) # Displays experimental layout
#'
#' @export
#' @seealso \code{\link{DesLayout}} for a layout with treatments in string format
#'

rcbd<- function(blocks, trt,rb,cb, Tr, Tc, coord=list(), regular=TRUE,plot=FALSE) {
  trt.list <- c(1:trt)
  trt.list <- as.numeric(as.factor(trt.list))
  Treatments <- c(replicate(blocks, sample(trt.list,trt,replace=FALSE)))

  if(blocks == 1){
    Reps <- rep(1:blocks,each=length(unique(trt.list)))
    cord <- cbind(Row=rep(1:Tr,each=Tc),Col=rep(1:Tc,Tr))
    matdf <- cbind(cord,Reps=c(rep(1,trt)),Treatments)
    row.names(matdf) <- NULL
  }
  if(cb == Tc & Tr > rb & regular==TRUE){
    Reps <- rep(1:blocks,each=length(unique(trt.list)))
    coord <- cbind(Row=rep(1:Tr,each=cb),Col=rep(1:Tc,rb))
    matdf <- cbind(coord,Reps,Treatments)
    row.names(matdf) <- NULL
  }
  if(rb == Tr & Tc > cb & regular==TRUE){
    Reps <- rep(1:blocks,each=length(unique(trt.list)))
    Row <- rep(rep(1:Tr,each=cb),Tc/cb)
    Col <- rep(split(1:Tc,cut(seq_along(1:Tc),blocks,labels=FALSE)),each=Tr)
    Col <- unlist(Col)
    coord <- cbind(Row,Col)
    matdf <- cbind(coord,Reps,Treatments)
    row.names(matdf) <- NULL
  }
  if(Tr > rb & Tc > cb & regular==TRUE){
    Reps <- rep(1:blocks,each=length(unique(trt.list)))
    Row <- rep(rep(1:Tr,each=cb),Tc/cb)
    Col <- rep(split(1:Tc,cut(seq_along(1:Tc),Tc/cb,labels=FALSE)),each=Tr)
    Col <- unlist(Col)
    cord <- cbind(Row,Col)
    matdf <- cbind(cord,Reps,Treatments)
    row.names(matdf) <- NULL
  }

  # Dealing with irregular trials
  if(regular==FALSE){
    if(length(coord) == 0){
      stop("No coordinates provided")
    }
    if(abs(nrow(coord) - length(Treatments)) > 0){
      stop("Coordinates provided do not match with size of experiment")
    }
    coord <- coord[order(coord[,3]),]  # Ordering by Reps
    matdf <- cbind(coord,Treatments)
    row.names(matdf)<-NULL
  }
  if(plot==TRUE){
    matdf = as.data.frame(matdf)
    P = ggplot2::ggplot(data=matdf)+geom_text(aes(x=Col,
                                              y=as.factor(Row),
                                              label=Treatments,
                                              col=as.factor(Reps)), size=6) +
      scale_y_discrete("Row coordinates") +
      scale_x_discrete(expand = c(0.055,0.055), "Column coordinates") +
      scale_color_discrete("Block") +
      ggtitle("Initial design") +
      theme(text = element_text(size=15,face="bold"),
            plot.title = element_text(size=20, face="bold",vjust=2),
            axis.text=element_text(size=17,face="bold"),
            axis.title=element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      coord_fixed() +
      guides(colour=FALSE)
    print(P)
  }
  matdf[order(matdf[,"Row"],matdf[,"Col"]),]
}
