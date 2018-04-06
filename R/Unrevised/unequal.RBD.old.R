#' Generates a randomized block experiment
#'
#' \code{unequal.RBD} generates an initial design for unequally
#' replicated, incomplete block and and augmented designs with regular or irregular grids
#'
#' @param trt  number of treatments
#' @param blocks number of blocks (replicate)
#' @param max.u maximum replicates allowed per treatment
#' @param min.u minimum replicates allowed per treatment
#' @param rb number of rows per block
#' @param cb number of columns per block
#' @param Tr total number of rows
#' @param Tc total number of columns
#' @param regular a logical statement, if FALSE, user has to provide the x and y cords for the experiment.
#' @param plot a logical argument for obtaining plots
#'
#' @return A list of objects, with an experimental design among others
#'
#' @references
#' Mramba, Lazarus. K. and Gezan, Salvador. A. (2016), Improving Unequally Replicated,
#' Incomplete Block and Augmented Experimental Designs with Spatially
#' and Genetically Correlated Observations, Submitted to the Journal of
#' Theoretical and Applied Genetics
#'
#' @examples
#' trt = 30;criteria="A"
#' blocks = 3; rb=5;cb=6;Tr=5;Tc=18;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#' min.u = sample(1:3,trt,replace=TRUE)
#' max.u = sample(3:5,trt,replace=TRUE)
#' des1 = unequal.RBD(trt,blocks,max.u,min.u,rb,cb,Tr,Tc,plot=TRUE)
#' attributes(des1)
#' des1$sumFreq
#' head(des1$matdf)
#' head(des1$datam)
#' DesLayout(matdf=des1$matdf,trt, cb, rb, blocks)
#'
#' @export
#' @seealso \code{\link{rcbd}}
#'
unequal.RBD<-function(trt,blocks,max.u,min.u,rb,cb,Tr,Tc,regular=TRUE, plot=FALSE)
{
  if (all(max.u == 1))
    stop("Number of replicates have to differ !")
  stopifnot(max.u >= min.u, min.u > 0)
  trt = c(1:trt)
  stopifnot(is.numeric(blocks), blocks > 0, length(trt) > 1 )
  N <- rb * cb * blocks
  genot <- as.numeric(as.factor(trt))
  genot1  <- rep(genot,min.u,replace=FALSE)
  genot1  <- rep(genot,min.u,replace=FALSE)
  n = length(genot1)
  stopifnot(n <= N)
  if(n == N) Treatments <- sample(genot1 ,N, replace = FALSE)
  if(n < N){
    freq = as.numeric(table(genot1))
    Diff = max.u - min.u
    datam = cbind(genot,min.u,max.u,freq,Diff)
    u.extra = datam[freq< max.u,"genot"]
    extraG =rep(u.extra,Diff[Diff!=0])
    sampleG <- sample(extraG,N-n)
    Treatments<-c(genot1,sampleG)
    Treatments <- sample(Treatments,N,replace=FALSE)
  }
  Treatments <- factor(Treatments)
  Reps <- rep(1:blocks,each=N/blocks)
  freq = as.numeric(table(Treatments))
  datam = cbind(genot,min.u,max.u,freq)

  if(regular==TRUE & cb == Tc & Tr > rb){
    cord <- cbind(Row=rep(1:Tr,each=cb), Col=rep(1:Tc,rb))
    matdf <- cbind(cord,Reps, Treatments)
    row.names(matdf)<-NULL
  }
  if(regular==TRUE & rb == Tr & Tc > cb){
    Row=rep(rep(1:Tr,each=cb),Tc/cb)
    Col <- rep(split(1:Tc, cut(seq_along(1:Tc), blocks,
                               labels = FALSE)),each=Tr)
    Col <- unlist(Col)
    cord <- cbind(Row,Col)
    matdf <- cbind(cord,Reps, Treatments)
    row.names(matdf)<-NULL
  }
  if(regular==TRUE  & Tr > rb & Tc > cb){
    Row=rep(rep(1:Tr,each=cb),Tc/cb)
    Col <- rep(split(1:Tc, cut(seq_along(1:Tc), Tc/cb,
                               labels = FALSE)),each=Tr)
    Col <- unlist(Col)
    cord <- cbind(Row,Col)
    matdf <- cbind(cord,Reps, Treatments)
    row.names(matdf)<-NULL
  }
  if(regular==FALSE){
    cord <-cbind(Row=Row,Col=Col)
    matdf <- cbind(cord, Reps, Treatments)
    row.names(matdf)<-NULL
  }
  if(plot==TRUE){
    matdf = as.data.frame(matdf)
    P = ggplot2::ggplot(data=matdf)+geom_text(aes(x=Col,
                                                  y=as.factor(Row),
                                                  label=Treatments,
                                                  col=as.factor(Reps)), size=6) +
      scale_y_discrete("Row coordinates") +
      scale_x_discrete("Column coordinates") +
      scale_color_discrete("Block") +
      ggtitle("Initial design") +
      theme(text = element_text(size=15,face="bold"),
            plot.title = element_text(size=20, face="bold",vjust=2),
            axis.text=element_text(size=17,face="bold"),
            axis.title=element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      coord_fixed()
    print(P)
  }
  return(list(matdf = matdf[order(matdf[,"Row"],matdf[,"Col"]),],
              datam = datam, sumFreq = sum(freq)))
}
