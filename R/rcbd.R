#' Generates a Randomized Complete Block Design (RCBD)
#'
#' \code{rcbd} Randomly generates a randomized complete block design. Note that for regular
#' grids blocks are assumed to be stacked down. Regular and irregular grids
#' can be considered. For irregular grids or more flexible regular grids, the position of the
#' experimental units (x and y) needs to be provided.
#'
#' @param nblock number of blocks (full replicates)
#' @param ntrt number of  treatments per block
#' @param rb number of rows per block (in the field)
#' @param cb number of columns per block (in the field)
#' @param coord matrix with coordinates (x and y) and blocks for each experimental unit,
#' in the form: Row, Col, Rep
#'
#' @return a data frame with the RCB design with 4 columns: Row, Col, Rep, Treatment
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2016), Generating experimental designs for spatially and
#' genetically correlated data using mixed models, Submitted to Australian and New Zealand
#' Journal of Statistics.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @examples
#' # Example 1: Generates a regular-grid RCB design with 4 blocks and 30 treatments
#' matdf <- rcbd(nblock=4, ntrt=30, rb=5, cb=6)
#' head(matdf)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Example 2: Reads an irregular-grid RCB design with 6 blocks and 30 treatments.
#' head(rcbd_user)
#' matdf <- rcbd(nblock=6, ntrt=30, coord=rcbd_user)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)

rcbd<- function(nblock=0, ntrt=0, rb=0, cb=0, coord=list()) {

  # Checks
  if(nblock==0 || ntrt==0){
    stop("Designs has parameters zero or not declared for nblock and/or ntrt.")
  }

  trt.list <- c(1:ntrt)
  trt.list <- as.numeric(as.factor(trt.list))
  Treatment <- c(replicate(nblock, sample(trt.list,ntrt,replace=FALSE)))

  # Dealing with regular trials
  if (length(coord) == 0){

    if(rb==0 || cb==0){
      stop("Designs has parameters zero or not declared for rb and/or cb")
    }

    Tr <- nblock*rb   # total number of rows in the experiment
    Tc <- cb           # total number of columns in the experiment
    if(nblock == 1){
      Rep <- rep(1:nblock,each=length(unique(trt.list)))
      coord <- cbind(Row=rep(1:Tr,each=Tc),Col=rep(1:Tc,Tr))
      matdf <- cbind(coord,Rep=c(rep(1,ntrt)),Treatment)
      row.names(matdf) <- NULL
    }
    if(nblock > 1){
      Rep <- rep(1:nblock,each=length(unique(trt.list)))
      coord <- cbind(Row=rep(1:Tr,each=cb),Col=rep(1:Tc,rb))
      matdf <- cbind(coord,Rep,Treatment)
      row.names(matdf) <- NULL
    }
  } else {
    coord <- coord[order(coord[,3]),]  # Ordering by Reps
    if(nblock == 1){
      nrep <- nrow(coord)/ntrt
      Treatment <- c(replicate(nrep, sample(trt.list,ntrt,replace=FALSE)))
    }
    matdf <- cbind(coord,Treatment)
    row.names(matdf)<-NULL
  }

  matdf[order(matdf[,"Row"],matdf[,"Col"]),]
  matdf<-data.frame(matdf)
  return(matdf)
}

# Comments
# - Some traps and checks, particularly for irregular designs
