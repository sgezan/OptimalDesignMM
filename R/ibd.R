#' Generates an Incomplete Block Design (IBD) (resolvable)
#'
#' \code{ibd} Randomly generates an incomplete block design. Note that for regular
#' grids incomplete blocks (and full resolvable blocks) are assumed to be stacked down.
#' Regular and irregular grids can be considered. For irregular grids or more flexible
#' regular grids, the position of the experimental units (x and y) needs to be provided.
#'
#' @param nblock number of full resolvable blocks (also number of replicates per treatment)
#' @param ntrt number of  treatments per block
#' @param rb number of rows per incomplete block
#' @param cb number of columns per incomplete block
#' @param coord matrix with coordinates (x and y) and blocks for each experimental unit,
#' in the form: Row, Col, FullRep, Rep
#'
#' @return a data frame with the IB design with 5 columns: Row, Col, FullRep, Rep, Treatment
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
#' # Example 1: Generates a regular-grid IB design with 2 full blocks and 36 treatments
#' with incomplete block size of 9 (3*3).
#' matdf <- ibd(nblock=2, ntrt=36, rb=3, cb=3)
#' head(matdf)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Example 2: Reads an user-provided IB design with 2 full blocks and 90 treatments
#' with incomplete block size of 18.
#' head(ibd_user)
#' matdf <- ibd(nblock=2, ntrt=36, rb=3, cb=3, coord=ibd_user)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, out1=FullRep, show.key=FALSE, main=NULL)

ibd <- function(nblock=0, ntrt=0, rb=0, cb=0, coord=list()) {

  nunits <- rb*cb  # Also k dimension of incomplete block
  nincblock <- ntrt/nunits  # number of incomplete blocks per full replicate
  test <- ntrt
  nrep.test <- nblock
  test.names <- c(1:test)
  test.names <- sample(test.names,size=test,replace=FALSE) # Randomize test.names
  test.plot <- matrix(rep(test.names,nrep.test),ncol=nblock,byrow=FALSE)
  Trts <- apply(test.plot, 2, sample)
  Treatment <- matrix(Trts,ncol=1,byrow = FALSE)

  # Dealing with regular trials
  if (length(coord) == 0){

     FullRep <- rep(1:nblock,each=length(Treatment)/nblock)
     Iblock <- rep(1:(nblock*nincblock),each=nunits)
     check <- cbind(FullRep,Iblock,Treatment)
     r <- rep(1:rb,each=cb)
     Row <- r
     for (i in 1:(nblock*nincblock-1)) {
       r <- r+rb
       Row <- c(Row,r)
     }
     Col <- c(rep(c(rep(1:cb,rb)),nblock*nincblock))  # this is fine
     matdf <- data.frame(cbind(Row,Col,FullRep,Iblock,Treatment))
  } else {
     # Dealing with user-provided grids
     #coord <- coord
    matdf <- data.frame(cbind(coord,Treatment))
  }

  colnames(matdf)<-c("Row","Col","FullRep","Rep", "Treatment")
  matdf[order(matdf[,"Row"],matdf[,"Col"]),]
  row.names(matdf)<-NULL

  return(matdf)
}

# Comments
# - Some traps and checks, particularly for irregular designs
