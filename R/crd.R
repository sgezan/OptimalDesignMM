#' Generates a Complete Rndomized Design (CRD))
#'
#' \code{crd} Randomly generates a complete randomized design. Regular and irregular grids
#' can be considered. For irregular grids or more flexible regular grids, the position of the
#' experimental units (x and y) needs to be provided. Note that for regular
#' grids blocks are assumed to be stacked down.
#'
#' @param nrep number of replicates per treatment
#' @param ntrt number oftreatments
#' @param rb number of rows in the field)
#' @param cb number of columns in the field)
#' @param coord matrix with coordinates (x and y) and for each experimental unit,
#' in the form: Row, Col
#'
#' @return a data frame with the CR design with 4 columns: Row, Col, Rep, Treatment. Rep will be
#' a single number 1.
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
#' # Example 1: Generates a  regular-grid CRD design with 12 treatments and 3 replicates.
#' matdf <- crd(nrep=3, ntrt=12, rb=6, cb=6)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Example 2: Reads an irregular-grid CRB design with 10 treatments and 6 replicates.
#' head(crd_user)
#' matdf <- crd(nrep=6, ntrt=10, coord=crd_user)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)

crd<- function(nrep=1, ntrt=0, rb=0, cb=0, coord=list()) {

  if (length(coord) == 0){
     N1 <- nrep*ntrt
     N2 <- rb*cb
     if(N1 != N2){
       stop("Designs has problems on its parameters, verify your parameters")
     }
     matdf <- rcbd(nblock=1, ntrt=ntrt, rb=rb, cb=cb)
  } else {
    Rep <- matrix(data=1,nrow=length(coord),ncol=1)
    coord <- cbind(coord,Rep)
    matdf <- rcbd(nblock=1, ntrt=ntrt, coord=coord)
  }

     matdf<-data.frame(matdf)
     return(matdf)
}

# Comments
# - Some traps and checks, particularly for irregular designs
