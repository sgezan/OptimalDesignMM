#' Generate an Augmented/Unreplicated/P-replicated experimental design
#'
#' \code{augmented} Generates an initial RCB design with replicated control treatments and
#' adds unreplicated/partially replicated test treatments. Note that blocks are assumed to be
#' stacked down. (Warining: function only available for regular designs)
#'
#' @param nblock number of blocks (full replicate)
#' @param test.trt is a list for the test treatments with the first number corresponding to
#' the number of levels and the second to how many replications of each level are desired (defaul = 1).
#' @param check.trt is a list for the check (control) treatments with the first number corresponding to
#' the number of levels and the second to how many replications of each level are desired (defaul = 1).
#' @param rb number of rows per block (in the field)
#' @param cb number of columns per block (in the field)
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
#' # Example 1: Generates a RCB design with 4 blocks and 36 treatments + 3 checks
#' matdf <- augmented(nblock=4, test.trt=c(36,1), check.trt=c(3,1), rb=4, cb=3)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
#'
#' # Example 2: Generates a RCB design with 2 blocks, 16 treatments (double replicated)
#' # and 3 checks (triple replication)
#' matdf <- augmented(nblock=2, test.trt=c(16,2), check.trt=c(3,3), rb=5, cb=5)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
#'
#' # Example 3: Generates a CRB design, 136 treatments + 3 checks (with 8 replications each)
#' matdf <- augmented(nblock=1, test.trt=c(136,1), check.trt=c(3,8), rb=16, cb=10)
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
#'
#' @export
#' @seealso \code{\link{rcbd}}, \code{\link{XXX}}

augmented <- function(nblock=0, test.trt=c(0,1), check.trt=c(0,1), rb=0, cb=0) {

  test <- test.trt[1]
  nrep.test <- test.trt[2]
  control <- check.trt[1]
  nrep.control <- check.trt[2]

  # Checking input
  N1 <- rb*cb
  N2 <- control*nrep.control + test*nrep.test/nblock
  if(N1!=N2){
    stop("Designs has problems on its parameters, verify the dimension of the blocks")
  }

  test.names <- c(1:test)
  test.names <- sample(test.names,size=test,replace=FALSE) # Randomize test.names
  control.names <- 1000+c(1:control)
  test.plot <- matrix(rep(test.names,nrep.test),ncol=nblock,byrow=FALSE)
  chk.plot <- matrix(rep(rep(control.names,nrep.control),nblock),ncol=nblock)
  newmat <- rbind(chk.plot,test.plot)

  Trts <- apply(newmat, 2, sample)
  Treatment <- matrix(Trts,ncol=1,byrow = FALSE)
  #N <- length(Treatment)
  Rep <- rep(1:nblock,each=length(Treatment)/nblock)

  r <- rep(1:rb,each=cb)
  Row <- r
  if (nblock > 1) {
     for (i in 1:(nblock-1)) {
        r <- r + rb   #rb or cb?????
        Row <- c(Row,r)
     }
  }
  Col <- c(rep(c(rep(1:cb,rb)),nblock))
  cord <- cbind(Row=Row,Col=Col)

  matdf <- data.frame(cbind(cord,Rep=Rep,Treatment=Treatment))
  colnames(matdf)<-c("Row","Col","Rep","Treatment")
  matdf[order(matdf[,"Row"],matdf[,"Col"]),]
  row.names(matdf)<-NULL

  return(matdf)
}

# Comments
# - Do we want to have coord here?
