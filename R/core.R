#' Core module for OptimalDesignMM
#'
#' \code{core} Functions that manages the generation and optimization of designs.
#'
#' @param nblock number of blocks (full replicates)
#' @param ntrt number of  treatments per block
#' @param rb number of rows per block (in the field)
#' @param cb number of columns per block (in the field)
#' @param coord matrix with coordinates (x and y) and blocks for each experimental unit,
#' in the form: Row, Col, Rep
#'
#' @return a data frame with ...
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
#' # Example 1: rcbd
#' mydesign<-core_module(design='rcbd', nblock=5, ntrt=30, rb=5, cb=6,
#'                       VarE=0.7, rhox=0.9, rhoy=0.9, VarG=0.3, ng=30, 
#'                       criteria="A", p=5000, plot=TRUE)
#' mydesign$ODEp
#' head(mydesign$matdf)
#' graphics.off()
#' desplot(Rep~Col+Row, mydesign$matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#' 
#' # Example 2: ???

core_module <- function(design=NA, nblock=NA, ntrt=NA, rb=NA, cb=NA,
                        regular=TRUE, matdf=NA, 
                        VarE=NA, rhox=NA, rhoy=NA, 
                        VarG=NA, ng=NA, 
                        criteria="A", p=1000,
                        swap.algorithm='pairwise', swap.method="within", pairs=2,
                        plot=FALSE) {

  # Many more to implement
  if(is.na(design)==T){
    stop("Warning - Please provide the experimental design.")
  }        
  if(p==0 | is.na(p)==T){
    print("Warning - Using default number of iteration p=1000.")
  }
   
  if (design=='crd') {
    print('runing crd')  # for later
  }

  if (design=='rcbd') {

    # Conditions to optimize rcbd
    swap.method<-"within"
    pairs<-2
    
    # Generating design and plotting original design
    Oldmatdf <- rcbd(nblock=nblock, ntrt=ntrt, rb=rb, cb=cb)

    # Building Matrices
    Rinv <- Rmatrix(matdf=Oldmatdf, VarE=VarE, rhox=rhox, rhoy=rhoy, regular=regular)
    Ginv <- Gmatrix(ng=ng, VarG=VarG)   
    resD <- VarCov_bFtR(matdf=Oldmatdf, criteria=criteria, Ginv=Ginv, Rinv=Rinv)  # K is not provided
    K<-resD$K
    
    InitScore<-resD$OptimScore
    OldScore<-resD$OptimScore
    for (i in 1:p) {
      newmatdf <- SwapMethods(Oldmatdf,pairs=pairs,swapmethod=swap.method)
      NewScore <- VarCov_bFtR(matdf=newmatdf, criteria=criteria, Ginv=Ginv, Rinv=Rinv, K=K)$OptimScore
      if(NewScore < OldScore){
        OldScore <- NewScore
        Oldmatdf <- newmatdf
        print(OldScore)
      }
    }
    ODEp<-100*(InitScore-OldScore)/InitScore
    
  }

  if (design=='ibd') {
    print('runing ibd')  # for later
  }

  return(list(matdf=Oldmatdf, Score=OldScore, ODEp=ODEp))

}
