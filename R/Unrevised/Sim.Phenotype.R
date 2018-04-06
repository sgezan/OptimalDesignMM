#' Simulates a phenotype of interest (Y)
#'
#' \code{Sim.Phenotype} Simulates a quantitative response variable (Y) such as height or weight or yield for the experimental design for further analysis
#'
#' @import Matrix
#'
#' @param matdf an experimental layout
#' @param G a variance-covariance matrix
#' @param Tr total number of rows
#' @param Tc total number of columns
#' @param rhox spatial correlation along the rows
#' @param rhoy spatial correlation along the columns
#' @param VarG variance of the treatment effects.
#' @param nugget  un-structured residual error (nugget effect)
#'
#' @details This function simulates a Y variable assumed to be a quantitative trait such as height, weight, volume, yield to be part of the generated experimental design for further analysis of estimation of heritabilities and prediction of genetic values.
#'
#' @references
#' Mramba, L. K. and Gezan, S. A. (2016), Generating experimental designs for spatially and genetically correlated data using mixed models, Submitted to Australian and New Zealand Journal of Statistics.
#'
#' @examples
#' # Example 1 # without pedigree
#' trt = length(c(1:30));criteria="A"
#' blocks = 4; rb=5;cb=6;Tr=10;Tc=12;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#' matdf = rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#'
#'
#'  ## for Independent treatments
#'  m = length(unique(matdf[,"Treatments"]))
#'  G <- VarG * Matrix::Diagonal(m)
#'  ##G <- as(G, "sparseMatrix")
#'
#' Data0 = Sim.Phenotype(matdf,G,Tr,Tc)
#' head(Data0)
#'
#' # Example2 # with pedigree
#' VarG = 0.3; rhox = 0.9; rhoy = 0.9; nugget = 0; criteria="A"
#' blocks = 3; rb = 5; cb = 6; Tr = 15; Tc = 6; trt= length(c(1:30))
#' matdf<- rcbd(blocks,trt,rb,cb, Tr, Tc,plot=TRUE)
#'
#' data("ped30fs")
#' Amat <- GenA(male=ped30fs[,"male"], female = ped30fs[,"female"])
#' Amat <- as.matrix(Amat[-c(1:5), -c(1:5)])
#' G <- VarG * as.matrix(Amat)
#' ##G <- as(G, "sparseMatrix")
#'
#' Data2 <- Sim.Phenotype(matdf, G,Tr,Tc)
#' head(Data2)
#'
#' @export
#' @seealso  \code{\link{rcbd}}, \code{\link{DesLayout}} \code{\link{Error.Surfaces}}
#'

Sim.Phenotype<-function(matdf,G,Tr,Tc,rhox=0,rhoy=0,VarG=0.3,nugget=0)
{
  m = length(unique(matdf[,"Treatments"]))
  G <- as(G, "sparseMatrix")
  L = chol(G)
  q<-matrix(rnorm(m),ncol=1)
  Ghat<-t(L)%*%q
  gg<-cbind(seq(1:m),Ghat)
  colnames(gg)<-c("Treatments","Ghat")
  dat2<-data.matrix(merge(matdf,as.matrix(gg),by="Treatments"))
  dat2<-dat2[order(dat2[,"Row"],dat2[,"Col"]),]
  dat2<-dat2[,c(2,3,4,1,5)]
  zst<-Error.Surfaces(Tr,Tc,rhox,rhoy,VarG,nugget)
  colnames(zst)<-"Esms"
  data.matrix(data.frame(dat2,zst,Y=c(10+dat2[,"Ghat"]+zst),row.names=NULL))
}

