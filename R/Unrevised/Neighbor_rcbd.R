#' Generate a new RCB Design using a genetic nearest neighbor algorithm
#'
#' \code{Neighbor_rcbd} Generates an experimental design using a genetic nearest neighbor algorithm. A numerator relationship matrix is required.
#'
#' @param matdf  an  experimental layout
#' @param Amat a numerator relationship matrix
#'
#' @return a new experimental design after swapping treatments that are genetically correlated and in close proximity.
#'
#' @examples
#' # Example
#' h2 = 0.3; rhox = 0.6; rhoy = 0.6; s20 = 0; criteria="A"
#' blocks = 3; rb = 5; cb = 6; Tr = 5; Tc = 18; trt = length(c(1:30))
#' set.seed(100)
#' matdf<- rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)

#' data(ped30fs)
#' Amat <- GenA(male=ped30fs[,"male"], female = ped30fs[,"female"])
#' Amat <- as.matrix(Amat[-c(1:5), -c(1:5)])

#' newmat <- Neighbor_rcbd(matdf,Amat)
#' which(matdf[,"Treatments"] != newmat[,"Treatments"])
#'
#' @export
#'
#' @seealso \code{\link{rcbd}} for design of initial experiments

Neighbor_rcbd<-function(matdf,Amat)
{
  bl<-sample(matdf[,"Reps"],1)
  temp1<- matdf[matdf[, "Reps"] == bl,]
  rb <-length(unique(temp1[,"Row"]))
  cb <-length(unique(temp1[,"Col"]))
  mat<-matrix(c(temp1[,"Treatments"]),nrow=rb,ncol=cb,byrow=TRUE)
  x<-sample(temp1[,"Treatments"],1)
  m2<-cbind(NA,rbind(NA,mat,NA),NA)
  cord <- expand.grid(Row = 1:rb, Col = 1:cb)
  ret<-c()
  for(i in 1:-1)
    for(j in 1:-1)
      if(i!=0 || j !=0)
        ret<-rbind(ret,m2[cord$Row+i+1+nrow(m2)*(cord$Col+j)])
  neigh<-ret[,which(mat==x)]
  neigh<-c(x,neigh[!is.na(neigh)])  # add the genotype plus its neighbors
  test<-t(combn(neigh,2))
  temp<-cbind(test[,1],test[,2],Amat[test])
  swapsData<-subset(temp, temp[,3] >= 0.25)
  samp1 <- setdiff(as.vector(temp1[,"Treatments"]),c(x,as.vector(temp[,2])))
  if(nrow(swapsData)!=0 && length(samp1) !=0)
  {
    val1<-c()
    val2<-c()
    for(i in 1:nrow(swapsData))
    {
      val1[i]<-swapsData[i,1]
      val2[i]<-sample(samp1,1)
      matdf<-Swap_Specific(matdf,g1=val1[i],g2=val2[i],bl=bl)
    }
    return(matdf[order(matdf[,"Row"],matdf[,"Col"]),])
  }
  if(nrow(swapsData)!=0 && length(samp1)==0)
  {
    val1<-c()
    val2<-c()
    for(i in 1:nrow(swapsData))
    {
      val1[i]=swapsData[i,1]
      val2[i]=swapsData[i,2]
      matdf<-Swap_Specific(matdf,g1=val1[i],g2=val2[i],bl=bl)
    }
    return(matdf[order(matdf[,"Row"],matdf[,"Col"]),])
  }
  if(nrow(swapsData)==0) return(SwapPair(matdf))
}
