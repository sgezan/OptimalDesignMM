#' A procedure to replace some treatments with others
#'
#' \code{ReplaceTrts} replaces some treatments using a criteria.
#'
#' @param matdf an experimental design
#' @param trt number of treatments
#' @param max.u the maximum number treatments available
#' @param min.u the minimum number of treatments allowed
#'
#' @return  a list of objects with a new experimental design
#'
#' @references
#' Mramba, Lazarus. K. and Gezan, Salvador. A. (2016), Improving Unequally Replicated,
#' Incomplete Block and Augmented Experimental Designs with Spatially
#' and Genetically Correlated Observations, Submitted to the Journal of
#' Theoretical and Applied Genetics
#'
#' @examples
#' ## Example
#' trt = length(1:30);criteria="A"
#' blocks = 3; rb=5;cb=6;Tr=5;Tc=18;rhox=0.6;rhoy=0.6;VarG=0.3;nugget=0
#' min.u = sample(1:3,trt,replace=TRUE)
#' max.u = sample(3:5,trt,replace=TRUE)
#' des1 = unequal.RBD(trt,blocks,max.u,min.u,rb,cb,Tr,Tc,plot=TRUE)
#' des2 = ReplaceTrts(matdf=des1$matdf,trt,max.u,min.u)
#' which(des1$matdf[,"Treatments"] != des2$matdf[,"Treatments"])
#'
#' @export
#' @seealso \code{\link{unequal.RBD}}, \code{\link{rcbd}}
#'
#'
ReplaceTrts <- function(matdf,trt,max.u,min.u)
{
  trt = c(1:trt)
  freq = as.numeric(table(matdf[,"Treatments"]))
  Dif = max.u - freq
  genot <- as.numeric(as.factor(trt))
  datam = cbind(genot,min.u,max.u,freq,Dif)
  u.extra = datam[freq < max.u,"genot"]
  numT <- c(Dif[Dif != 0])
  extraG =rep(u.extra,numT)
  sam <- sample(genot[freq > min.u],1,replace=TRUE)
  sam <- which(matdf[,"Treatments"]==sam)[1]
  Lg <- sample(extraG[extraG!=sam],1,replace=TRUE)
  Lg <- which(matdf[,"Treatments"]==Lg)[1]
  matdf[,"Treatments"][sam] <- matdf[,"Treatments"][Lg]
  freq = as.numeric(table(matdf[,"Treatments"]))
  datam = cbind(genot,min.u,max.u,freq)
  list(matdf=matdf,datam=datam, sumFreq = sum(freq),changed_index_dataframe = sam)
}
