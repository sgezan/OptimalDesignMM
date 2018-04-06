#' Swap two specified treatments in a defined block
#'
#' Randomly selects a block and a pair of treatments, swaps them
#'  and creates a new experimental layout
#'
#' @param matdf an xperimental design
#' @param g1 a specific treatment 1
#' @param g2 a specific treatment 2 not the same as treatment 1
#' @param bl which block they belong to
#'
#' @return an experimental design with swapped pairs of treatments
#'
#' @references
#' Mramba, L. K. and Gezan, S. A. (2015), Optimal Randomized Complete Block Designs for Spatially and Genetically Correlated Data using Mixed Models, Journal of Agricultural, Biological and Environmental Statistics, 150, 1-32.
#'
#' @examples
#' blocks = 3; trt = length(c(1:9)); Tr = 3; Tc = 9; rb = 3; cb = 3
#' matdf <- rcbd(blocks,trt,rb,cb,Tr,Tc,plot=TRUE)

#' newmat <- Swap_Specific(matdf,g1=1,g2=4,bl=2)
#' which(matdf[,"Treatments"] != newmat[,"Treatments"])
#' DesLayout(matdf, trt, cb, rb, blocks)
#' DesLayout(matdf=newmat, trt, cb, rb, blocks)
#'
#' @export
#' @seealso \code{\link{DesLayout}} for a proper physical layout
#'
#'
#'
Swap_Specific <- function(matdf,g1,g2,bl) {
  mat <- data.frame(matdf, row.names = NULL)
  gg1 <- mat$Treatments[mat$Reps == bl]
  temp <- gg1
  temp[gg1 == g1] <- g2
  temp[gg1 == g2] <- g1
  mat$Treatments[mat$Reps == bl] <- temp
  mat[order(mat[,"Row"],mat[,"Col"]),]
}
