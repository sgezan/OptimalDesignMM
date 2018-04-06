#' Randomly swap a pair of treatments
#'
#' Randomly selects a block and a pair of treatments, swaps them
#'  and creates a new experimental layout
#'
#' @param matdf an xperimental design
#'
#' @return an experimental design with swapped pairs of treatments
#'
#' @references
#' Mramba, L. K. and Gezan, S. A. (2015), Optimal Randomized Complete Block Designs for Spatially and Genetically Correlated Data using Mixed Models, Journal of Agricultural, Biological and Environmental Statistics, 150, 1-32.
#'
#' @examples
#' blocks = 2; trt = length(LETTERS[1:9]); Tr = 3; Tc = 6; rb = 3; cb = 3
#' matdf <- rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#' newmat <- SwapPair(matdf)
#' which(matdf[,"Treatments"] != newmat[,"Treatments"])
#' DesLayout(matdf, trt, cb, rb, blocks)
#' DesLayout(matdf=newmat, trt, cb, rb, blocks)
#'
#' @export
#' @seealso \code{\link{DesLayout}} for a proper physical layout
#'
#'
SwapPair <- function(matdf) {
  mat <- data.frame(matdf, row.names = NULL)
  b1 <- sample(mat$Reps, 1, replace = TRUE)
  gg1 <- mat$Treatments[mat$Reps == b1]
  g1 <- sample(gg1, 2)
  temp <- gg1
  temp[gg1 == g1[1]] <- g1[2]
  temp[gg1 == g1[2]] <- g1[1]
  mat$Treatments[mat$Reps == b1] <- temp
  mat[order(mat[,"Row"],mat[,"Col"]),]
}
