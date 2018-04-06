#' Randomly swap several pairs of treatments
#'
#' Randomly selects a block and a number of pairs of treatments, swaps them
#'  and creates a new experimental layout
#'
#' @param matdf an xperimental design
#' @param gsize number of treatments to be swapped per iteration. This should be an even number
#'
#' @return an experimental design with swapped pairs of treatments
#'
#' @references
#' Mramba, L. K. and Gezan, S. A. (2016), Optimal Randomized Complete Block Designs for Spatially and Genetically Correlated Data using Mixed Models, Submitted to Journal of Theoretical and Applied Genetics.
#'
#' @examples
#' blocks = 2; trt = length(LETTERS[1:9]); Tr = 3; Tc = 6; rb = 3; cb = 3
#' matdf <- rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#' newmat <- SwapGreedy(matdf,gsize=4)
#' which(matdf[,"Treatments"] != newmat[,"Treatments"])
#' DesLayout(matdf, trt, cb, rb, blocks)
#' DesLayout(matdf=newmat, trt, cb, rb, blocks)
#'
#' @export
#' @seealso \code{\link{DesLayout}} for a proper physical layout
#'

SwapGreedy <- function(matdf,gsize) {
  stopifnot(gsize %% 2 == 0)
  mat <- data.frame(matdf, row.names = NULL)
  b1 <- sample(mat$Reps, 1, replace = TRUE)
  gg1 <- mat$Treatments[mat$Reps == b1]
  g1 <- sample(gg1, gsize)
  temp <- gg1
  for (i in seq(1,length(g1),2))
  {
    temp[gg1 == g1[i]] <- g1[i+1]
    temp[gg1 == g1[i+1]] <- g1[i]
  }
  mat$Treatments[mat$Reps == b1] <- temp
  mat[order(mat[,"Row"],mat[,"Col"]),]
}
