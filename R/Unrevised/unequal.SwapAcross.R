#' Randomly Swap A Pairs of Treatments across blocks
#'
#' \code{unequal.SwapAcross} Randomly swaps any pairs of treatments across blocks
#'
#' @param matdf an xperimental design
#'
#' @return an experimental design with swapped pairs of treatments
#'
#' @references
#' Mramba, Lazarus. K. and Gezan, Salvador. A. (2016), Improving Unequally Replicated,
#' Incomplete Block and Augmented Experimental Designs with Spatially
#' and Genetically Correlated Observations, Submitted to the Journal of
#' Theoretical and Applied Genetics
#'
#' @examples
#'  blocks = 2; trt = length(LETTERS[1:9]); Tr = 3; Tc = 6; rb = 3; cb = 3
#'  matdf <- rcbd(blocks, trt,rb,cb, Tr, Tc,plot=TRUE)
#'  newdes <- unequal.SwapAcross(matdf)
#'  which(matdf[,"Treatments"] != newdes[,"Treatments"])
#'
#' @export
#' @seealso \code{\link{SwapPair}}, \code{\link{ReplaceTrts}}
#'
unequal.SwapAcross <- function(matdf) {
  mat <- data.frame(matdf, row.names = NULL)
  blks <- sample(unique(mat$Reps),2, replace = FALSE)
  gg1 <- mat$Treatments[mat$Reps == blks[1]]
  gg2 <- mat$Treatments[mat$Reps == blks[2]]
  g1 <- sample(gg1, 1)
  g2 <- sample(gg2, 1)
  temp1 <- gg1
  temp2 <- gg2
  temp1[gg1 == g1] <- g2
  temp2[gg2 == g2] <- g1
  mat$Treatments[mat$Reps == blks[1]] <- temp1
  mat$Treatments[mat$Reps == blks[2]] <- temp2
  mat[order(mat[,"Row"],mat[,"Col"]),]
}
