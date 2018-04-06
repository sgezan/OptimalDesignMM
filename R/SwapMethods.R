#' Randomly swap pairs of treatments for optimization
#'
#' \code{SwapMethods} Randomly selects a block and a single or several pairs of treatments, swaps them
#' and creates a new experimental layout.
#'
#' @param matdf an experimental design (layout) based on a randomized complete block designs (RCBD),
#' or incomplete block designs (IBD).
#' @param pairs number of pairs of treatments to be swapped in the run.
#' @param swapmethod is the selected method to be used. This can be "within" which the default
#' (that should be used for randomized complete block designs) or "across" or "any" which can be used for
#' incomplete blocks or unbalanced designs.
#'
#' @return an new experimental design with swapped pairs of treatments with the 4 columns:
#' Row, Col, Rep, Treatment
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2015), Optimal Randomized Complete Block Designs for Spatially
#' and Genetically Correlated Data using Mixed Models, Journal of Agricultural, Biological and
#' Environmental Statistics, 150, 1-32.
#'
#' @author
#' Lazarus Mramba & Salvador Gezan
#'
#' @examples
#' # Example: Regular-grid experiment with independent random effects
#' matdf <- rcbd(nblock=2, ntrt=9, rb=3, cb=3)  # Original design
#' desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping within a single pair
#' newmatdf <- SwapMethods(matdf=matdf, pairs=1, swapmethod="within")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping within 2 pairs
#' newmatdf <- SwapMethods(matdf=matdf, pairs=2, swapmethod="within")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping across a single pair
#' newmatdf <- SwapMethods(matdf=matdf, pairs=1, swapmethod="across")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' # Swapping any (random for within or across) for a single pair
#' newmatdf <- SwapMethods(matdf=matdf, pairs=1, swapmethod="any")
#' desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
#'
#' @export
#' @seealso \code{\link{rcbd}}

SwapMethods <- function(matdf,pairs=1,swapmethod="within") {

  trt <- max(matdf[,"Treatment"])
  gsize <- pairs*2
  stopifnot(gsize %% 2 == 0)
  if(gsize > trt){
    stop("Number of swaps is larger than the number of Treatments")
  }

  # Swaping  any blocks any pairs
  if(swapmethod=="any"){
    list <- c('within','across')
    sel <- sample(list,1)
    swapmethod <- sel
  }

  # Swapping only within a block for any pair
  if(swapmethod=="within"){
    mat <- as.data.frame(matdf, row.names = NULL)
    b1 <- sample(mat$Rep, 1, replace = TRUE)  # Selects a block at random

    g1<- sample(1:length(mat$Rep[mat$Rep==1]),gsize,replace=FALSE)
    gg1 <- mat$Treatment[mat$Rep == b1]
    temp <- gg1

    for (i in seq(1,length(g1),2)) {
      temp[g1[i+1]] <- gg1[g1[i]]
      temp[g1[i]] <- gg1[g1[i+1]]
    }
    mat$Treatment[mat$Rep == b1] <- temp
  }

  # Swaping  across blocks for any pairs
  if(swapmethod=="across"){
    mat <- as.data.frame(matdf, row.names = NULL)
    blks <- sample(unique(mat$Rep),2, replace = FALSE)

    gg1 <- mat$Treatment[mat$Rep == blks[1]]
    gg2 <- mat$Treatment[mat$Rep == blks[2]]
    temp1 <- gg1
    temp2 <- gg2

    g1<- sample(1:length(mat$Rep[mat$Rep==blks[1]]),gsize/2,replace=FALSE)
    g2<- sample(1:length(mat$Rep[mat$Rep==blks[2]]),gsize/2,replace=FALSE)

    for (i in seq(1,length(g1),1)) {
      temp2[g2[i]] <- gg1[g1[i]]
      temp1[g1[i]] <- gg2[g2[i]]
    }
    mat$Treatment[mat$Rep == blks[1]] <- temp1
    mat$Treatment[mat$Rep == blks[2]] <- temp2
  }
  return(mat[order(mat[,"Row"],mat[,"Col"]),])
}
