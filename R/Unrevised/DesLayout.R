#' Print the Experimental Layout
#'
#' \code{DesLayout} prints an experimental layout as an array of blocks.
#'
#' @details the experimental layout generated can display the treatments in string format if the user entered the vector of treatments as strings.
#'
#' @note if the vector of treatments was numerical, the layout will always display the treatments in numerical form
#'
#' @param matdf an xperimental layout (design)
#' @param trt number of treatments per block
#' @param cb number of columns per block
#' @param rb number of rows per block
#' @param blocks number of blocks (replicate)
#'
#' @return returns an experimental layout with blocks arrays
#'
#' @examples
#' blocks = 2; trt = length(LETTERS[1:9]); Tr = 3; Tc = 6; rb = 3; cb = 3
#' matdf <- rcbd(blocks, trt, rb,cb,Tr, Tc,plot=TRUE)
#' DesLayout(matdf, trt, cb, rb, blocks)
#'
#' @export
#' @seealso \code{\link{rcbd}} for generating the \code{matdf}


DesLayout <- function(matdf, trt, cb, rb, blocks) {
  matdf <- matdf[order(matdf[, "Reps"]), ]
  xx = matdf[, "Treatments"]
  trt.list <- c(1:trt)
  yy = trt.list[xx]
  Des_char = aperm(array(yy, dim = c(cb=cb, rb=rb,
                                     blocks=blocks)), perm = c(2, 1, 3))
  Des_nums = aperm(array(matdf[, "Treatments"],
                         dim = c(cb=cb, rb=rb, blocks=blocks)),
                   perm = c(2, 1, 3))
  if(!is.numeric(trt.list)){
     return(list(Layout_char = Des_char, Layout_Nums = Des_nums))
  }else{
    return(Layout_Nums = Des_nums)
  }
}
