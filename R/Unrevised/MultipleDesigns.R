#' Generate Multiple Initial RCB Designs
#'
#' A choice for the number of initial RCB designs is available
#' in order to generate more than one design and select the best.
#'
#' @details
#' Generates multiple experimental designs, their traces and  determinants
#' where no genetic relationship information exists or where the treatments
#' are known to be genetically unrelated.
#' The best design out of the \code{DesN} initial ones is selected to be
#' optimized based on a choice of optimality criteria.

#' @import Matrix
#' @import ggplot2
#'
#' @param DesN number of initial designs to be generated
#' @param blocks number of blocks (replicate)
#' @param trt number of  treatments per block
#' @param rb number of rows per block
#' @param cb number of columns per block
#' @param Tr total number of rows
#' @param Tc total number of columns
#' @param criteria either ``A" or ``D"
#' @param rhox  spatial correlations along the x-coordinates (rows)
#' @param rhoy  spatial correlations along the y-coordinates (columns)
#' @param nugget  un-structured residual error (nugget effect)
#' @param regular is logical, if FALSE, a different procedure is used to calculate AR1
#' @param plot a logical argument for obtaining plots
#' @param VarG variance of the treatment effects.
#' @param Ginv a variance-covariance matrix from pedigree or molecular data
#' @param Rinv an inverse of the error variance-covariance matrix if available
#' @param K a matrix calculated from original dataset
#' @param Update a logical statement, if TRUE, it updates the criterion values after swaps
#'
#' @return A vector of all traces and determinants from the \code{DesN}
#'  initial experimental designs, the  best experimental layout that
#'  had the smallest trace or log(determinant) depending
#'  on the choice of optimality criteria.
#'
#' @examples
#' # Example
#' VarG = 0.3; rhox = 0.3; rhoy = 0.6; nugget = 0; criteria="A"
#' blocks = 2; rb = 3; cb = 3; Tr = 3; Tc = 6; trt = length(c(1:9))
#'
#' # Treatments are independent
#' Ginv <- round((1/VarG) * Matrix::Diagonal(trt),7)
#' Ginv <- as.matrix(Ginv)
#' Ginv <- as(Ginv, "sparseMatrix")
#'
#' s2e <- 1-VarG-nugget
#' N <-Tr*Tc
#' Rinv <- round((1/(s2e+nugget))*Matrix::Diagonal(N),7)
#'
#' res1 = MultipleDesigns(DesN=4, blocks, trt, rb,cb, Tr, Tc,Ginv,Rinv,
#' criteria="A",plot=TRUE)
#' attributes(res1)
#' DesLayout(matdf=res1$newmatdf, trt, cb, rb, blocks)
#'
#' @export
#'
MultipleDesigns <- function(DesN, blocks, trt, rb,cb, Tr, Tc,Ginv,Rinv,
                            criteria="A", VarG=0.3, rhox=0,rhoy=0,K=FALSE,
                            nugget=0,regular=TRUE,plot=FALSE,Update=FALSE) {
  matrix0 <- list()
  initialValues1 <- c()
  initialValues2 <- c()
  for (i in 1:DesN) {
    print(sprintf("generating initial design: %d", i, "complete\n",
                  sep = ""))
    flush.console()
    matrix0[[i]] <- rcbd(blocks, trt,rb,cb, Tr, Tc, coord=list(), regular)
    initialValues1[i] <- VarCov.rcbd(matdf = matrix0[[i]],
                      Ginv,Rinv,rhox,rhoy,VarG,nugget,
                      criteria,regular,K, Update)[[1]]
    CRITERIA <- ifelse(criteria == "A", "D", "A")
    initialValues2[i] <- VarCov.rcbd(matdf = matrix0[[i]],Ginv,Rinv,
                      rhox,rhoy,VarG,nugget,regular,
                      K, Update,criteria = CRITERIA)[[1]]
    a <- which.min(initialValues1)
    d <- which.min(initialValues2)
    newmatdfA <- matrix0[a][[1]]
    newmatdfD <- matrix0[d][[1]]
    min_initialValues1 <- initialValues1[a][[1]]
    min_initialValues2 <- initialValues2[d][[1]]
  }
  if(plot==TRUE){
    newmatdfA = as.data.frame(newmatdfA)
    P = ggplot2::ggplot(data=newmatdfA)+geom_text(aes(x=newmatdfA$Col,
                                                  y=as.factor(newmatdfA$Row),
                                                  label=newmatdfA$Treatments,
                                                  col=as.factor(newmatdfA$Reps)), size=6) +
      scale_y_discrete("Row coordinates") +
      scale_x_discrete("Column coordinates") +
      scale_color_discrete("Block") +
      ggtitle("Initial design") +
      theme(text = element_text(size=15,face="bold"),
            plot.title = element_text(size=20, face="bold",vjust=2),
            axis.text=element_text(size=17,face="bold"),
            axis.title=element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      coord_fixed()
    print(P)
  }
  return(list(newmatdf = newmatdfA, initialValues1 = initialValues1,
              min_initialValues1 = min_initialValues1, initialValues2 = initialValues2,
              min_initialValues2 = min_initialValues2))
}
