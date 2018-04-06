#'  Generates the Numerator Relationship Matrix from pedigree information
#'
#' \code{GenA} uses a sorted pedigree file to generate a numerator relationship matrix (A). All individuals in pedigree need to be defined. If parents are missing then they are specified as 0 or NA.
#' @import Matrix
#' @import nadiv
#' @import MASS
#'
#' @param male a vector of males (sires)
#' @param female a vector of females (dams)
#' @return a sparse numerator relationship matrix of class "dgCMatrix"
#'
#' @references
#' Mramba, L.K. and Gezan, S.A. (2015), Optimal Randomized Complete Block Designs for Spatially and Genetically Correlated Data using Mixed Models, Journal of Agricultural, Biological and Environmental Statistics, 150, 1-32.
#'
#' Mrode, R.A. (2014), Linear Models for the Prediction of Animal Breeding Values, 3rd ed., CABI, Oxfordshire, UK.
#'
#' @author
#' Lazarus Mramba
#'
#' @examples
#' ## Example 1: Generates a matrix from simple pedigree
#' indiv <- c(1,2,3,4,5,6,7,8)
#' sire <- c(0,0,0,0,1,1,2,2);
#' dam <- c(0,0,0,0,3,4,3,4)
#' Amatrix <- GenA(male=sire,female=dam)
#' Amatrix
#'
#' ## Example 2: Pedigree from Mrode (2014)
#' require(nadiv)
#' Amrode <- GenA(male=Mrode2$sire,female=Mrode2$dam);
#' Amrode
#'
#' ## Example 3: Using Full-Sib pedigree
#' data(ped30fs)
#' head(ped30fs)
#' AFS <- GenA(male=ped30fs$male,female=ped30fs$female)
#' head(AFS)
#'
#' @export

GenA <- function(male, female) {
    if (nargs() == 1) {
        stop("require male and female entries")
    }
    if (length(male) != length(female)) {
        stop("length of male and female differ")
    }
    male[is.na(male)] <- 0  # convert all NA to zeros
    female[is.na(female)] <- 0  # convert all NA to zeros
    n <- length(male)
    N <- n + 1
    A <- matrix(0, ncol = N, nrow = N)
    male <- (male == 0) * (N) + male
    female <- (female == 0) * N + female
    for (i in 1:n) {
        A[i, i] <- 1 + A[male[i], female[i]]/2
        for (j in (i + 1):n) {
            if (j > n)
                break
            A[i, j] <- (A[i, male[j]] + A[i, female[j]])/2
            A[j, i] <- A[i, j]
        }
    }
    A <- as(A, "sparseMatrix")
    return(A[1:n, 1:n])
}


# Note:
# - It needs a check that all indivs in male and female are defined (for which it needs indiv)
