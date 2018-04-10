
## Equation 16

# correlated observations
# fixed nuisance covariate
# Random nuisance block effects
# Treatments are fixed effects

# y = Wg + Xv + Zb + e

# W is dsign matrix for fixed treatments
# g is a vector of fixed treatment effects
# X is a design matrix for all fixed effects other than treatment effects
# v is a vector of all fixed effects other than treatment
# Z is a design matrix of random block effects
# b is a vector of random block effects
# e is a vector of random error terms

# b ~ N(0,B), e ~ N(0,R), b and e are uncorrelated

VarCov.rcbd.blocksRandom <- function(matdf,Tr, Tc, rhox=0, rhoy=0, h2=0.3, s20=0,
                                     criteria="A",sigBl=FALSE,irregular=FALSE) {
  if(nrow(matdf)==length(unique(matdf[,"Treatments"]))){
    X <- as.matrix(matdf[, "Treatments"])
    colnames(X)<-NULL
  }
  if(nrow(matdf) > length(unique(matdf[,"Treatments"]))){
    X <- Matrix::sparse.model.matrix(~as.factor(matdf[, "Treatments"])-1)
    colnames(X)<-NULL
  }
  
  bb <- length(unique(matdf[,"Reps"]))
  if(is.numeric(sigBl))
  {
    Binv <- (1/sigBl)*Matrix::Diagonal(bb)
  }else{
    sigBl <- 0.2*(1 - h2)
    Binv <- (1/sigBl)*Matrix::Diagonal(bb)
    Binv <- as(Binv, "sparseMatrix")
  }
  
  s2e <- (1 - s20) * (1 - h2 - sigBl)
  
  stopifnot(s2e > 0)
  
  Z<- Matrix::sparse.model.matrix(~as.factor(matdf[, "Reps"]) - 1)
  colnames(Z)<-NULL
  
  # calculating R and its inverse for spatial analysis
  matdf <- matdf[order(matdf[,"Row"],matdf[,"Col"]),]
  if(irregular==TRUE){
    R <- Matrix::Diagonal(nrow(matdf))
    for(i in 1:(nrow(matdf)-1)) {
      x1  <- matdf[,"Col"][i]
      y1  <- matdf[,"Row"][i]
      for (j in (i+1):nrow(matdf)){
        x2 <- matdf[,"Col"][j]
        y2  <- matdf[,"Row"][j]
        R[i,j]<-(rhox^abs(x2 -x1))*(rhoy^abs(y2 -y1))
      }
    }
    R = as.matrix(s2e*R)
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
    R <- as(R, "sparseMatrix")
    Rinv <- chol2inv(chol(R))
    Rinv <- as(Rinv, "sparseMatrix")
  }
  
  if(irregular==FALSE){
    sigx <- Matrix::Diagonal(Tc)
    sigx <- rhox^abs(row(sigx) - col(sigx))
    sigy <- Matrix::Diagonal(Tr)
    sigy <- rhoy^abs(row(sigy) - col(sigy))
    R <- s2e * kronecker(sigy, sigx)
    R <- as(R, "sparseMatrix")
    Rinv <- chol2inv(chol(R))
    Rinv <- as(Rinv, "sparseMatrix")
  }
  
  X =  as.matrix(X)
  Z =  as.matrix(Z)
  Rinv = as.matrix(Rinv)
  
  C11 <- Matrix::crossprod(X, Rinv) %*% X
  C11  <- as(C11, "sparseMatrix")
  
  C22 <- t(Z) %*% Rinv  %*% Z + Binv
  C22inv <- solve(C22)
  
  K <- Rinv %*% Z %*% C22inv %*% Matrix::crossprod(Z, Rinv)
  
  K = Matrix::drop0(K)
  X = as.matrix(X)
  
  temp0 <- C11 - t(X) %*% K %*% X
  M <- solve(temp0)
  M  <- as(M, "sparseMatrix")

  X = Matrix::drop0(X)
  Rinv = Matrix::drop0(Rinv)
  
  if (criteria == "A") {
    return(c(traceI = sum(Matrix::diag(M)), K = K, Rinv = Rinv))
  }
  if (criteria == "D") {
    deTm = Matrix::det(M)
    return(c(doptimI = log(deTm), K = K, Rinv = Rinv))
  }
}


#####
#' ## Example 1: Generates a regular-grid RCB design with 4 blocks and 30 treatments
 library(ggplot2)
library(Matrix)
 blocks <- 4; trt <- 30
 rb <- 5; cb <- 6; Tr <- 10; Tc <- 12
 matdf <- rcbd(blocks,trt,rb,cb,Tr,Tc,plot=TRUE)

 ans = VarCov.rcbd.blocksRandom(matdf,Tr, Tc, rhox=0.6, rhoy=0.6, h2=0.3, s20=0,
                                      criteria="A",sigBl=0.00000001,irregular=FALSE)  
   

 ls(ans)
 ans$traceI
 ###
ans2 =  VarCov.rcbd.fixedAll(matdf, Tr, Tc, rhox=0.6, rhoy=0.6,
                                  h2=0.3, s20=0, criteria="A", irregular=FALSE)

ans2$traceI
