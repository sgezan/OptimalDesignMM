## TEST SAG EXAMPLE

matdf <- rcbd(nblock=5, ntrt=30, rb=5, cb=6)
head(matdf)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)

Rinv <- Rmatrix(matdf=matdf, VarE=0.7, rhox=0.9, rhoy=0.9, regular=TRUE)
Ginv <- Gmatrix(ng=30, VarG=0.3)   # Independent random effects for a heritability of 0.3
resD <- VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv)  # K is not provided but calculated
resD$OptimScore
K<-resD$K
VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv, K=K)$OptimScore # K is provided

(InitScore<-resD$OptimScore)
(OldScore<-resD$OptimScore)
Oldmatdf<-matdf
print(OldScore)
iter<-2000
for (i in 1:iter) {
 newmatdf <- SwapMethods(matdf,pairs=1,swapmethod="within")
 #newmatdf <- SwapMethods(matdf,pairs=4,swapmethod="across")
 #newmatdf <- SwapMethods(Oldmatdf,pairs=2,swapmethod="any")
 NewScore <- VarCov_bFtR(matdf=newmatdf,criteria="A",Ginv,Rinv,K)$OptimScore # K is provided
 if(NewScore < OldScore){
   OldScore <- NewScore
   Oldmatdf <- newmatdf
   print(OldScore)
 }
}

(Eff<-100*(InitScore-NewScore)/InitScore)
head(newmatdf)
desplot(Rep~Col+Row, newmatdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)


# Designing an IBD

# Example 1: Generates a regular-grid IB design with 2 full blocks and 36 treatments
# with incomplete block size of 9 (3*3).
matdf <- ibd(nblock=2, ntrt=36, rb=3, cb=3)
head(matdf)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)

# Example 2: Reads an user-provided IB design with 2 full blocks and 90 treatments
# with incomplete block size of 18.
head(ibd_user)
matdf <- ibd(nblock=2, ntrt=36, rb=3, cb=3, coord=ibd_user)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, out1=FullRep, show.key=FALSE, main=NULL)


# Designing an unreplicated trial

# Example 1: Generates a RCB design with 4 blocks and 36 treatments + 3 checks
matdf <- augmented(nblock=4, test.trt=c(36,1), check.trt=c(3,1), rb=4, cb=3)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')

# Example 2: Generates a RCB design with 2 blocks, 16 treatments (double replicated)
# and 3 checks (triple replication)
matdf <- augmented(nblock=2, test.trt=c(16,2), check.trt=c(3,3), rb=5, cb=5)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')

# Example 3: Generates a CRB design, 136 treatments + 3 checks (with 8 replications each)
matdf <- augmented(nblock=1, test.trt=c(136,1), check.trt=c(3,8), rb=16, cb=10)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL, shorten='no')
