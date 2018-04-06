rm(list=ls()) # It removes ALL objects

library(ggplot2)
library(desplot)
library(nadiv)

source('R/rcbd.R')          # I made some changes for crd
source('R/Rmatrix.R')
source('R/Gmatrix.R')
source('R/VarCov_bFtR.R')   # blocks Fixed,  treatments Random - This is good!
source('R/VarCov_bRtR.R')   # blocks Random, treatments Random - This is good!
source('R/VarCov_bFtF.R')   # blocks Fixed,  treatments Fixed  - Some issues!!
source('R/VarCov_bRtF.R')   # blocks Random, treatments Fixed  - Some issues!!
source('R/SwapMethods.R')
source('R/ibd.R')
source('R/augmented.R')     # I made some small changes on this one... for a single block
source('R/crd.R')           # Created based on rcbd but easier to use...

# Evalute VarCov_bRtF with an incomplete block design...
# Need software to genrete the IBD
# Expand also to consider CRD not only rcbd (i.e. nblock=1)

rcbd_user <- read.csv(file = 'R/data/rcbd_user.csv')
ibd_user <- read.csv(file = 'R/data/ibd_user.csv')
crd_user <- read.csv(file = 'R/data/crd_user.csv')


# Name Variables

# nblock number of blocks (full replicates)
# ntrt number of  treatments per block
# rb number of rows per block (in the field)
# cb number of columns per block (in the field)
# coord matrix with coordinates (x and y)
# rhox spatial correlation between experimental units along the rows. Default value is 0.
# rhoy spatial correlation between experimental units along the columns. Default value is 0.
# VarE variance of the residuals. Default value is 1.
# nugget spatial nugget error. Default value is 0.

# Most functions will need see also, and update on reference.
# - The roundings can make some inconsistencies in some data...


###############
# A check of the functions

matdf <- rcbd(nblock=16, ntrt=10, rb=1, cb=10)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, show.key=FALSE, main=NULL)
Rinv <- Rmatrix(matdf=matdf, VarE=1, rhox=0.2, rhoy=0.2, regular=TRUE)
Ginv <- Gmatrix(ng=10, VarG=0.3)   # Independent random effects for a heritability of 0.3

# These should be similar - problems!
VarCov_bRtF(matdf=matdf, criteria="A", s2Bl=6000000, Rinv=Rinv)$OptimScore
VarCov_bRtF(matdf=matdf, criteria="A", s2Bl=0.000001, Rinv=Rinv)$OptimScore
VarCov_bFtF(matdf=matdf, criteria="A", Rinv=Rinv)$OptimScore

# These should be similar - perfect!
VarCov_bFtR(matdf=matdf, criteria="A", Ginv=Ginv, Rinv=Rinv)$OptimScore
VarCov_bRtR(matdf=matdf, criteria="A", s2Bl=6000000, Ginv=Ginv, Rinv=Rinv)$OptimScore
VarCov_bRtR(matdf=matdf, criteria="A", s2Bl=0.000001, Ginv=Ginv, Rinv=Rinv)$OptimScore


###########################################################

# Eaxmple ibd
matdf <- ibd(nblock=2, ntrt=36, rb=3, cb=3, coord=ibd_user)
head(matdf)
desplot(Rep~Col+Row, matdf, text=Treatment, cex=1, out1=FullRep, show.key=FALSE, main=NULL)
Rinv <- Rmatrix(matdf=matdf, VarE=1, rhox=0.2, rhoy=0.2, regular=TRUE)
Ginv <- Gmatrix(ng=36, VarG=0.3)   # Independent random effects for a heritability of 0.3
VarCov_bRtF(matdf=matdf, criteria="A", s2Bl=0.1, Rinv=Rinv)$OptimScore
VarCov_bRtR(matdf=matdf, criteria="A", s2Bl=0.1, Ginv=Ginv, Rinv=Rinv)$OptimScore


