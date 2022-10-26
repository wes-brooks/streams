library("SSN")
library('dplyr')
library('sp')
library('rgdal')
library('stringi')
library('Matrix')
library('MASS')
library('spatstat')
library('TMB')

#dyn.load(dynlib('src/tmb'))
#load("prep-brown-trout.RData")
outfile <- "data/cox-resamples-brown/points-brown-cox"


source('code/cox.variational.R')
source('code/spatstat-override/linnet.R')
source('code/SSN-override/createDistMatInMemory.R')
source('code/spatstat-override/countends2.R')
source('code/spatstat-override/linearK2.R')
source('code/spatstat-override/linearKengine2.R')
source('code/spatstat-override/linearKinhom2.R')
source('code/spatstat-override/project2segment2.R')
source('code/spatstat-override/ppllengine2.R')
source('code/spatstat-override/lppm.ssnpp.R')
source('code/quad.dwpr.R')
source('code/spatstat-override/getlambda2.lpp.R')
source('code/spatstat-override/simulate.lppm.R')
source('code/spatstat-override/rpoislpp2.R')
source('code/spatstat-override/datagen.rpoisppOnLines.R')
source('code/spatstat-override/pointsOnLines2.R')
source('code/spatstat-override/datagen.runifpointOnLines.R')
source('code/spatstat-override/datagen.runifpoisppOnLines.R')
source('code/spatstat-override/runifpoisppOnLines.R')


seed <- as.integer(Sys.getenv("id"))
set.seed(seed)


#load("brown-image.rdata")
load("data/brown-cox.rdata")
#rid.lookup$br[is.na(rid.lookup$br)] <- -1

#lam <- lambda.cox.fun.factory(model.terms, covarfuns, beta, rid.lookup, mu, Sigma)



obj <- list()
Q <- list()

Q$covars <- covarfun
obj$X <- domain
obj$Q <- Q





#################################
# find the maximum value of lambda
# Link the random effects to the branches in the random effect basis
std <- exp(logV / 2)
u.realized <- rnorm(length(std), mu, std) - std^2 / 2  #random effects corrected for log transformation bias

# identify which random effect goes with which entry
# branches with negative id get no random effect
# brange <- rep(NA, max(rid.lookup$br))
# valid <- unique(rid.lookup$br)
# valid <- sort(valid[valid > 0])
# brange[valid] <- 1:length(valid)

# match segment ids to branch ids
# indx <- match(XX$seg, rid.lookup$seg)
# braw <- rid.lookup$br[indx]
# bflag <- braw > 0
# indx <- cbind((1:length(XX$seg))[bflag], brange[braw[bflag]])

# populate the random design matrix
# S <- Matrix(0, nrow=nrow(XX), ncol=length(u.realized))
# S[indx] <- 1

fixed <- as.vector(X %*% beta)
noeffort <- as.vector(X[,1:5] %*% beta[1:5])
ranef <- as.vector(Z %*% u.realized)
lambda <- exp(noeffort + ranef - mean(XX$offset) * eee)
lmax <- max(exp(fixed + ranef))
#####################################


lamdat <- data.frame(x=XX$NEAR_X, y=XX$NEAR_Y, seg=XX$seg, tp=XX$seg, lambda=lambda)
lam <- linfun.factory(lamdat, 'lambda')
sim <- simulate2.lppm(obj, lambda=lam, lmax=lmax, effort=eff, nzw.rids=obs.rids)

save(sim, file=paste(outfile, seed, 'RData', sep='.'))






