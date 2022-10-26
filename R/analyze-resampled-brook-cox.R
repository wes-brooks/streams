# use DWPR to analyze a resampled linear point process
library("SSN")
library('dplyr')
library('sp')
library('rgdal')
library('stringi')
library('Matrix')
library('MASS')
library('spatstat')
library('TMB')


dyn.load(dynlib("tmb_ICAR"))
# load("prep-brook-trout.RData")


# get the id of this simulation
id <- Sys.getenv("id")


ssnpath <- "data/maris_ssn/maris4.ssn"
spp <- "Brook trout"
epsilon <- 1e-10





########################################################
# import SSN data and the supplement
ssnpath <- "data/maris_ssn/maris4.ssn"
maris <- importSSN(ssnpath, predpts='preds')
supplement <- read.csv("data/quadpts.csv")
supplement$OBSPRED_ID <- 1:nrow(supplement) + max(maris@predpoints@SSNPoints[[1]]@point.data$OBSPRED_ID)
colnames(supplement)[23:24] <- c('NEAR_X', 'NEAR_Y')
supplement$NEAR_ANGLE <- 0
supplement$NEAR_DIST <- 0
supplement$pid <- 1:nrow(supplement) + max(maris@predpoints@SSNPoints[[1]]@point.data$pid)
supplement$locID <- supplement$siteID <- as.factor(1:nrow(supplement) + max(as.integer(levels(maris@predpoints@SSNPoints[[1]]@point.data$locID))))
colnames(supplement)[colnames(supplement) == 'ELEV_NSI'] <- 'Elev_NSI'
supplement$ratio <- 0.5
supplement$NEAR_FID <- -1
supplement$netID <- as.factor(supplement$netID)

supplement$FID <- NULL
supplement$ORIG_FID <- NULL
supplement$quad <- NULL
supplement$COMID <- NULL
supplement$GNIS_NAME <- NULL
supplement$FTYPE <- NULL
supplement$areaPI <- NULL
supplement$Length <-  NULL

# join the supplement and the original maris data
maris@predpoints@SSNPoints[[1]]@point.data <- rbind(maris@predpoints@SSNPoints[[1]]@point.data, supplement)
maris@predpoints@SSNPoints[[1]]@point.coords <- rbind(maris@predpoints@SSNPoints[[1]]@point.coords, as.matrix(supplement[,c('NEAR_X', 'NEAR_Y')]))
colnames(maris@predpoints@SSNPoints[[1]]@point.coords) <- c('coords.x1', 'coords.x2')
npc <- data.frame(DistanceUpstream = c(maris@predpoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream, supplement$upDist))
npc$NetworkID <- as.factor(c(as.character(maris@predpoints@SSNPoints[[1]]@network.point.coords$NetworkID), as.character(supplement$netID)))
npc$SegmentID <- as.factor(c(as.character(maris@predpoints@SSNPoints[[1]]@network.point.coords$SegmentID), as.character(supplement$rid)))
maris@predpoints@SSNPoints[[1]]@network.point.coords <- npc[, c('NetworkID', 'SegmentID', 'DistanceUpstream')]







################################################
# import the data and the resampled points
load("data/brook-etc.rdata")
load(paste("data/cox-resamples-brook/points-brook-cox", id, "RData", sep='.'))
outfile <- "output/brook-cox/brook-cox-results"

indx <- match(sim[[1]]$data$seg, rid.lookup$seg)
sim[[1]]$data$rid <- rid.lookup$rid[indx]
sim[[1]]$data$netID <- rid.lookup$netID[indx]





############################################
# source the scripts that overrive built-in spatstat functions
source('code/kernel.R')
source('code/bisquare.R')
source('code/effort-wrapper.R')
source('code/cox.variational.ICAR.R')
source('code/spatstat-override/linnet.R')
source('code/SSN-override/createDistMatInMemory.R')
source('code/spatstat-override/countends2.R')
source('code/spatstat-override/linearK2.R')
source('code/spatstat-override/linearKengine2.R')
source('code/spatstat-override/linearKinhom2.R')
source('code/spatstat-override/project2segment2.R')
source('code/spatstat-override/pairdist2.lpp.R')
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
source('code/k-correction/net.R')




# kill any presence points that are on a rid with no quadrature points
obs.rids <- unique(maris@obspoints@SSNPoints[[1]]@point.data$rid)
rid.noquad <- obs.rids[!(obs.rids %in% maris@predpoints@SSNPoints[[1]]@point.data$rid)]
indx.noquad <-  maris@obspoints@SSNPoints[[1]]@point.data$rid %in% rid.noquad
maris@obspoints@SSNPoints[[1]]@point.data <- maris@obspoints@SSNPoints[[1]]@point.data[!indx.noquad,]
obs.rids <- unique(maris@obspoints@SSNPoints[[1]]@point.data$rid)

# New GMRF precision matrix based on individual stream segments:
# First, let's consider only rids on networks with at least one point.
obs.nets <- as.numeric(as.character(unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)))
indx.valid.net <- edges@data$netID %in% obs.nets
gmrf.lookup <- edges@data[indx.valid.net, c('netID', 'rid')]
gmrf.lookup$gmrfID <- 1:nrow(gmrf.lookup)





################################################
# get the order of the eigenvalues
netlookup <- vector()
for (nn in obs.nets) {
  len <- sum(gmrf.lookup$netID == nn)
  netlookup <- c(netlookup, rep(nn, len))
}



##################################################
# reconstitute the X.cox data.frame using the simulated data
X.cox <- sim[[1]]$data

# match.indx <- match(X.cox$rid, basis$rid)
# X.cox$br <- basis$br[match.indx]
# X.cox$br[is.na(X.cox$br)] <- -1
colnames(X.cox)[colnames(X.cox)=='x'] <- "NEAR_X"
colnames(X.cox)[colnames(X.cox)=='y'] <- "NEAR_Y"
X.cox$ratio <- 1 - X.cox$tp

# filter out quadrature knots on networks with no points
orig.nets <- unique(X.q$netID)
indx.q <- X.q$netID %in% unique(X.cox$netID)
drop.nets <- unique(X.q$netID[!indx.q])
X.q <- X.q[indx.q,]
wt.q <- wt.q[indx.q]

# create a matrix of spatial random effects
Z.q <- tail(Z, n.q)[indx.q,]
n.q <- nrow(Z.q)

#log the flow
X.q$MS_Hist[X.q$MS_Hist <= 0] <- min(X.q$MS_Hist[X.q$MS_Hist > 0])
X.q$logMS_Hist <- log(X.q$MS_Hist)

# use only columns that are in both the obs and pred data
common <- names(X.cox)[which(sapply(names(X.cox), function(z) z %in% names(X.q)))]
XX <- rbind(X.cox[,common], X.q[,common])

# add some metadata to the data
match.indx <- match(XX$rid, rid.lookup$rid)
XX$seg <- rid.lookup$seg[match.indx]
XX$tp <- 1 - XX$ratio

# count the points
n.cox <- nrow(X.cox)
n <- n.cox + n.q

# quadrature weights
presence <- c(rep(1, n.cox), rep(0, n.q))
p.wt <- c(rep(epsilon, n.cox), wt.q)


################################################
# generate the data for fitting the model
# 
# # create a data frame for the model
dat <- data.frame(resp=presence/p.wt, wt=p.wt, rid=XX$rid, seg=XX$seg, tp=XX$tp, x=XX$NEAR_X, y=XX$NEAR_Y)
dat <- within(dat, {
    logMS_Hist <- poly(XX$logMS_Hist, 1);
    #W95_Hist <- poly(XX$W95_Hist, 1);
    #Elev_NSI <- poly(XX$Elev_NSI, 1);
    S1_93_11 <- poly(XX$S1_93_11, 2);
    SLOPE <- poly(XX$SLOPE, 1);
    #CFM_Hist <- poly(XX$CFM_Hist, 1);
    offset <- poly(XX$offset, 1);
})

dat$logMS_Hist <- XX$logMS_Hist
dat$S1_93_11 <- XX$S1_93_11
dat$SLOPE <- XX$SLOPE
dat$offset <- XX$offset

rownames(dat) <- 1:n
dat$id <- rownames(dat)

X <- model.matrix(~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1), data=XX)
resp <- presence / p.wt


# populate the random effects design matrix. 
Z.cox <- Matrix(0, nrow=nrow(X.cox), ncol=ncol(Z.q))
indx.orig <- gmrf.lookup$netID %in% orig.nets
rid.indx <- match(X.cox$rid, gmrf.lookup$rid[indx.orig])
colID <- rid.indx
vecID <- nrow(Z.cox) * (colID-1) + (1:nrow(Z.cox))
Z.cox[vecID] <- 1
Z <- rbind(Z.cox, Z.q)
indx.z <- gmrf.lookup$netID[indx.orig] %in% XX$netID
Z <- Z[,indx.z]


# subset the big precision matrix to only include the networks that appear in this data.
indx.sigma <- gmrf.lookup$netID %in% XX$netID
Sigma.evens <- Sigma.mask[indx.sigma, indx.sigma]
eigenvalues.evens <- eigenvalues[netlookup %in% unique(XX$netID)]



# the netmeans.matrix is used to impose a Gaussian distribution on the means of each network's random effects
# this is necessary because the random effect has an intrinsic specification, and therefore the mean is otherwise allowed to
# deviate far from zero with no penalty.
sindx <- apply(Z, 1, function(r) ifelse(sum(r)==1, which(r==1)-1, -1))
A.lookup <- unique(XX$netID)
A.indx <- match(gmrf.lookup$netID, A.lookup)
A.indx <- A.indx[!is.na(A.indx)]
netmeans.matrix <- Matrix(0, nrow=length(A.lookup), ncol=ncol(Z))
rowID <- A.indx
netmeans.lookup <- as.data.frame(rowID) %>% group_by(rowID) %>% summarize(ct=n())
vecID <- rowID + ((1:ncol(netmeans.matrix)) - 1) * nrow(netmeans.matrix)
netmeans <- 1 / netmeans.lookup$ct[rowID]
netmeans.matrix[vecID] <- netmeans


wes.time <- system.time(cox <- cox.variational.ICAR(resp, X, Z,
  Sigma.evens, sindx, p.wt, beta.start=vari2$beta,
  means.matrix=netmeans.matrix, u.start=vari2$M[indx.z], ltau.start=vari2$ltau,
  logDetQ=sum(log(eigenvalues.evens[eigenvalues.evens > sqrt(.Machine$double.eps)])),
  logV.start = vari2$logV[indx.z], hess=FALSE, sd=FALSE, tol=.Machine$double.eps,
  maxit=1e7))


# cox <- cox.variational.ICAR(resp, X, Z, Sigma.mask, sindx, logV.start = logV, p.wt, beta, hess=FALSE, sd=FALSE, tol=epsilon)
cox <- cox.variational.ICAR(resp, X, Z, Sigma.evens, sindx, p.wt,
        means.matrix=netmeans.matrix, beta.start=vari2$beta,
        logDetQ=sum(log(eigenvalues.evens[eigenvalues.evens > sqrt(.Machine$double.eps)])),
        u.start=vari2$M[indx.z], ltau.start = vari2$ltau, logV.start =vari2$logV[indx.z],
        hess=FALSE, sd=FALSE, tol=.Machine$double.eps, maxit=1e7)



#################################
# find the maximum value of lambda
# Link the random effects to the branches in the random effect basis
std <- exp(cox$logV / 2)
u.realized <- rnorm(length(std), cox$M, std) - std^2 / 2


fixed <- as.vector(X %*% cox$beta)
ranef <- as.vector(Z %*% u.realized)
lambda <- exp(fixed + ranef)
#####################################




# compute the inhomogeneous K function
indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$netID <- edges@data$netID[indx]
rid.lookup$Length <- edges@data$Length[indx]
rid.lookup$upDist <- edges@data$upDist[indx]
mylpp <- lpp(X.cox, net)
inhom <- linearKinhom2(mylpp, net, lambda, lookup=rid.lookup, ssn=maris, edges=edges@data)



save(inhom, cox, dat, file=paste(outfile, id, "RData", sep='.'))

