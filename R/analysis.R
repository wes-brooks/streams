# attach the necessary libraries
library("SSN")
library('dplyr')
library('sp')
library('rgdal')
library('stringi')
library('Matrix')
library("TMB")
library('MASS')
library('spatstat')

# these scripts define functions that we use later in the analysis.
source('code/kernel.R')
source('code/bisquare.R')
source('code/effort-wrapper.R')
source('code/subset.wrapper.R')
source('code/cox.variational.ICAR.R')
source('code/spatstat-override/linnet.R')
source('code/SSN-override/createDistMatInMemory.R')
source('code/spatstat-override/countends2.R')
source('code/spatstat-override/linearK2.R')
source('code/spatstat-override/linearKengine2.R')
source('code/spatstat-override/linearKinhom2.R')
source('code/spatstat-override/pairdist2.lpp.R')
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
source('code/k-correction/net.R')

# load some stored data products
load('data/intermediate-products/maris-ssn.rdata')
load('data/intermediate-products/stream-network-structure.rdata')
load('data/intermediate-products/eigenvalues.rdata')
load('data/intermediate-products/kernel.wt.pt.rdata')
load('data/intermediate-products/kernel.wt.q.rdata')



# attach the TMB library that is used to estimate the model parameters, compiling it first if necessary.
# compile('src/tmb_ICAR.cpp')
dyn.load(dynlib('src/tmb_ICAR'))

# define some constants
ssnpath <- "data/maris_ssn/maris4.ssn"
spp <- "Brook trout"
spp.short <- list("Brown trout"='brown', "Brook trout"='brook')
epsilon <- 1e-6

# the sampling intensity is defined by a bisquare kernel density with bandwidth 2000 meters
# centered at every location where any species was observed
bws <- c(1e3, 2e3, 5e3, 10e3, 20e3, 50e3) # kernel bandwidths in meters
bw = 2000
offset.cox <- log(kernel.wt.pt[,bws==bw])
offset.q <- log(kernel.wt.q[,bws==bw])







########################################################################################
# identify which networks and reaches have observation points.
obs.rids <- unique(maris@obspoints@SSNPoints[[1]]@point.data$rid)
obs.nets <- as.numeric(as.character(unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)))



#########################################################################################
# establish the order of reaches in the GMRF specification:
indx.valid.net <- edges@data$netID %in% obs.nets
gmrf.lookup <- edges@data[indx.valid.net, c('netID', 'rid')]
gmrf.lookup$gmrfID <- 1:nrow(gmrf.lookup)



########################################################################################
# recall which eigenvalues are connected with which networks (eigenvalues are sorted by network, not in the same order as gmrf.lookup)
netlookup <- vector()
for (nn in obs.nets) {
  len <- sum(gmrf.lookup$netID == nn)
  netlookup <- c(netlookup, rep(nn, len))
}


########################################################################################
# kernel weight obspoints
sampling.point.data <- maris@obspoints@SSNPoints[[1]]@point.data %>% group_by(siteID) %>% distinct(.keep_all=TRUE) %>% as.data.frame

# adjust the weighting of individual reaches and quadrature knots because some reaches
# aren't full covered within one bandwidth of an observation point.
# coverage is the proportion of the reach that is covered by nonzero sampling intensity
coverage <- vector()

for (rid in obs.rids) {
  edgelength <- edges$Length[edges$rid == rid]
  
  if (edgelength < bw) {
    covered <- 1
  } else {
    
    frac <- bw / edgelength # what proportion of the rid length is covered by 1 bw?
    pts <- sort(sampling.point.data$ratio[sampling.point.data$rid == rid]) # 
    
    # turn count active while there's an observation point within a bandwidth
    covered <- 0
    start <- 0
    frontier <- 0
    
    for (pt in pts) {
      if (pt - frac > frontier) {
        # There is a gap between this point and the one prior, so add the last chunk to the proportion covered.
        covered <- covered + (frontier - start)
        start <- max(0, pt - frac)
        frontier <- min(1, pt + frac)
      } else {
        # there is no gap between this point and the one prior so just extend the frontier
        frontier <- min(1, pt + frac)
      }
    }
    
    # increase the coverage proportion by the newly discovered length
    covered <- covered + (frontier - start)
  }
  
  # concatenate the last reach's coverage to the vector
  coverage <- c(coverage, covered)
}
indx <- match(obs.rids, edges$rid)
edges@data$factor <- 1
edges@data$factor[indx] <- 1 / coverage



###########################################################################################
# define rid.lookup, which rid to seg
rid <- edges$rid
rid.lookup <- data.frame(rid=rid, seg=1:length(rid), netID=as.numeric(as.character(edges@data$netID)))
colnames(rid.lookup) <- c('rid', 'seg')



###########################################################################################
# Create a data frame for point data (X.cox) and one for quadrature data (X.q)
# get data from the observation locations and the quarature knots
X.cox <- maris@obspoints@SSNPoints[[1]]@point.data
X.cox$offset <- offset.cox
X.cox$tp <- 1 - X.cox$ratio
X.q <- maris@predpoints@SSNPoints[[1]]@point.data
X.q$offset <- offset.q
X.q$tp <- 1 - X.q$ratio

# convert netword IDs from factors to integers
X.cox$netID <- as.numeric(as.character(X.cox$netID))
X.q$netID <- as.numeric(as.character(X.q$netID))

# filter out rid 39142
X.q <- X.q[X.q$rid %in% edges@data$rid,]

# only nets where points were observed
X.q <- X.q[X.q$netID %in% obs.nets, ]

##############################################
# subsetting to only the even-numbered nets:
# even.nets <- obs.nets[(obs.nets %% 2) == 0]
# X.q <- X.q[X.q$netID %in% even.nets,]
# X.cox <- X.cox[X.cox$netID %in% even.nets, ]
##############################################

# Drop any networks where the only remaining data is a presence point (no quadrature points drives random effect estimation haywire).
indx.orphan <- sapply(obs.nets, function(nn) !(nn %in% unique(X.q$netID)))
if (sum(indx.orphan) > 0) {
  X.cox <- X.cox[!(X.cox$netID %in% obs.nets[indx.orphan]), ]
  obs.nets <- obs.nets[!indx.orphan]
}

# specify the sequential segments
indx <- match(X.q$rid, rid.lookup$rid)
X.q$seg <- rid.lookup$seg[indx]
indx <- match(X.cox$rid, rid.lookup$rid)
X.cox$seg <- rid.lookup$seg[indx]

# extract just the locations where spp was observed.
X.cox <- X.cox %>% filter(comm_name == spp) %>% distinct(rid, ratio, .keep_all=TRUE)

# filter out quadrature knots on networks with no points
X.q <- X.q[X.q$netID %in% unique(X.cox$netID),]

# use only columns that are in both the obs and pred data
common <- names(X.cox)[which(sapply(names(X.cox), function(z) z %in% names(X.q)))]
XX <- rbind(X.cox[,common], X.q[,common])
# XX$netID <- as.integer(levels(XX$netID)[XX$netID])
loc <- data.frame(netID=XX$netID, rid=XX$rid, upDist=XX$upDist)

# -9999 and -9998 are how NAs are coded. Convert these back to NAs.
X.cox[X.cox == -9999 | X.cox == -9998 | X.cox == Inf | X.cox == -Inf] <- NA
X.q[X.q == -9999 | X.q == -9998 | X.q == Inf | X.q == -Inf] <- NA
XX[XX == -9999 | XX == -9998 | XX == Inf | XX == -Inf] <- NA

# Use only the complete cases
X.cox <- X.cox[complete.cases(X.cox),]
X.q <- X.q[complete.cases(X.q),]
XX <- XX[complete.cases(XX),]

# indicator of whether each line segment has a complete data case:
indx <- match(net$lines$rid, XX$rid)
net$lines$known <- !is.na(indx)

# count the observations and quadrature knots
n.q <- nrow(X.q)
n.cox <- nrow(X.cox)
n <- n.cox + n.q

# count the remaining quadrature knots on each remaining reach
seg.q <- data.frame(seg = X.q$seg)
seg.q <- seg.q %>% group_by(seg) %>% summarize(count = n())
count <- rep(0, nrow(rid.lookup))
count[seg.q$seg] <- seg.q$count
rid.lookup$count <- count

# calculate the quadrature weights: Length of a reach (edges$Length),
# multiplied by the proportion with nonzero weight (1 / edges$factor),
# divided by the number of quadrature knots on the reach (rid.lookup$count)
match.indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$wt <- edges@data$Length[match.indx] / edges@data$factor / rid.lookup$count

# this is the key that matches rid to seg
indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$netID <- edges@data$netID[indx]
rid.lookup$Length <- edges@data$Length[indx]
rid.lookup$upDist <- edges@data$upDist[indx]

# match rid to seg on the linear network
indx <- match(net$lines$rid, rid.lookup$rid)
net$lines$seg <- rid.lookup$seg[indx]
net$lines$Length <- rid.lookup$Length[indx]

# indicator of whether each line segment has a complete data case:
indx <- match(net$lines$rid, XX$rid)
net$lines$known <- !is.na(indx)

# match the quadrature weights to the corresponding quadrature knots
match.indx <- match(X.q$rid, rid.lookup$rid)
X.q$wt <- rid.lookup$wt[match.indx]

# We want to log the streamflow, so can't have zeroes.
# Any zero streamflow measurements are set to the minimum nonzero streamflow.
XX$MS_Hist[XX$MS_Hist <= 0] <- min(XX$MS_Hist[XX$MS_Hist > 0])
X.cox$MS_Hist[X.cox$MS_Hist <= 0] <- min(XX$MS_Hist[XX$MS_Hist > 0])
XX$logMS_Hist <- log(XX$MS_Hist)
X.cox$logMS_Hist <- log(X.cox$MS_Hist)

# set the vectors presence and quadrature weights for DWPR
presence <- c(rep(1, n.cox), rep(0, n.q))
p.wt <- c(rep(epsilon, n.cox), X.q$wt)

# create a data frame for the model
dat <- data.frame(resp=presence/p.wt, wt=p.wt)
dat <- within(dat, {
  logMS_Hist <- poly(XX$logMS_Hist, 1);
  S1_93_11 <- poly(XX$S1_93_11, 2);
  SLOPE <- poly(XX$SLOPE, 1);
  offset <- poly(XX$offset, 1);
})
rownames(dat) <- 1:n
dat$id <- rownames(dat)

# estimate the DWPR model
# dwpr <- glm(resp ~ logMS_Hist + S1_93_11 + SLOPE + offset, family=poisson(), weights=wt, data=dat)
# dwpr$loglik <- sum(log(dwpr$fitted.values[1:n.cox])) - n.cox
# dwpr$pres <- dat$pres

# use the centered and scaled coefficients for the Cox model.
X <- model.matrix(~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1), data=XX)
formula.dwpr <- resp ~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1)




###############################################################################################
# estimate the IPP model (no random effect) via the spatstat package
# parts of the resulting object are needed later when calculating the inhomogeneous K function of the Cox model.
obj <- list()
obj$window <- net$window
obj$wt <- p.wt 
obj$loc.cox <- data.frame(x=X.cox$NEAR_X, y=X.cox$NEAR_Y)
obj$loc.q <- data.frame(x=X.q$NEAR_X, y=X.q$NEAR_Y)
obj$coxdata <- head(XX, n.cox)
obj$quaddata <- tail(XX, n.q)
obj$data <- XX
obj$covar.names <- c('logMS_Hist', 'S1_93_11', 'SLOPE', 'offset')
obj$formula <- formula.dwpr

# define the quadrature object as used for DWPR
qq <- quad.dwpr(obj)
qq$covars$offset <- kernel.wrap(kern, rid.lookup, sampling.point.data, edges, bw, nzw.rids=obs.rids, logresult=TRUE)

# add the segment ID as used by spatstat functions to the Brown trout object
cat("About to create the spatstat-type IPP model. This will take a few minutes.\n")
bt.lpp <- lpp(X.cox, net)
obj$dwpr <- lppm.ssnpp(bt.lpp, qq)
cat("Finished!\n")




###############################################################################
# now we can estimate the Cox process model
# create a matrix of spatial random effects
Z <- Matrix(0, nrow=nrow(XX), ncol=nrow(gmrf.lookup))

# populate the random effects design matrix
rid.indx <- match(XX$rid, gmrf.lookup$rid)
colID <- gmrf.lookup$gmrfID[rid.indx]
vecID <- nrow(Z) * (colID-1) + (1:nrow(Z))
Z[vecID] <- 1

# we're only going to fit a random effect to the subset of the networks that have at least one observation point (of species spp)
# identify which elements of the gmrf satisfy the criterion:
indx.live <- gmrf.lookup$netID %in% XX$netID

# subset the random effect design matrix, the random effect precision matrix and its eigenvalues:
Z <- Z[,indx.live]
Sigma.live <- Sigma.mask[indx.live, indx.live]
eigenv.live <- eigenvalues[netlookup %in% unique(XX$netID)]
sindx <- apply(Z, 1, function(r) ifelse(sum(r)==1, which(r==1)-1, -1))

# generate a matrix that, when multiplied by the vector of random effect coefficients, will recover
# the mean of the random effect on each stream network:
A.lookup <- unique(XX$netID)
A.indx <- match(gmrf.lookup$netID, A.lookup)
A.indx <- A.indx[!is.na(A.indx)]
netmeans.matrix <- Matrix(0, nrow=length(A.lookup), ncol=ncol(Z))
rowID <- A.indx
netmeans.lookup <- as.data.frame(rowID) %>% group_by(rowID) %>% summarize(ct=n())
vecID <- rowID + ((1:ncol(netmeans.matrix)) - 1) * nrow(netmeans.matrix)
netmeans <- 1 / netmeans.lookup$ct[rowID]
netmeans.matrix[vecID] <- netmeans

# set starting values for the estimation algorithm.
logV <- rep(0.01, ncol(Z))
M <- rep(0.01, ncol(Z))
resp <- presence / p.wt
beta <- coef(obj$dwpr)
sindx <- apply(Z, 1, function(r) ifelse(sum(r)==1, which(r==1)-1, -1))

# precision matrix is defined intrinsically so there are some eigenvalues equal to zero
# for calculating the log-determinant, ignore these eigenvalues
logDetPrecision <- sum(log(eigenv.live[eigenv.live > sqrt(.Machine$double.eps)]))

# estimate the Cox model
cat("About to estimate the Cox process model parameters. This will take several  minutes.\n")
cox.time <- system.time(cox <-
                          cox.variational.ICAR(resp, X, Z,  Sigma.live, sindx, p.wt, beta.start = beta,
                                               means.matrix = netmeans.matrix, u.start = M, ltau.start=0.01,
                                               logDetQ = logDetPrecision,
                                               logV.start = logV, hess=FALSE, sd=FALSE,
                                               tol=.Machine$double.eps, maxit=1e7))
cox$time <- cox.time
obj$cox <- cox
cat("Finished!\n")




#############################################################################################
# generate the matrices that can back-transform the Cox process coefficients and covariance to the scale of the raw covariates
# Generate the matrix that converts the coefficients from orthogonal polynomials to raw polynomials
lincomb <- diag(ncol(X))
quad.coefs <- c('S1_93_11')
lin.coefs <- c('SLOPE', 'logMS_Hist', 'offset')
for (cc in quad.coefs) {
  locname <- paste("poly(", cc, ", 2)", sep='')
  lindx <- which(names(coef(obj$dwpr)) == paste(locname, '1', sep=''))
  quadx <- which(names(coef(obj$dwpr)) == paste(locname, '2', sep=''))
  
  pars <- attr(dat[[cc]], 'coefs')
  
  # Orthogonal polynomials are a complicated transformation.
  # Add some terms to the intercept
  lincomb[1, lindx] <- -pars$alpha[1] / sqrt(pars$norm2[3])
  lincomb[1, quadx] <- -mean(XX[[cc]]^2) / sqrt(pars$norm2[4]) +
    mean(XX[[cc]]) * sum(dat[[cc]][,1]*XX[[cc]]^2) / sqrt(pars$norm2[4]) / sqrt(pars$norm2[3])
  
  # Raw linear coefficient is a linear combination of the orthogonal polynomial linear and quadratic coefficients
  lincomb[lindx, lindx] <- 1 / sqrt(pars$norm2[3])
  lincomb[lindx, quadx] <- -sum(dat[[cc]][,1]*XX[[cc]]^2) / sqrt(pars$norm2[4]) / sqrt(pars$norm2[3])
  
  # Raw quadratic coefficient is just scaled from the orthogonal polynomial quadratic coefficient
  lincomb[quadx, quadx] <- 1 / sqrt(pars$norm2[4])
  
  # lincomb[1, quadx] <- -(mean(XX[[cc]]) * lincomb[lindx, quadx] + mean(XX[[cc]]^2) * lincomb[quadx, quadx])
}
# back-tranfsorm the coefficients that were only centered and scaled:
for (cc in lin.coefs) {
  locname <- paste("poly(", cc, ", 1)", sep='')
  lindx <- which(names(coef(obj$dwpr)) == locname)
  
  pars <- attr(dat[[cc]], 'coefs')
  
  # Orthogonal polynomials are a complicated transformation.
  # Add some terms to the intercept
  lincomb[1, lindx] <- -pars$alpha[1] / sqrt(pars$norm2[3])
  
  # Raw linear coefficient is a linear combination of the orthogonal polynomial linear and quadratic coefficients
  lincomb[lindx, lindx] <- 1 / sqrt(pars$norm2[3])
}


#################################################################################
# Back-transform the regression coefficients back to the scale of the raw data for interpretation:
# constrain optimization to only the fixed effects and calculate the Hessian
paramselect <- c(TRUE, rep(FALSE, length(obj$cox$M)), rep(FALSE, length(obj$cox$M)), rep(TRUE, length(beta)))
fn <- subset.wrapper(obj$cox$obj$fn, paramselect, obj$cox$res$par)
gr <- subset.wrapper(obj$cox$obj$gr, paramselect, obj$cox$res$par, vector.result=TRUE)
constrained <- nlminb(c(obj$cox$ltau, obj$cox$beta), objective=fn, gradient=gr, control=list(iter.max=1e7, eval.max=1e7, rel.tol=.Machine$double.eps))
constrained.hess <- Matrix(optimHess(constrained$par, fn=fn, gr=gr, control=list(reltol=.Machine$double.eps)))

# use linear combination matrix to transform the coefficients and covariances
beta <- as.vector(lincomb %*% obj$cox$beta)
covmat <- lincomb %*% solve(constrained.hess[2:7,2:7]) %*% t(lincomb)

# report the back-transformed regression coefficients and their standard errors
coefficients <- data.frame(beta=beta, sd=sqrt(tail(diag(covmat), 6)))
print("Estimates of the regression parameters: ")
print(cbind(colnames(X), round(coefficients, 3)))

# report the back-transpormed variance component for the random field
ltau.samples <- rnorm(1000, mean=obj$cox$ltau, sd = 1 / sqrt(diag(constrained.hess)[1]))
cat("mean of the estimated variance component: ", mean(1/exp(ltau.samples) ), '\n')
cat("standard deviation of the estimated variance component: ", sd(1/exp(ltau.samples)), '\n')









#########################################################################################
# recover the intensity estimates for use in the inhomogeneous K function and parametric bootstrap simulations:
std <- exp(obj$cox$logV / 2)
u.realized <- rnorm(length(std), obj$cox$M, std) - std^2 / 2  #a realization of the random effects, corrected for log transformation bias

# calculate the intensity at all observation & quadrature locations in order to get the (approximate) max intensity, lmax:
fixed <- as.vector(X %*% cox$beta)
ranef <- as.vector(Z %*% u.realized)
lambda <- exp(fixed + ranef)



##############################################################################
# compute the inhomogeneous K function
indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$netID <- edges@data$netID[indx]
rid.lookup$Length <- edges@data$Length[indx]
rid.lookup$upDist <- edges@data$upDist[indx]
mylpp <- lpp(X.cox, net)
cat("About to calculate the inhomogeneous K function. This will take several minutes (~20ish).\n")
inhom <- linearKinhom2(mylpp, net, lambda, lookup=rid.lookup, ssn=maris, edges=edges@data)

obj$cox$inhomogeneousK <- inhom
modelname <- paste('obj', spp.short[[spp]], sep='.')
assign(modelname, obj)

save(list=c('X', 'XX', 'X.cox', 'X.q', 'n.cox', 'n.q', 'n', 'Sigma.mask', 'netmeans', 'lincomb', modelname),
     file=paste('output/', spp.short[[spp]], '-cox/', spp.short[[spp]], '.model.rdata', sep=''))

stop("Stopping before simulating - if you want to try simulation, run the last section of the analysis.R script manually.")





###############################################################################
# simulate a new data set:
# hold out effort's effect because we need to calculate it at each simulated location:
noeffort <- as.vector(X[,1:5] %*% cox$beta[1:5])
lambda.noeffort <- exp(noeffort + ranef - mean(XX$offset) * beta[6])

# generate a lookup function for the estimated intensity, holding out effort.
lamdat <- data.frame(x=XX$NEAR_X, y=XX$NEAR_Y, seg=XX$seg, tp=XX$seg, lambda=lambda.noeffort)
lam <- linfun.factory(lamdat, 'lambda')

# wrap up a fuction that returns the effort at any point on the river network
lmax <- max(lambda)
kf <- kernel.wrap(kern, rid.lookup, ssn.kern@obspoints@SSNPoints[[1]]@point.data, edges, bw, nzw.rids=obs.rids, distpath='~/Dropbox/streams/data/maris_ssn/maris4.ssn/distance/obs/')
eff <- effort.wrapper.kernel(edges@data, rid.lookup, beta[6], kf)

# run the simulation script:
cat("About to simulate a point process This will take about half an hour.\n")
sim <- simulate2.lppm(obj$dwpr, lambda=lam, lmax=lmax, effort=eff, nzw.rids=obs.rids, loglambda=FALSE)