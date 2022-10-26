library("SSN")
library('dplyr')
library('sp')
library('rgdal')
library('stringi')
library('Matrix')
library("TMB")
library('MASS')
library('spatstat')

compile('src/tmb_ICAR.cpp')
dyn.load(dynlib('src/tmb_ICAR'))

ssnpath <- "data/maris_ssn/maris4.ssn"
spp <- "Brook trout"
epsilon <- 1e-6

source('code/cox.variational.ICAR.R')
source('code/kernel.R')
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

# import the data
maris <- importSSN(ssnpath, predpts='preds')


# import the data supplement that gives full coverage for the quadrature knots
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

s2 <- rbind(maris@predpoints@SSNPoints[[1]]@point.data, supplement)

maris@predpoints@SSNPoints[[1]]@point.data <- s2
maris@predpoints@SSNPoints[[1]]@point.coords <- as.matrix(s2[,c('NEAR_X', 'NEAR_Y')])
colnames(maris@predpoints@SSNPoints[[1]]@point.coords) <- c('coords.x1', 'coords.x2')
maris@predpoints@SSNPoints[[1]]@network.point.coords <- s2[,c('netID', 'rid', 'upDist')]
colnames(maris@predpoints@SSNPoints[[1]]@network.point.coords) <- c('NetworkID', 'SegmentID', 'DistanceUpstream')


# Get the IDs of distinct networks in the prediction dataset
netids.pred <- unique(maris@predpoints@SSNPoints[[1]]@point.data$netID)
netids.pred <- as.integer(levels(netids.pred)[netids.pred])

# Get the IDs of distinct networks in the observation data
netids.obs <- unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)
netids.obs <- as.integer(levels(netids.obs)[netids.obs])

# import stream shapefiles
edges <- readOGR(paste(ssnpath, 'edges.shp', sep='/'), layer='edges')
netids <- unique(edges@data$netID)

# correct some small errors in the data
edges@data$upDist[edges@data$netID == 30] <- edges@data$upDist[edges@data$netID == 30] + edges@data$upDist[edges@data$rid == 12473]
edges@data$upDist[edges@data$netID == 84] <- edges@data$upDist[edges@data$netID == 84] + edges@data$upDist[edges@data$rid == 20199]
edges@data$upDist[edges@data$netID == 87] <- edges@data$upDist[edges@data$netID == 87] + edges@data$upDist[edges@data$rid == 20212]
edges@data$upDist[edges@data$netID == 91] <- edges@data$upDist[edges@data$netID == 91] + edges@data$upDist[edges@data$rid == 22498]
edges@data$netID[edges@data$netID == 30] <- 82
edges@data$netID[edges@data$netID == 84] <- 82
edges@data$netID[edges@data$netID == 87] <- 86
edges@data$netID[edges@data$netID == 91] <- 90

# drop the reach with rid 39142 because it connects two points that are also connected by 39143
indx <- edges@data$rid != 39142
edges@lines <- edges@lines[indx]
edges@data <- edges@data[indx,]

# import the binary codes that indicate flow connections in a stream network, and remove whitespace
connect <- dbConnect(RSQLite:::SQLite(), paste(ssnpath, "binaryID.db", sep='/'))
bincodes <- list()
meta <- matrix(NA, 0, 2)
run <- 0
for (netID in netids) {
    bincodes[[netID]] <- dbReadTable(connect, paste('net', netID, sep=''))
    bincodes[[netID]]$binaryID <- sapply(bincodes[[netID]]$binaryID, function(z) stri_replace_all_charclass(z, "\\p{WHITE_SPACE}", ""))
    run <- run + sum(edges@data$Length[edges@data$netID == netID])
    
    meta <- rbind(meta, cbind(bincodes[[netID]], netID))
}
meta <- as.data.frame(meta)

# Fix some minor data errors:

# network 30 is actually an offshoot of network 82
bincodes[[30]]$binaryID <- paste0(bincodes[[82]]$binaryID[bincodes[[82]]$rid == 12473], bincodes[[30]]$binaryID)
bincodes[[82]] <- rbind(bincodes[[82]], bincodes[[30]])
bincodes[[30]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 30] <- 82
meta$binaryID[meta$netID == 82] <- bincodes[[82]]$binaryID

# network 84 is actually an offshoot of network 82
bincodes[[84]]$binaryID <- paste0(bincodes[[82]]$binaryID[bincodes[[82]]$rid == 20199], bincodes[[84]]$binaryID)
bincodes[[82]] <- rbind(bincodes[[82]], bincodes[[84]])
bincodes[[84]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 84] <- 82
meta$binaryID[meta$netID == 82] <- bincodes[[82]]$binaryID

# network 87 is actually an offshoot of network 86
bincodes[[87]]$binaryID <- paste0(bincodes[[86]]$binaryID[bincodes[[86]]$rid == 20212], c('2', '20'))
bincodes[[86]] <- rbind(bincodes[[86]], bincodes[[87]])
bincodes[[87]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 87] <- 86
meta$binaryID[meta$netID == 86] <- bincodes[[86]]$binaryID

# network 91 is actually an offshoot of network 90
bincodes[[91]]$binaryID <- paste0(bincodes[[90]]$binaryID[bincodes[[90]]$rid == 22498], bincodes[[91]]$binaryID)
bincodes[[90]] <- rbind(bincodes[[90]], bincodes[[91]])
bincodes[[91]] <- data.frame(rid=vector(), binaryID=vector())
meta$netID[meta$netID == 91] <- 90
meta$binaryID[meta$netID == 90] <- bincodes[[90]]$binaryID

# get data on the points from the SSN object
ptdata <- maris@obspoints@SSNPoints[[1]]@point.data
ptdata$netID <- as.integer(levels(ptdata$netID)[ptdata$netID])
ptdata <- ptdata %>% filter(comm_name == spp) %>% distinct(.keep_all=TRUE)

# fix some data errors
ptdata$netID[ptdata$netID == 30] <- 82
ptdata$netID[ptdata$netID == 84] <- 82
ptdata$netID[ptdata$netID == 87] <- 86
ptdata$netID[ptdata$netID == 91] <- 90

# This is an attempt to efficiently establish a basis of step functions on the stream network.
# Begin with one branch that covers the entirety of each network with at least one observation of the target species
basis <- edges@data[, c('netID', 'rid', 'Length', 'upDist')]
# basis <- edges@data[edges@data$netID %in% unique(ptdata$netID), c('netID', 'rid', 'Length', 'upDist')]
indx <- match(basis$rid, meta$rid)
basis$binaryID <- meta$binaryID[indx]
basis$binlen <- nchar(basis$binaryID)
#basis <- basis[basis$netID %in% unique(ptdata$netID), ]

indx <- match(edges@data$rid, meta$rid)
edges@data$binaryID <- meta$binaryID[indx]

# # Make a separate list of branches that summarizes all the runs making up the branch
# branches <- data.frame(netID = unique(basis$netID[basis$netID %in% unique(ptdata$netID)]))
# branches <- branches[order(branches$netID),, drop=FALSE]
# branches$br <- 1:nrow(branches)
# branches$outlet <- NA
# 
# # Copy branch ids to the basis dataframe
# match.indx <- match(basis$netID, branches$netID)
# basis$br <- branches$br[match.indx]
# # basis <- left_join(basis, branches, by='netID')
# 
# # Find the length of each branch
# branches$Length <- sapply(branches$br, function(br) sum(basis$Length[basis$br == br], na.rm=TRUE))
# branches$ptcount <- sapply(branches$netID, function(z) sum(ptdata$netID == z, na.rm=TRUE))  
# 
# # Maximum number of branches to have at the end
# br.max <- 20000
# 
# # Iterate until we've generated br.max branches
# while (nrow(branches) < br.max) {
#     cat(paste(nrow(branches), '\n', sep=''))
#     
#     # We will break the longest branch in (roughly) half
#     # But only break branches that contain at least one point
#     pdx <- branches$ptcount > 0
#     longest <- which.max(branches$Length[pdx])
#     indx.branch <- !is.na(basis$br) & basis$br == branches$br[pdx][longest]
#     branch <- basis[indx.branch, ]
#     
#     # Joints will contain all the possible splits, with the info needed to pick one
#     joints <- matrix(NA, 0, 2)
#     
#     for (run in branch$rid) {
#         # Get the runs that are at least as high upstream as the junction (distance upstream measured by number of branch points)
#         row <- which(branch$rid == run)
#         above <- branch[branch$binlen >= branch$binlen[row] & branch$upDist >= branch$upDist[row], ]
#         
#         
#         # check whether there are any runs above the current one
#         if (nrow(above) > 1) {
#             # identify all the runs that split left (right) from the current run
#             id <- above$rid[substr(above$binaryID, 1, branch$binlen[row]) == branch$binaryID[row]]
#             
#             # sum up the lengths of all runs that split left (right) from the current run
#             up <- sum(above$Length[above$rid %in% id]) + branch$Length[row]
#         } else {
#             up <- branch$Length[row]
#         }
#         
#         # add the possible left and right splits to the list of joints
#         joints <- rbind(joints, c(run, up))
#         
#         # Stop early if we're within 10% of target
#         if (abs(up * 2 / branches$Length[pdx][longest] - 1) < 0.1)
#             break()
#     }
#     
#     # post-processing the joints data
#     colnames(joints) <- c('rid', 'Length')
#     joints <- as.data.frame(joints)
#     
#     # cut at the joint that most nearly divides the length in two.
#     cut <- which.min(abs(joints$Length - branches$Length[pdx][longest]/2))
#     
#     # Cut the branch in two, and add the newly created branch to our data.frame:
#     run <- joints$rid[cut]
#     row <- which(branch$rid == run)
#     indx.above <- branch$binlen >= branch$binlen[row] & branch$upDist >= branch$upDist[row]
#     above <- branch[indx.above, ]
#     indx.new <- substr(above$binaryID, 1, branch$binlen[row]) == branch$binaryID[row]
#     
#     # Add the new branch
#     newbranch <- branches[pdx,][longest,]
#     newbranch$br <- max(branches$br) + 1
#     newbranch$Length <- joints$Length[cut]
#     basis$br[indx.branch][indx.above][indx.new] <- newbranch$br
#     newbranch$ptcount <- sum(ptdata$rid %in% basis$rid[!is.na(basis$br) & basis$br == newbranch$br])
# 
#     # Find the outlet of the new branch
#     indx.newbranch <- !is.na(basis$br) & (basis$br == newbranch$br)
#     b1 <- basis$binaryID[indx.newbranch & basis$binlen == min(basis$binlen[indx.newbranch])]
#     b0 <- substr(b1, 1, nchar(b1)-1)
#     outlet.id <- which(basis$netID == newbranch$netID & basis$binaryID == b0)
#     newbranch$outlet <- ifelse(length(outlet.id)==1, basis$rid[outlet.id], NA)
#     branches <- rbind(branches, newbranch)
# 
#     # Adjust the cut branch
#     branches$Length[pdx][longest] <- branches$Length[pdx][longest] - joints$Length[cut]
#     branches$ptcount[pdx][longest] <- sum(ptdata$rid %in% basis$rid[basis$br == branches$br[pdx][longest]])
# }




# estimate the sampling effort by the density of all fish observations
# pointlocs <- maris@obspoints@SSNPoints[[1]]@point.coords
# pointlocs <- as.data.frame(pointlocs[!duplicated(paste0(pointlocs[,1], pointlocs[,2])),])
# colnames(pointlocs) <- c('X', 'Y')
# 
# # estimate the sampling effort at the point locations
# h <- c(bandwidth.nrd(pointlocs$X), bandwidth.nrd(pointlocs$Y)) / 4
# gx <- maris@obspoints@SSNPoints[[1]]@point.coords[,1]
# gy <- maris@obspoints@SSNPoints[[1]]@point.coords[,2]
# ax <- outer(gx, pointlocs$X, "-") / h[1L]
# ay <- outer(gy, pointlocs$Y, "-") / h[2L]
# z <- rowSums(dnorm(ax) * dnorm(ay)) / (nrow(pointlocs) * h[1L] * h[2L])
# offset.cox <- log(z)
# 
# 
# # estimate the sampling effort on the quarature knots
# stepsize <- 20000
# indx.begin <- 1
# indx.end <- min(stepsize, nrow(maris@predpoints@SSNPoints[[1]]@point.coords))
# z <- vector()
# 
# while(indx.begin <= nrow(maris@predpoints@SSNPoints[[1]]@point.coords)) {
#   gx <- maris@predpoints@SSNPoints[[1]]@point.coords[indx.begin:indx.end, 1]
#   gy <- maris@predpoints@SSNPoints[[1]]@point.coords[indx.begin:indx.end, 2]
#   ax <- outer(gx, pointlocs$X, "-") / h[1L]
#   ay <- outer(gy, pointlocs$Y, "-") / h[2L]
#   z <- c(z, rowSums(dnorm(ax) * dnorm(ay)) / (nrow(pointlocs) * h[1L] * h[2L]))
#   
#   indx.begin <- indx.end + 1
#   indx.end <- min(indx.end + stepsize, nrow(maris@predpoints@SSNPoints[[1]]@point.coords))
# }
# offset.q <- log(z)













# this is the key that matches rid to seg
rid <- basis$rid
rid.lookup <- data.frame(rid=rid, seg=1:length(rid))
colnames(rid.lookup) <- c('rid', 'seg')







# adjust the weighting of individual reaches and quadrature knots because some reaches aren't full covered within one bandwidth of an observation point.

# coverage is the proportion of the reach that is covered by nonzero sampling intensity
coverage <- vector()

for (rid in obs.rids) {
  edgelength <- edges$Length[edges$rid == rid]
  
  if (edgelength < bw) {
    covered <- 1
  } else {
    
    frac <- bw / edgelength
    pts <- sort(ssn.kern@obspoints@SSNPoints[[1]]@point.data$ratio[ssn.kern@obspoints@SSNPoints[[1]]@point.data$rid == rid])
    
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
    
    covered <- covered + (frontier - start)
  }
  
  coverage <- c(coverage, covered)
}

indx <- match(obs.rids, edges$rid)
edges@data$factor <- 1
edges@data$factor[indx] <- 1 / coverage






# get data from the observation locations and the quarature knots
X.cox <- maris@obspoints@SSNPoints[[1]]@point.data
X.cox$offset <- offset.cox
X.cox$tp <- 1 - X.cox$ratio
X.q <- maris@predpoints@SSNPoints[[1]]@point.data
X.q$offset <- offset.q
X.q$tp <- 1 - X.q$ratio

# filter out rid 39142
X.q <- X.q[X.q$rid %in% edges@data$rid,]

# only nets where points were observed
X.q <- X.q[X.q$netID %in% obs.nets, ]


##############################################
# subsetting to only the even-numbered nets:
# even.nets <- obs.nets[(obs.nets %% 2) == 0]
# X.q <- X.q[X.q$netID %in% even.nets,]
# X.cox <- X.cox[X.cox$netID %in% even.nets, ]

# Drop any networks where the only remaining data point is a presence point.
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


# filter out a single fish species and collapse locations with multiple points
X.cox <- X.cox %>% filter(comm_name == spp) %>% distinct(rid, ratio, .keep_all=TRUE)

# filter out quadrature knots on networks with no points
X.q <- X.q[X.q$netID %in% unique(X.cox$netID),]
# 
# indx <- match(X.cox$rid, basis$rid)
# X.cox$br <- basis$br[indx]
# 
# indx <- match(X.q$rid, basis$rid)
# X.q$br <- basis$br[indx]
# X.cox <- left_join(X.cox, basis[,c('rid', 'br')], by='rid')
# X.q <- left_join(X.q, basis[,c('rid', 'br')], by='rid')

# Code the branch of quadrature knots on networks without any points as -1
# This will cause these networks to be assigned no random effect, since the
# random effects are assigned from 1 to max(XX$br)
# X.q$br[is.na(X.q$br)] <- -1

# use only columns that are in both the obs and pred data
common <- names(X.cox)[which(sapply(names(X.cox), function(z) z %in% names(X.q)))]
XX <- rbind(X.cox[,common], X.q[,common])
XX$netID <- as.integer(levels(XX$netID)[XX$netID])
loc <- data.frame(netID=XX$netID, rid=XX$rid, upDist=XX$upDist)

# -9999 and -9998 are how NAs are coded. Convert these back to NAs.
X.cox[X.cox == -9999 | X.cox == -9998 | X.cox == Inf | X.cox == -Inf] <- NA
X.q[X.q == -9999 | X.q == -9998 | X.q == Inf | X.q == -Inf] <- NA
XX[XX == -9999 | XX == -9998 | XX == Inf | XX == -Inf] <- NA

# Use only the complete cases
X.cox <- X.cox[complete.cases(X.cox),]
X.q <- X.q[complete.cases(X.q),]
XX <- XX[complete.cases(XX),]

# count the observations and quadrature knots
n.q <- nrow(X.q)
n.cox <- nrow(X.cox)
n <- n.cox + n.q

# add tp and seg to the data and create a temporary object where we can change the names of x and y coordinates
# wyom.q <- X.q
# colnames(wyom.q)[colnames(wyom.q) == 'NEAR_X'] <- 'x'
# colnames(wyom.q)[colnames(wyom.q) == 'NEAR_Y'] <- 'y'
# qpoints <- lpp(X=wyom.q, L=net)

seg.q <- data.frame(seg = X.q$seg)
seg.q <- seg.q %>% group_by(seg) %>% summarize(count = n())
count <- rep(0, nrow(rid.lookup))
count[seg.q$seg] <- seg.q$count
rid.lookup$count <- count

match.indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$wt <- edges@data$Length[match.indx] / edges@data$factor / rid.lookup$count

match.indx <- match(X.q$rid, rid.lookup$rid)
X.q$wt <- rid.lookup$wt[match.indx]



# We want to log the streamflow, so can't have zeroes.
# Any zero streamflow measurements are set to the minimum nonzero streamflow.
XX$MS_Hist[XX$MS_Hist <= 0] <- min(XX$MS_Hist[XX$MS_Hist > 0])
XX$logMS_Hist <- log(XX$MS_Hist)

X.cox$MS_Hist[X.cox$MS_Hist <= 0] <- min(X.cox$MS_Hist[X.cox$MS_Hist > 0])
X.cox$logMS_Hist <- log(X.cox$MS_Hist)


# set the presence and quadrature weights
presence <- c(rep(1, n.cox), rep(0, n.q))
p.wt <- c(rep(epsilon, n.cox), X.q$wt)

# create a data frame for the model
dat <- data.frame(resp=presence/p.wt, wt=p.wt)
dat <- within(dat, {
              logMS_Hist <- poly(XX$logMS_Hist, 1);
              #W95_Hist <- poly(XX$W95_Hist, 1);
              #Elev_NSI <- poly(XX$Elev_NSI, 1);
              S1_93_11 <- poly(XX$S1_93_11, 2);
              SLOPE <- poly(XX$SLOPE, 1);
              #CFM_Hist <- poly(XX$CFM_Hist, 1);
              offset <- poly(XX$offset, 1);
              })

rownames(dat) <- 1:n
dat$id <- rownames(dat)

# estimate the DWPR model
dwpr <- glm(resp ~ logMS_Hist + S1_93_11 + SLOPE + offset, family=poisson(), weights=wt, data=dat)
dwpr$loglik <- sum(log(dwpr$fitted.values[1:n.cox])) - n.cox
dwpr$pres <- dat$pres



#Identify networks that have no points for this species
#X.cox <- X.cox %>% filter(comm_name == spp) %>% group_by(locID) %>% distinct

# 
# # create a matrix of spatial random effects, using 'breaks' which are intended to have roughly equal total stream length
# Z <- Matrix(0, nrow=nrow(XX), ncol=max(XX$br))
# for (i in 1:ncol(Z)) {
#     Z[, i] <- ifelse(XX$br == i, 1, 0)
# }
# 
# # Convert the flow-connection matrix to the initial covariance matrix
# # identify the potentially nonzero entries in the covariance matrix
# # covariance can be nonzero for branches that are on the same network
# # Initialize the diagonal to be larger than the other entries to ensure positive-definiteness
# Sigma.mask <- Matrix(0, nrow=ncol(Z), ncol=ncol(Z))
# for (bid in unique(XX$br[XX$br>=0])) {
#   br.loc <- basis[!is.na(basis$br) & basis$br == bid,]
#   nid <- br.loc$netID[1]
#   
#   #allow covariance between this branch and its downstream outlet
#   outlet <- branches$outlet[branches$br == bid]
#   downstream <- ifelse(is.na(outlet), NA, basis$br[!is.na(basis$br) & basis$rid == outlet])
#   if (!is.na(downstream))
#     Sigma.mask[bid, downstream] <- -1
# }
# 
# Sigma.mask <- Sigma.mask + t(Sigma.mask)
# diag(Sigma.mask) <- -colSums(Sigma.mask)
# 
# 
# # Drop branches that don't appear in the final data set (because all of their rows were incomplete data cases)
# keep.indx <- colSums(Z) > 0
# Z <- Z[,keep.indx]
# Sigma.mask <- Sigma.mask[keep.indx, keep.indx]
# to.drop <- which(!keep.indx)
# if (length(to.drop) > 0) {
#     branches$br[branches$br %in% to.drop] <- -1
#     basis$br[basis$br %in% to.drop] <- -1
#     XX$br[XX$br %in% to.drop] <- -1
# }

# convert to "triples" form and then get the coordinates of nonzero entries:
# Sigma.mask <- as(Sigma.mask, 'dgTMatrix')
# xindx <- Sigma.mask@i
# yindx <- Sigma.mask@j
# logV <- rep(0, ncol(Z))
# sigma <- max(diag(Sigma.mask)) / 2
# cov <- max(Sigma.mask[upper.tri(Sigma.mask)])
#logV <- log(Sigma.mask@x)

#Z <- model.matrix(~br, data=XX, contrasts.arg = list(br='contr.treatment', contrasts=FALSE, sparse=TRUE))[,-1]
X <- model.matrix(~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1), data=XX)
# X <- as.matrix(cbind(1, dat[,3:6]))
# X <- X[,c('1', 'MS_Hist', 'S1_93_11.1', 'S1_93_11.2', 'SLOPE', 'offset')]

formula.dwpr <- resp ~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1)






# extract the endpoints of the stream segments
ends <- sapply(edges@lines, function(x) c(head(x@Lines[[1]]@coords, 1), tail(x@Lines[[1]]@coords, 1))) %>% t
ends <- cbind(ends, edges@data$rid, edges@data$netID)
ends <- as.data.frame(rbind(ends[,c(1, 2, 5, 6)], ends[,c(3, 4, 5, 6)]))
colnames(ends) <- c('x', 'y', 'rid', 'netID')

# fix some errors in the data
ends$netID[ends$netID == 30] <- 82
ends$netID[ends$netID == 84] <- 82
ends$netID[ends$netID == 87] <- 86
ends$netID[ends$netID == 91] <- 90



# create a ppp object holding the vertices
win <- owin(xrange = range(ends$x), yrange = range(ends$y))
pts <- ends %>% group_by(x, y) %>% distinct(.keep_all = TRUE)
pts$pid <- 1:nrow(pts)
pts$rid <- NULL
ends <- right_join(pts, ends, by=c('x', 'y', 'netID'))
verts <- ppp(x=pts$x, y=pts$y, window=win)


# create the linnet object defining the network
m <- Matrix(FALSE, nrow(pts), nrow(pts)) # matrix indicating which points are connected
m.rid <- Matrix(0, nrow(pts), nrow(pts)) # matrix indicating the reach id of connected points

# identify which reaches touch each other
for (pp in unique(ends$pid)) {
  cat(paste(pp, '\n'))
  tangent <- ends$rid[ends$pid == pp]
  connected <- ends$pid[ends$rid %in% tangent & ends$pid != pp]
  m[pp, connected] <- m[connected, pp] <- TRUE
  
  m.rid[pp, ends$pid[ends$rid %in% tangent & ends$pid != pp]] <- ends$rid[ends$rid %in% tangent & ends$pid != pp]
}




# create a spatstat linear network object using sparse matrices
# create a linear network object using sparse matrices
net <- MyLinnet(vertices=verts, m=m, m.rid=m.rid, sparse=TRUE, warn=FALSE)

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







# estimate the DWPR model
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

bw <- 2000
qq <- quad.dwpr(obj)
qq$covars$offset <- kernel.wrap(kern, rid.lookup, ssn.kern@obspoints@SSNPoints[[1]]@point.data, edges, bw, nzw.rids=obs.rids, logresult=TRUE)

# add the segment ID as used by spatstat functions to the Brown trout object
bt.lpp <- lpp(X.cox, net)
obj$dwpr <- lppm.ssnpp(bt.lpp, qq)


resp <- presence / p.wt
beta <- dwpr$coefficients


# create a matrix of spatial random effects
Z <- Matrix(0, nrow=nrow(XX), ncol=nrow(gmrf.lookup))

# populate the random effects design matrix
rid.indx <- match(XX$rid, gmrf.lookup$rid)
colID <- gmrf.lookup$gmrfID[rid.indx]
# rowID <- XX$rid == edges@data$rid[i]
vecID <- nrow(Z) * (colID-1) + (1:nrow(Z))
Z[vecID] <- 1


indx.evens <- gmrf.lookup$netID %in% XX$netID
Z <- Z[,indx.evens]
Sigma.evens <- Sigma.mask[indx.evens, indx.evens]
eigenvalues.evens <- eigenvalues[netlookup %in% unique(XX$netID)]

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


# for (i in 1:length(A.indx)) {
#   lagrange.matrix[i,] <- ifelse(A.indx == i-1, 1, 0)
# }

# estimate the Cox process model
# wes.time <- system.time(vari2 <- cox.variational(resp, X, Z, xindx, yindx, p.wt, beta, logV.start = logV, rho.start = 0.01, sigma.start=sigma, hess=TRUE, sd=FALSE))
logV <- rep(0.01, ncol(Z))
# A <- rep(0.01, length(A.lookup))
# wes.time <- system.time(vari2 <- cox.variational.ICAR(resp, X, Z,  Sigma.evens, sindx, p.wt, beta.start=vari2$beta, u.start=vari2$M, ltau.start=vari2$ltau, logDetQ=sum(log(eigenvalues.evens[eigenvalues.evens > sqrt(.Machine$double.eps)])), logV.start = vari2$logV, hess=FALSE, sd=FALSE, tol=.Machine$double.eps, maxit=1e6))
 
wes.time <- system.time(vari2 <- cox.variational.ICAR(resp, X, Z,  Sigma.evens, sindx, p.wt, beta.start=dwpr$coefficients, means.matrix=netmeans.matrix, u.start=rep(0.01, ncol(Z)), ltau.start=obj.2k$cox$ltau, logDetQ=sum(log(eigenvalues.evens[eigenvalues.evens > sqrt(.Machine$double.eps)])), logV.start = rep(0.01, ncol(Z)), hess=FALSE, sd=FALSE, tol=.Machine$double.eps, maxit=1e7))


# trace <- rbind(trace, c(br.max, vari2$object$objective))
# print(trace)

vari2$time <- wes.time

obj$cox <- vari2





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

# use linear combination matrix to transform the coefficients and covariances
eee <- as.vector(lincomb %*% obj$cox$beta)[6]



#################################
# find the maximum value of lambda
std <- exp(obj$cox$logV / 2)
u.realized <- rnorm(length(std), obj$cox$M, std) - std^2 / 2  #random effects corrected for log transformation bias
# 
# # identify which random effect goes with which entry
# # branches with negative id get no random effect
# brange <- rep(NA, max(rid.lookup$br))
# valid <- unique(rid.lookup$br)
# valid <- sort(valid[valid > 0])
# brange[valid] <- 1:length(valid)
# 
# # match segment ids to branch ids
# indx <- match(XX$seg, rid.lookup$seg)
# braw <- rid.lookup$br[indx]
# bflag <- braw > 0
# indx <- cbind((1:length(XX$seg))[bflag], brange[braw[bflag]])
# 
# # populate the random design matrix
# S <- Matrix(0, nrow=nrow(XX), ncol=length(u.realized))
# S[indx] <- 1

# calculate the intensity at aa observation & quadrature locations in order to get the lmax
fixed <- as.vector(X %*% vari2$beta)
noeffort <- as.vector(X[,1:5] %*% vari2$beta[1:5])
ranef <- as.vector(Z %*% u.realized)
lambda.noeffort <- exp(noeffort + ranef - mean(XX$offset) * eee)
lambda <- exp(fixed + ranef)
lmax <- max(lambda)

# This version doesn't include the intensity because we want to calculate that at the 
#####################################

lamdat <- data.frame(x=XX$NEAR_X, y=XX$NEAR_Y, seg=XX$seg, tp=XX$seg, lambda=lambda.noeffort)
lam <- linfun.factory(lamdat, 'lambda')
    

# compute the inhomogeneous K function
indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$netID <- edges@data$netID[indx]
rid.lookup$Length <- edges@data$Length[indx]
rid.lookup$upDist <- edges@data$upDist[indx]
mylpp <- lpp(X.cox, net)
inhom <- linearKinhom2(mylpp, net, lambda, lookup=rid.lookup, ssn=maris, edges=edges@data)

kf <- kernel.wrap(kern, rid.lookup, ssn.kern@obspoints@SSNPoints[[1]]@point.data, edges, bw, nzw.rids=obs.rids, distpath='~/streams/data/maris_ssn/maris4.ssn/distance/obs/')
eff <- effort.wrapper.kernel(edges@data, rid.lookup, eee, kf)
sim <- simulate2.lppm(obj$dwpr, lambda=lam, lmax=lmax, effort=eff, nzw.rids=obs.rids, loglambda=FALSE)


 
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


# find the hessian of the fixed-effect coefficients by wrapping the objective and 
# gradient functions so that only the fixed effets can vary.
# training: allow the regression coefficients to vary. parameters masked TRUE can vary.
paramselect <- c(TRUE, rep(FALSE, length(vari2$M)), rep(FALSE, length(vari2$M)), rep(TRUE, length(beta)))

subsetwrapper <- function(fun, mask, masked.pars, vector.result=FALSE) {
  paramvec <- masked.pars
  
  out <- function(x=NULL) {
    if (!is.null(x))
      paramvec[mask] <- x
    
    if (vector.result)
      fun(paramvec)[,mask, drop=FALSE]
    else
      fun(paramvec)
  }
}

fn <- subsetwrapper(obj$cox$obj$fn, paramselect, obj$cox$res$par)
gr <- subsetwrapper(obj$cox$obj$gr, paramselect, obj$cox$res$par, vector.result=TRUE)
constrained <- nlminb(c(obj$cox$ltau, obj$cox$beta), objective=fn, gradient=gr, control=list(iter.max=1e7, eval.max=1e7, rel.tol=.Machine$double.eps))
constrained.hess <- Matrix(optimHess(constrained$par, fn=fn, gr=gr, control=list(reltol=.Machine$double.eps)))


# use linear combination matrix to transform the coefficients and covariances
beta <- lincomb %*% obj$cox$beta
covmat <- lincomb %*% solve(constrained.hess[2:7,2:7]) %*% t(lincomb)

# report the back-transformed regression coefficients and their standard errors
coefficients <- data.frame(beta=lincomb %*% obj$cox$beta, sd=sqrt(tail(diag(covmat), 6)))
print(cbind(colNames(X), round(coefficients, 3)))

# report the back-transpormed variance component for the random field
print("mean of the estimated variance component: ", mean(1/exp(rnorm(1000, mean=obj$cox$ltau, sd=1/diag(constrained.hess)[1]**0.5))))
print("mean of the estimated variance component: ", sd(1/exp(rnorm(1000, mean=obj$cox$ltau, sd=1/diag(constrained.hess)[1]**0.5)))

lincomb %*% coef(obj$dwpr$fit)
sqrt(diag(lincomb %*% summary(obj$dwpr$fit$internal$glmfit)$cov.unscaled %*% t(lincomb))) %>% round(4)
