library("SSN")
library('dplyr')
library('sp')
library('rgdal')
library('stringi')
library('Matrix')
library('MASS')
library('spatstat')

ssnpath <- "data/maris_ssn/maris4.ssn"
spp <- "Brown trout"
epsilon <- 1e-6

# import the data
maris <- importSSN(ssnpath, predpts='preds')

NO_OFFSET_HACK <- TRUE

# import the data supplement that gives full coverage for the quadrature knots
supplement <- read.csv("data/quadpts.csv")
supplement$OBSPRED_ID <- 1:nrow(supplement) + max(maris@predpoints@SSNPoints[[1]]@point.data$OBSPRED_ID)
colnames(supplement)[colnames(supplement) %in% c('POINT_X', 'POINT_Y')] <- c('NEAR_X', 'NEAR_Y')
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


# drop the reach with rid 39142 because it connects two points that are also connected by 39143
indx <- edges@data$rid != 39142
edges@lines <- edges@lines[indx]
edges@data <- edges@data[indx,]
supplement <- supplement[supplement$rid != 39142, ]
bincodes[[253]] <- bincodes[[253]][bincodes[[253]]$rid != 39142,]
meta <- meta[meta$rid != 39142, ]


# join the supplement and the original maris data
maris@predpoints@SSNPoints[[1]]@point.data <- rbind(maris@predpoints@SSNPoints[[1]]@point.data, supplement)
maris@predpoints@SSNPoints[[1]]@point.coords <- rbind(maris@predpoints@SSNPoints[[1]]@point.coords, as.matrix(supplement[,c('NEAR_X', 'NEAR_Y')]))
colnames(maris@predpoints@SSNPoints[[1]]@point.coords) <- c('coords.x1', 'coords.x2')
npc <- data.frame(DistanceUpstream = c(maris@predpoints@SSNPoints[[1]]@network.point.coords$DistanceUpstream, supplement$upDist))
npc$NetworkID <- as.factor(c(as.character(maris@predpoints@SSNPoints[[1]]@network.point.coords$NetworkID), as.character(supplement$netID)))
npc$SegmentID <- as.factor(c(as.character(maris@predpoints@SSNPoints[[1]]@network.point.coords$SegmentID), as.character(supplement$rid)))
maris@predpoints@SSNPoints[[1]]@network.point.coords <- npc[, c('NetworkID', 'SegmentID', 'DistanceUpstream')]


# Get the IDs of distinct networks in the prediction dataset
netids.pred <- unique(maris@predpoints@SSNPoints[[1]]@point.data$netID)
netids.pred <- as.integer(levels(netids.pred)[netids.pred])

# Get the IDs of distinct networks in the observation data
netids.obs <- unique(maris@obspoints@SSNPoints[[1]]@point.data$netID)
netids.obs <- as.integer(levels(netids.obs)[netids.obs])


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
ptdata <- ptdata %>% filter(comm_name == spp) %>% group_by(locID) %>% distinct()

# fix some data errors
ptdata$netID[ptdata$netID == 30] <- 82
ptdata$netID[ptdata$netID == 84] <- 82
ptdata$netID[ptdata$netID == 87] <- 86
ptdata$netID[ptdata$netID == 91] <- 90

# This is an attempt to efficiently establish a basis of step functions on the stream network.
# Begin with one branch that covers the entirety of each network with at least one observation of the target species
basis <- edges@data[edges@data$netID %in% unique(ptdata$netID), c('netID', 'rid', 'Length', 'upDist')]
indx <- match(basis$rid, meta$rid)
basis$binaryID <- meta$binaryID[indx]
basis$binlen <- nchar(basis$binaryID)
basis <- basis[basis$netID %in% unique(ptdata$netID), ]

# Make a separate list of branches that summarizes all the runs making up the branch
branches <- data.frame(netID = unique(basis$netID))
branches <- branches[order(branches$netID),, drop=FALSE]
branches$br <- 1:nrow(branches)

# Copy branch ids to the basis dataframe
basis <- left_join(basis, branches, by='netID')

# Find the length of each branch
branches$Length <- sapply(branches$br, function(br) sum(basis$Length[basis$br == br]))
branches$ptcount <- sapply(branches$netID, function(z) sum(ptdata$netID == z))

# Maximum number of branches to have at the end
br.max <- 200

# Iterate until we've generated br.max branches
while (nrow(branches) < br.max) {
    cat(paste(nrow(branches), '\n', sep=''))
    
    # We will break the longest branch in (roughly) half
    # But only break branches that contain at least one point
    pdx <- branches$ptcount > 0
    longest <- which.max(branches$Length[pdx])
    indx.branch <- basis$br == branches$br[pdx][longest]
    branch <- basis[indx.branch, ]
    
    # Joints will contain all the possible splits, with the info needed to pick one
    joints <- matrix(NA, 0, 2)
    
    for (run in branch$rid) {
        # Get the runs that are at least as high upstream as the junction (distance upstream measured by number of branch points)
        row <- which(branch$rid == run)
        above <- branch[branch$binlen >= branch$binlen[row] & branch$upDist >= branch$upDist[row], ]

        
        # check whether there are any runs above the current one
        if (nrow(above) > 1) {
            # identify all the runs that split left (right) from the current run
            id <- above$rid[substr(above$binaryID, 1, branch$binlen[row]) == branch$binaryID[row]]

            # sum up the lengths of all runs that split left (right) from the current run
            up <- sum(above$Length[above$rid %in% id]) + branch$Length[row]
        } else {
            up <- branch$Length[row]
        }
        
        # add the possible left and right splits to the list of joints
        joints <- rbind(joints, c(run, up))
        
        # Stop early if we're within 10% of target
        if (abs(up * 2 / branches$Length[pdx][longest] - 1) < 0.1)
            break()
    }
    
    # post-processing the joints data
    colnames(joints) <- c('rid', 'Length')
    joints <- as.data.frame(joints)
    
    # cut at the joint that most nearly divides the length in two.
    cut <- which.min(abs(joints$Length - branches$Length[pdx][longest]/2))
    
    # Cut the branch in two, and add the newly created branch to our data.frame:
    run <- joints$rid[cut]
    row <- which(branch$rid == run)
    indx.above <- branch$binlen >= branch$binlen[row] & branch$upDist >= branch$upDist[row]
    above <- branch[indx.above, ]
    id.new <- which(substr(above$binaryID, 1, branch$binlen[row]) == branch$binaryID[row])
    
    # Add the new branch
    newbranch <- branches[pdx,][longest,]
    newbranch$br <- max(branches$br) + 1
    newbranch$Length <- joints$Length[cut]
    basis$br[indx.branch][indx.above][id.new] <- newbranch$br
    newbranch$ptcount <- sum(ptdata$rid %in% basis$rid[basis$br == newbranch$br])
    branches <- rbind(branches, newbranch)
    
    # Adjust the cut branch
    branches$Length[pdx][longest] <- branches$Length[pdx][longest] - joints$Length[cut]
    branches$ptcount[pdx][longest] <- sum(ptdata$rid %in% basis$rid[basis$br == branches$br[pdx][longest]])
}


if (!NO_OFFSET_HACK) {
    # estimate the sampling effort by the density of all fish observations
    pointlocs <- maris@obspoints@SSNPoints[[1]]@point.coords
    pointlocs <- as.data.frame(pointlocs[!duplicated(paste0(pointlocs[,1], pointlocs[,2])),])
    colnames(pointlocs) <- c('X', 'Y')
    
    # estimate the sampling effort at the point locations
    h <- c(bandwidth.nrd(pointlocs$X), bandwidth.nrd(pointlocs$Y)) / 4
    gx <- maris@obspoints@SSNPoints[[1]]@point.coords[,1]
    gy <- maris@obspoints@SSNPoints[[1]]@point.coords[,2]
    ax <- outer(gx, pointlocs$X, "-") / h[1L]
    ay <- outer(gy, pointlocs$Y, "-") / h[2L]
    z <- rowSums(dnorm(ax) * dnorm(ay)) / (nrow(pointlocs) * h[1L] * h[2L])
    offset.cox <- log(z)
    
    
    # estimate the sampling effort on the quarature knots
    gx <- maris@predpoints@SSNPoints[[1]]@point.coords[,1]
    gy <- maris@predpoints@SSNPoints[[1]]@point.coords[,2]
    
    
    ax <- outer(gx[1:75000], pointlocs$X, "-") / h[1L]
    ay <- outer(gy[1:75000], pointlocs$Y, "-") / h[2L]
    z <- rowSums(dnorm(ax) * dnorm(ay)) / (nrow(pointlocs) * h[1L] * h[2L])
    offset.q <- log(z)
    
    
    ax <- outer(gx[75001:length(gx)], pointlocs$X, "-") / h[1L]
    ay <- outer(gy[75001:length(gx)], pointlocs$Y, "-") / h[2L]
    z <- rowSums(dnorm(ax) * dnorm(ay)) / (nrow(pointlocs) * h[1L] * h[2L])
    offset.q <- c(offset.q, log(z))
} else {
    # Computing effort is v. slow so this hack just assigns unit effort everywhere
    pointlocs <- maris@obspoints@SSNPoints[[1]]@point.coords
    pointlocs <- as.data.frame(pointlocs[!duplicated(paste0(pointlocs[,1], pointlocs[,2])),])
    offset.cox <- rep(1, nrow(pointlocs))
    offset.cox <- rep(1, nrow(maris@obspoints@SSNPoints[[1]]@point.coords))
}



# get data from the observation locations and the quadrature knots
X.cox <- maris@obspoints@SSNPoints[[1]]@point.data
X.cox$offset <- offset.cox
X.q <- maris@predpoints@SSNPoints[[1]]@point.data
X.q$offset <- offset.q
X.q <- X.q[X.q$rid != 39142,]

X.cox$netID <- as.integer(levels(X.cox$netID)[X.cox$netID])
X.q$netID <- as.integer(levels(X.q$netID)[X.q$netID])

# Fix some errors in the data
X.cox$netID[X.cox$netID == 30] <- 82
X.cox$netID[X.cox$netID == 84] <- 82
X.cox$netID[X.cox$netID == 87] <- 86
X.cox$netID[X.cox$netID == 91] <- 90

X.q$netID[X.q$netID == 30] <- 82
X.q$netID[X.q$netID == 84] <- 82
X.q$netID[X.q$netID == 87] <- 86
X.q$netID[X.q$netID == 91] <- 90

# log transform the flow, adding one to account for segments with no flow.
X.cox$logMS_Hist <- log(X.cox$MS_Hist)
X.q$logMS_Hist <- log(X.q$MS_Hist)
X.q$logMS_Hist[X.q$logMS_Hist == -Inf] <- X.cox$logMS_Hist[X.cox$logMS_Hist == -Inf] <- min(X.q$logMS_Hist[is.finite(X.q$logMS_Hist)]) 

# filter out a single fish species and collapse locations with multiple points
X.cox <- X.cox %>% filter(comm_name == spp) %>% group_by(locID) %>% distinct

match.indx <- match(X.cox$rid, basis$rid)
X.cox$br <- basis$br[match.indx]

match.indx <- match(X.q$rid, basis$rid)
X.q$br <- basis$br[match.indx]

# Code the branch of quadrature knots on networks without any points as -1
# This will cause these networks to be assigned no random effect, since the
# random effects are assigned from 1 to max(XX$br)
X.q$br[is.na(X.q$br)] <- -1

# use only columns that are in both the obs and pred data
common <- names(X.cox)[which(sapply(names(X.cox), function(z) z %in% names(X.q)))]
XX <- rbind(X.cox[,common], X.q[,common])
loc <- data.frame(netID=XX$netID, rid=XX$rid, upDist=XX$upDist)

# -9999 and -9998 are how NAs are coded. Convert these back to NAs.
X.cox[X.cox == -9999 | X.cox == -9998] <- NA
X.q[X.q == -9999 | X.q == -9998] <- NA
XX[XX == -9999 | XX == -9998] <- NA

# Use only the complete cases
X.cox <- X.cox[complete.cases(X.cox),]
X.q <- X.q[complete.cases(X.q),]
XX <- XX[complete.cases(XX),]

# count the observations and quadrature knots
n.q <- nrow(X.q)
n.cox <- nrow(X.cox)
n <- n.cox + n.q

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
pts <- ends %>% group_by(x, y) %>% distinct
pts$pid <- 1:nrow(pts)
pts$rid <- NULL
ends <- right_join(pts, ends, by=c('x', 'y', 'netID'))
verts <- ppp(x=pts$x, y=pts$y, window=win)


# create the linnet object defining the network
m <- Matrix(FALSE, nrow(pts), nrow(pts))
m.rid <- Matrix(0, nrow(pts), nrow(pts))

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
rid.lookup <- data.frame(rid=net$lines$rid, seg=1:net$lines$n)
indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$netID <- edges@data$netID[indx]
net$lines$Length <- rid.lookup$Length <- edges@data$Length[indx]
rid.lookup$upDist <- edges@data$upDist[indx]

# indicator of whether each line segment has a complete data case:
indx <- match(net$lines$rid, XX$rid)
net$lines$known <- !is.na(indx)

# add some metadata to the data
match.indx <- match(XX$rid, rid.lookup$rid)
XX$seg <- rid.lookup$seg[match.indx]
XX$tp <- 1 - XX$ratio

# get the quadrature weights
seg.q <- data.frame(seg = tail(XX$seg, n.q))
seg.q <- seg.q %>% group_by(seg) %>% summarize(count = n())
count <- rep(0, nrow(rid.lookup))
count[seg.q$seg] <- seg.q$count
rid.lookup$count <- count

match.indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$wt <- edges@data$Length[match.indx] / rid.lookup$count

match.indx <- match(X.q$rid, rid.lookup$rid)
wt.q <- rid.lookup$wt[match.indx]
presence <- c(rep(1, n.cox), rep(0, n.q))
p.wt <- c(rep(epsilon, n.cox), wt.q)

# create a data frame for the model
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
formula.dwpr <- resp ~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1)

# estimate the DWPR model
# dwprmodel <- glm(formula.dwpr, family=poisson(), weights=wt, data=dat)
   
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

qq <- quad.dwpr(obj)

# add the segment ID as used by spatstat functions to the Brown trout object
bt <- X.cox
bt$tp <- 1 - bt$ratio
indx <- match(bt$rid, rid.lookup$rid)
bt$seg <- rid.lookup$seg[indx]

bt.lpp <- lpp(bt, net)

obj$dwpr <- lppm.ssnpp(bt.lpp, qq)
lam <- lambda.fun.factory(obj$dwpr$fit$internal$glmfit$terms, qq$covars, obj$dwpr$fit$coef)














# Generate the matrix that converts the coefficients from orthogonal polynomials to raw polynomials
lincomb <- diag(length(obj$dwpr$fit$coef))
quad.coefs <- c('S1_93_11')
lin.coefs <- c('SLOPE', 'logMS_Hist', 'offset')
for (cc in quad.coefs) {
    locname <- paste("poly(", cc, ", 2)", sep='')
    lindx <- which(names(coef(obj$dwpr$fit)) == paste(locname, '1', sep=''))
    quadx <- which(names(coef(obj$dwpr$fit)) == paste(locname, '2', sep=''))
    
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
}

# back-tranfsorm the coefficients that were only centered and scaled:
for (cc in lin.coefs) {
    locname <- paste("poly(", cc, ", 1)", sep='')
    lindx <- which(names(coef(obj$dwpr$fit)) == locname)
    
    pars <- attr(dat[[cc]], 'coefs')
    
    # Orthogonal polynomials are a complicated transformation.
    # Add some terms to the intercept
    lincomb[1, lindx] <- -pars$alpha[1] / sqrt(pars$norm2[3])
    
    # Raw linear coefficient is a linear combination of the orthogonal polynomial linear and quadratic coefficients
    lincomb[lindx, lindx] <- 1 / sqrt(pars$norm2[3])
}



lincomb %*% coef(obj$dwpr$fit)
sqrt(diag(lincomb %*% summary(obj$dwpr$fit$internal$glmfit)$cov.unscaled %*% t(lincomb))) %>% round(4)
