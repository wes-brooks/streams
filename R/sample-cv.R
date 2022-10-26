# sample repeatedly from the predicted intensity on the networks not used to fit the model.

# need to draw the left-out random effects and use them to calculate the predictive likelihood
# first, get the precision matrix for the left-out networks.
left <- list()
network.order <- vector()
rid.order <- vector()
# cv.nets <- c(253, 23, 185, 265, 163, 95, 167, 165, 197, 187, 179, 243, 79, 83, 233, 251, 255)
for (nn in cv.nets) {
  indx <- gmrf.lookup$netID == nn
  # local.precision <- Sigma.mask[indx, indx] * exp(vari2$ltau)
  # left[[nn]] <- solve(chol(local.precision))
  
  
  left[[nn]] <- eigenstuff[[nn]]$vectors[,1:(sum(indx)-1)] %*% diag(1 / sqrt(eigenstuff[[nn]]$values[1:(sum(indx)-1)]))
  network.order <- c(network.order, rep(nn, sum(indx)))
  rid.order <- c(rid.order, gmrf.lookup$rid[indx])
}




bws <- c(50, 100, 200, 500, 1000, 2000, 5000, 10000)



# get data from the observation locations and the quarature knots
X.cox.cv <- maris@obspoints@SSNPoints[[1]]@point.data
X.cox.cv$offset <- offset.cox
X.cox.cv$tp <- 1 - X.cox.cv$ratio
X.q.cv <- maris@predpoints@SSNPoints[[1]]@point.data
X.q.cv$offset <- offset.q
X.q.cv$tp <- 1 - X.q.cv$ratio

# filter out rid 39142
X.q.cv <- X.q.cv[X.q.cv$rid %in% edges@data$rid,]


################################################
# Here we extract the cross-validation networks:
X.q.cv <- X.q.cv[X.q.cv$netID %in% cv.nets,]
X.cox.cv <- X.cox.cv[X.cox.cv$netID %in% cv.nets, ]
################################################


# specify the sequential segments
indx <- match(X.q.cv$rid, rid.lookup$rid)
X.q.cv$seg <- rid.lookup$seg[indx]
indx <- match(X.cox.cv$rid, rid.lookup$rid)
X.cox.cv$seg <- rid.lookup$seg[indx]


# filter out a single fish species and collapse locations with multiple points
X.cox.cv <- X.cox.cv %>% filter(comm_name == spp) %>% distinct(rid, ratio, .keep_all=TRUE)


# use only columns that are in both the obs and pred data
common <- names(X.cox.cv)[which(sapply(names(X.cox.cv), function(z) z %in% names(X.q.cv)))]
XX.cv <- rbind(X.cox.cv[,common], X.q.cv[,common])
XX.cv$netID <- as.integer(levels(XX.cv$netID)[XX.cv$netID])

# -9999 and -9998 are how NAs are coded. Convert these back to NAs.
X.cox.cv[X.cox.cv == -9999 | X.cox.cv == -9998 | X.cox.cv == Inf | X.cox.cv == -Inf] <- NA
X.q.cv[X.q.cv == -9999 | X.q.cv == -9998 | X.q.cv == Inf | X.q.cv == -Inf] <- NA
XX.cv[XX.cv == -9999 | XX.cv == -9998 | XX.cv == Inf | XX.cv == -Inf] <- NA

# Use only the complete cases
X.cox.cv <- X.cox.cv[complete.cases(X.cox.cv),]
X.q.cv <- X.q.cv[complete.cases(X.q.cv),]
XX.cv <- XX.cv[complete.cases(XX.cv),]

# count the observations and quadrature knots
n.q <- nrow(X.q.cv)
n.cox <- nrow(X.cox.cv)
n <- n.cox + n.q

seg.q <- data.frame(seg = X.q.cv$seg)
seg.q <- seg.q %>% group_by(seg) %>% summarize(count = n())
count <- rep(0, nrow(rid.lookup))
count[seg.q$seg] <- seg.q$count
rid.lookup$count <- count

match.indx <- match(rid.lookup$rid, edges@data$rid)
rid.lookup$wt <- edges@data$Length[match.indx] / edges@data$factor[match.indx] / rid.lookup$count

match.indx <- match(X.q.cv$rid, rid.lookup$rid)
X.q.cv$wt <- rid.lookup$wt[match.indx]



# We want to log the streamflow, so can't have zeroes.
# Any zero streamflow measurements are set to the minimum nonzero streamflow.
XX.cv$MS_Hist[XX.cv$MS_Hist <= 0] <- min(XX.cv$MS_Hist[XX.cv$MS_Hist > 0])
XX.cv$logMS_Hist <- log(XX.cv$MS_Hist)

X.cox.cv$MS_Hist[X.cox.cv$MS_Hist <= 0] <- min(X.cox.cv$MS_Hist[X.cox.cv$MS_Hist > 0])
X.cox.cv$logMS_Hist <- log(X.cox.cv$MS_Hist)


# set the presence and quadrature weights
presence.cv <- c(rep(1, n.cox), rep(0, n.q))
p.wt.cv <- c(rep(epsilon, n.cox), X.q.cv$wt)
resp.cv <- presence.cv / p.wt.cv


X.cv <- model.matrix(~ poly(logMS_Hist, 1, coefs=attr(dat$logMS_Hist, 'coefs')) +
                    poly(S1_93_11, 2, coefs=attr(dat$S1_93_11, 'coefs')) +
                    poly(SLOPE, 1, coefs=attr(dat$SLOPE, 'coefs')) +
                    poly(offset, 1, coefs=attr(dat$offset, 'coefs')), data=XX.cv)




# create a matrix of spatial random effects
Z.cv <- Matrix(0, nrow=nrow(XX.cv), ncol=nrow(gmrf.lookup))

# populate the random effects design matrix
rid.indx <- match(XX.cv$rid, gmrf.lookup$rid)
colID <- gmrf.lookup$gmrfID[rid.indx]
# rowID <- XX$rid == edges@data$rid[i]
vecID <- nrow(Z.cv) * (colID-1) + (1:nrow(Z.cv))
Z.cv[vecID] <- 1


indx.cv <- match(rid.order, gmrf.lookup$rid)
Z.cv <- Z.cv[,indx.cv]




# need to draw the left-out random effects and use them to calculate the predictive likelihood
# first, get the precision matrix for the left-out networks.
left.cv <- list()
for (nn in cv.nets) {
  indxCol <- gmrf.lookup$netID[gmrf.lookup$netID %in% cv.nets] == nn
  indxRow <- XX.cv$netID == nn
  left.cv[[nn]] <- as.matrix(Z.cv[indxRow, indxCol] %*% left[[nn]]) * exp(-obj$cox$ltau/2)
}





fix <- as.vector(X.cv %*% obj$cox$beta)

ll <- vector()
for (i in 1:10) {
  cat(i, '\n')
  N <- 1000
  
  # Sample the random field in parts. One part for each separate network.
  # z <- matrix(0, length(network.order), N)
  ran <- matrix(0, n, N)
  # u <- matrix(0, length(network.order), 1)
  
  # Gaussian log-likelihood of the random effect
  gll <- matrix(NA, 0, N)
  
  for (nn in cv.nets) {
    cat(nn, '\n')
    indx <- network.order == nn
    u <- matrix(rnorm(N * (sum(indx)-1)), nrow=sum(indx)-1, ncol=N)
    A <- rnorm(n=N, sd=exp(-obj$cox$ltau/2))
    
    indx2 <- XX.cv$netID == nn
    ran[indx2,] <- sweep(left.cv[[nn]] %*% u, 2, A, '+')
    gll <- rbind(gll, colSums(dnorm(u, log=TRUE)) + dnorm(A, sd=exp(-obj$cox$ltau/2)))
  }
  gll <- colSums(gll)
  
  
  # ran <- Z.cv %*% z
  eta <- apply(ran, 2, function(x) x+fix)
  ll <- c(ll, apply(eta, 2, function(y) sum(p.wt.cv*(resp.cv*y - exp(y)))) + gll)
  rm(ran, gll, u, eta)
  gc()
}


