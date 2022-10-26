MyLinnet <- function(vertices, m, m.rid, edges, sparse = FALSE, warn = TRUE) {
    if (missing(m) && missing(edges)) 
        stop("specify either m or edges")
    if (!missing(m) && !missing(edges)) 
        stop("do not specify both m and edges")
    stopifnot(is.ppp(vertices))
    nv <- npoints(vertices)
    if (nv <= 1) {
        m <- matrix(FALSE, nv, nv)
    } else if (!missing(m)) {
        if (!is.matrix(m) && !inherits(m, c("lgCMatrix", "lgTMatrix"))) 
            stop("m should be a matrix or sparse matrix")
        stopifnot(isSymmetric(m))
        if (nrow(m) != vertices$n) 
            stop("dimensions of matrix m do not match number of vertices")
        if (any(diag(m))) {
            warning("diagonal entries of the matrix m should not be TRUE; ignored")
            diag(m) <- FALSE
        }
        sparse <- !is.matrix(m)
        ij <- which(m, arr.ind = TRUE)
        ij <- ij[ij[, 1] < ij[, 2], , drop = FALSE]
        from <- ij[, 1]
        to <- ij[, 2]
        rid <- m.rid[ij]
    } else {
        stopifnot(is.matrix(edges) && ncol(edges) == 2)
        if (any((edges%%1) != 0)) 
            stop("Entries of edges list should be integers")
        if (any(self <- (edges[, 1] == edges[, 2]))) {
            warning("edge list should not join a vertex to itself; ignored")
            edges <- edges[!self, , drop = FALSE]
        }
        np <- npoints(vertices)
        if (any(edges > np)) 
            stop("index out-of-bounds in edges list")
        from <- edges[, 1]
        to <- edges[, 2]
        if (!sparse) {
            m <- matrix(FALSE, np, np)
            m[edges] <- TRUE
        } else
            m <- sparseMatrix(i = from, j = to, x = TRUE, dims = c(np, np))
        m <- m | t(m)
    }
    xx <- vertices$x
    yy <- vertices$y
    lines <- psp(xx[from], yy[from], xx[to], yy[to], window = vertices$window, check = FALSE)
    lines$rid <- rid
        
    out <- list(vertices = vertices, m = m, lines = lines, from = from, to = to, sparse = sparse, window = vertices$window)
    class(out) <- c("linnet", class(out))
    if (sparse)
        return(out)
    out$n <- n <- nrow(m)
    d <- Matrix(0, n, n)
    diag(d) <- 0
    d[m] <- pairdist(vertices)[m]
    out$dpath <- dpath <- dist2dpath(d)
    if (warn && any(is.infinite(dpath))) 
        warning("Network is not connected", call. = FALSE)
    out$circumradius <- circumradius(out)
    return(out)
}



# want a vectorized function that tells the value of the variable at the closest measured point
linfun.factory <- function(data, var) {
    myvar <- var
    
    function(x, y, seg, tp) {
        # create data objects
        out <- rep(NA, length(seg))
    
        # find the nearest quadrature knot (or point location) on the same segment,
        # and copy the covariate value observed there
        for (i in 1:length(seg)) {
            tmp.indx <- data$seg == seg[i]
            
            # if no quadrature knots are on this segment, then return NA
            if (any(is.na(tmp.indx))) {
              
              stop("NAs found while resolving functions on stream segments")
            }
            if (any(tmp.indx)) {
                out[i] <- data[[myvar]][tmp.indx][which.min(abs(data$tp[tmp.indx] - tp[i]))]
            }
        }
        
        # return the results
        out
    }
}





# want a vectorized function that tells the value of the variable at the closest measured point
lambda.cox.fun.factory <- function(object, covars, coef, lookup, mean, covariance) {
    datafun.list <- covars
    beta <- coef
    p <- length(datafun.list)
    varnames <- names(datafun.list)
    lookup.frame <- lookup
    mu <- mean
    Sigma <- covariance
    
    # generate the random effects
    u <- mvrnorm(1, mu=mu, Sigma=Sigma)
    
    function(x, y, seg, tp) {
        dat <- list()
        if (inherits(object, 'terms'))
            dat[[names(attr(object, 'dataClasses'))[1]]] <- 1
        
        for (i in seq_len(p))
            dat[[varnames[i]]] <- datafun.list[[i]](x, y, seg, tp)
        
        dat <- as.data.frame(dat)
        
        # the model.matrix function seems not to handle NAs.
        eta <- rep(NA, nrow(dat))
        indx <- apply(dat, 1, function(z) !any(is.na(z)))
        dat <- dat[indx,, drop=FALSE]
        
        # current.na.action <- options('na.action')
        options(na.action = na.pass)
        X <- model.matrix(object, dat)
        options(na.action = na.omit)
        
        # calculate the fixed-effect part of the output
        eta[indx] <- as.vector(X %*% coef)
        

        
        # identify which random effect goes with which entry
        # branches with negative id get no random effect
        brange <- rep(NA, max(lookup$br))
        valid <- unique(lookup$br)
        valid <- sort(valid[valid > 0])
        brange[valid] <- 1:length(valid)
        
        # match segment ids to branch ids
        indx <- match(seg, lookup$seg)
        braw <- lookup$br[indx]
        bflag <- braw > 0
        indx <- cbind((1:length(seg))[bflag], brange[braw[bflag]])
        
        # populate the random design matrix
        S <- Matrix(0, nrow=length(x), ncol=length(u))
        S[indx] <- 1
                
        # return the results
        eta <- eta + as.vector(S %*% u)
        exp(eta)
    }
}








# want a vectorized function that tells the value of the variable at the closest measured point
lambda.fun.factory <- function(object, covars, coef) {
    datafun.list <- covars
    beta <- coef
    p <- length(datafun.list)
    varnames <- names(datafun.list)
    
    
    function(x, y, seg, tp) {
        dat <- list()
        if (inherits(object, 'terms'))
            dat[[names(attr(object, 'dataClasses'))[1]]] <- 1
        
        for (i in seq_len(p))
            dat[[varnames[i]]] <- datafun.list[[i]](x, y, seg, tp)
        
        dat <- as.data.frame(dat)
        
        out <- rep(NA, nrow(dat))
        indx <- apply(dat, 1, function(z) !any(is.na(z)))
        dat <- dat[indx,, drop=FALSE]
        
        # current.na.action <- options('na.action')
        options(na.action = na.pass)
        X <- model.matrix(object, dat)
        options(na.action = na.omit)
        
        # return the results
        out[indx] <- as.vector(X %*% coef)
        exp(out)
    }
}





# extract the endpoints of the stream segments
# ends <- sapply(edges@lines, function(x) c(head(x@Lines[[1]]@coords, 1), tail(x@Lines[[1]]@coords, 1))) %>% t
# ends <- cbind(ends, edges@data$rid, edges@data$netID)
# ends <- as.data.frame(rbind(ends[,c(1, 2, 5, 6)], ends[,c(3, 4, 5, 6)]))
# colnames(ends) <- c('x', 'y', 'rid', 'netID')

# fix some errors in the data
# ends$netID[ends$netID == 30] <- 82
# ends$netID[ends$netID == 84] <- 82
# ends$netID[ends$netID == 87] <- 86
# ends$netID[ends$netID == 91] <- 90

# win <- owin(xrange=range(ends$x), yrange=range(ends$y))

# create a ppp object holding the vertices
# pts <- ends %>% group_by(x, y) %>% distinct
# pts$pid <- 1:nrow(pts)
# pts$rid <- NULL
# ends <- right_join(pts, ends, by=c('x', 'y', 'netID'))
# verts <- ppp(x=pts$x, y=pts$y, window=win)


# create the linnet object defining the network
# m <- Matrix(FALSE, nrow(pts), nrow(pts))
# m.rid <- Matrix(0, nrow(pts), nrow(pts))

# identify which reaches touch each other
# for (pp in unique(ends$pid)) {
#     cat(paste(pp, '\n'))
#     tangent <- ends$rid[ends$pid == pp]
#     connected <- ends$pid[ends$rid %in% tangent & ends$pid != pp]
#     m[pp, connected] <- m[connected, pp] <- TRUE
#     
#     m.rid[pp, ends$pid[ends$rid %in% tangent & ends$pid != pp]] <- ends$rid[ends$rid %in% tangent & ends$pid != pp]
# }

# create a linear network object using sparse matrices
# net <- MyLinnet(vertices=verts, m=m, m.rid=m.rid, sparse=TRUE, warn=FALSE)

# this is the key that matches rid to seg
# rid.lookup <- data.frame(rid=net$lines$rid, seg=1:net$lines$n)
# indx <- match(rid.lookup$rid, edges@data$rid)
# rid.lookup$netID <- edges@data$netID[indx]

# add the segment ID as used by spatstat functions to the Brown trout object
# bt <- X.cox
# bt$tp <- 1 - bt$ratio
# indx <- match(bt$rid, rid.lookup$rid)
# bt$seg <- rid.lookup$seg[indx]

# bt.lpp <- lpp(bt, net)

# add tp and seg to the data and change the names of x and y coordinates
# wyom.q <- X.q
# wyom.q$tp <- 1 - wyom.q$ratio
# wyom.q <- inner_join(wyom.q, rid.lookup, by='rid')
# colnames(wyom.q)[colnames(wyom.q) == 'NEAR_X'] <- 'x'
# colnames(wyom.q)[colnames(wyom.q) == 'NEAR_Y'] <- 'y'
# qpoints <- lpp(X=wyom.q, L=net)

# seg.q <- data.frame(seg = tail(wyom.q$seg, nrow(X.q)))
# seg.q <- seg.q %>% group_by(seg) %>% summarize(count = n())
# seg.q.count <- rep(0, nrow(rid.lookup))
# seg.q.count[seg.q$seg] <- seg.q$count
# rid.lookup$count <- seg.q.count


# rid.lookup$wt <- with(inner_join(edges@data, rid.lookup), Length / count)
# wyom.q$wt <- (left_join(wyom.q, rid.lookup))$wt
# wdum <- rep(epsilon, nrow(X.cox))

# wrap the data into functions of the type that linfun expects
# f.MS_Hist <- function(x, y, seg, tp) generic.linfun(x, y, seg, tp, streamdata, 'MS_Hist')
# f.log_MS_Hist <- function(x, y, seg, tp) log(generic.linfun(x, y, seg, tp, streamdata, 'MS_Hist'))
# f.S1_93_11 <- function(x, y, seg, tp) generic.linfun(x, y, seg, tp, streamdata, 'S1_93_11')
# f.SLOPE <- function(x, y, seg, tp) generic.linfun(x, y, seg, tp, streamdata, 'SLOPE')

# create linfuns for the covariates
# lf.MS_Hist <- linfun(f.MS_Hist, net)
# lf.log_MS_Hist <- linfun(f.log_MS_Hist, net)
# lf.S1_93_11 <- linfun(f.S1_93_11, net)
# lf.SLOPE <- linfun(f.SLOPE, net)

# plot(lf.MS_Hist)
# plot(lf.log_MS_Hist)
# plot(lf.S1_93_11)
# plot(lf.SLOPE)
