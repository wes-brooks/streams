datagen.rpoisppOnLines2 <- function (lambda, L, lmax = NULL, ..., effort=NULL, nzw.rids=NULL, correction.factor=NULL, check = TRUE, loglambda=FALSE) {
  # args <- list(...)
  # effort <- args$effort
  # 
  # # recover the effort function if it was supplied. If not, create a dummy function that returns zeroes.
  # if (is.null(effort)) {
  #   effort <- list(function(x, y, seg, tp) {rep(0, length(x))})
  # } else if (is.function(lambda) || is.im(lambda)) {
  #   effort <- list(effort)
  # }
  # m.effort <- length(effort)
  # 
  # argtype.eff <- if (all(unlist(lapply(effort, is.im)))) 
  #   "im"
  # else if (all(unlist(lapply(effort, is.function)))) 
  #   "function"
  # else stop(paste(sQuote("effort"), "must be a function, an image,", 
  #                 "a list of functions, or a list of images"))
  
  if (is.numeric(lambda)) 
    return(datagen.runifpoisppOnLines(lambda, L, ...))
  if (is.function(lambda) || is.im(lambda)) 
    lambda <- list(lambda)
  m <- length(lambda)
  argtype <- if (all(unlist(lapply(lambda, is.im)))) 
    "im"
  else if (all(unlist(lapply(lambda, is.function)))) 
    "function"
  else stop(paste(sQuote("lambda"), "must be a numeric vector, a function, an image,", 
                  "a list of functions, or a list of images"))
  
  # recover the effort function if it was supplied. If not, create a dummy function that returns zeroes.
  if (is.null(effort)) {
    if (argtype == "function") {
      effort <- sapply(1:m, function(i) return(function(x, y, seg, tp) {rep(0, length(x))}), simplify=FALSE)
    } else if (argtype == "im") {
      effort <- list(lambda)
      for (i in 1:m)
        effort[[i]]$v <- matrix(0, nrow(effort[[i]]$v), ncol(effort[[i]]$v))
    }
  } else if (is.function(effort) || is.im(effort)) {
    effort <- list(effort)
  }
  m.effort <- length(effort)
  
  argtype.eff <- if (all(unlist(lapply(effort, is.im)))) 
    "im"
  else if (all(unlist(lapply(effort, is.function)))) 
    "function"
  else stop(paste(sQuote("effort"), "must be a function, an image,", 
                  "a list of functions, or a list of images"))
  
  
  if (argtype != argtype.eff)
    stop(paste(sQuote("lambda"), "and", sQuote("effort"), "must be of the same type (image or function)"))
  
  if (m != m.effort)
    stop(paste(sQuote("lambda"), "and", sQuote("effort"), "must be of the same length"))

  
  if (argtype == "im") {
    for (j in seq_len(m)) {
      lamj <- lambda[[j]]
      if (!(lamj$type %in% c("real", "integer"))) 
        stop("lambda must be numeric-valued or integer-valued")
      lrange <- range(lamj)
      if (any(is.infinite(lrange))) 
        stop("Infinite pixel values not permitted")
      if (lrange[1] < 0) 
        stop("Negative pixel values not permitted")
    }
  }
  if (!is.null(lmax)) {
    stopifnot(is.numeric(lmax))
    if (length(lmax) != m) {
        if (length(lmax) == 1) {
          lmax <- rep.int(lmax, m)
        }
        else stop("Length of lmax does not match length of lambda")
    }
  }
  else {
    lmax <- numeric(m)
    for (j in seq_len(m)) {
      lamj <- lambda[[j]]
      effj <- effort[[j]]
      if (is.function(lamj)) {
        X <- pointsOnLines2(L, np = 10000, ...)
        lambdaX <- lamj(X$data$x, X$data$y, X$data$seg, X$data$tp, ...) *  effj(X$data$x, X$data$y, X$data$seg, X$data$tp, ...)
        lmax[j] <- max(lambdaX, na.rm = TRUE)
      }
      else if (is.im(lamj)) 
        lmax[j] <- max(lamj * effj)
    }
    if (!all(is.finite(lmax))) 
      stop("Infinite values of lambda obtained")
    if (any(lmax < 0)) 
      stop("Negative upper bound for lambda obtained")
    names(lmax) <- names(lambda)
  }
  Y <- datagen.runifpoisppOnLines(lmax, L, nzw.rids=nzw.rids, correction.factor=correction.factor, ...)
  n <- nrow(Y)
  if (n == 0) 
      return(Y)
  if (m == 1) {
    lambda <- lambda[[1]]
    effort <- effort[[1]]
    markindex <- 1
    if (is.function(lambda)) {
      if (loglambda) {
        lambdaY <- exp(lambda(Y$x, Y$y, seg=Y$seg, tp=Y$tp, ...) + effort(Y$x, Y$y, seg=Y$seg, tp=Y$tp, logresult=TRUE, ...))
      } else {
        lambdaY <- lambda(Y$x, Y$y, seg=Y$seg, tp=Y$tp, ...) * effort(Y$x, Y$y, seg=Y$seg, tp=Y$tp, ...)
      }
    } else lambdaY <- safelookup(lambda, as.ppp(Y, W = as.owin(L))) * safelookup(effort, as.ppp(Y, W = as.owin(L)))
  }
  else {
    lambdaY <- numeric(n)
    markindex <- as.integer(Y$marks)
    for (j in seq_len(m)) {
      lamj <- lambda[[j]]
      effj <- effort[[j]]
      jrows <- (markindex == j)
      Yj <- Y[jrows, , drop = FALSE]
      if (is.function(lamj)) 
        if (loglambda)
          lambdaY[jrows] <- exp(lamj(Yj$x, Yj$y, ...) + effj(Yj$x, Yj$y, logresult=TRUE, ...))
        else
          lambdaY[jrows] <- lamj(Yj$x, Yj$y, ...) * effj(Yj$x, Yj$y, ...) 
      else lambdaY[jrows] <- safelookup(lamj, as.ppp(Yj, W = as.owin(L))) * safelookup(effj, as.ppp(Yj, W = as.owin(L)))
    }
  }
    
    
  naprop <- sum(L$lines$Length[Y$seg[is.na(lambdaY)]]) / sum(L$lines$Length) / 2
  lambdaY[is.na(lambdaY)] <- 0
  pY <- lambdaY / lmax[markindex]
  if (check) {
    if (any(pY < 0)) 
      warning("Negative values of lambda obtained")
    if (any(pY > 1)) 
      warning("lmax is not an upper bound for lambda")
  }
  retain <- (runif(n) < pY)
  Y <- Y[retain, , drop = FALSE]
  return(Y)
}
