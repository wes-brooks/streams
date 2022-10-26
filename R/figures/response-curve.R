# Generate the matrix that converts the coefficients from orthogonal polynomials to raw polynomials
lincomb <- function(obj, dat) {
  XX <- rbind(obj$coxdata, obj$quaddata)
  
  res <- diag(ncol(as.matrix(dat)))
  quad.coefs <- c('S1_93_11')
  lin.coefs <- c('SLOPE', 'logMS_Hist', 'offset')
  for (cc in quad.coefs) {
    locname <- paste("poly(", cc, ", 2)", sep='')
    lindx <- which(names(coef(obj$dwpr)) == paste(locname, '1', sep=''))
    quadx <- which(names(coef(obj$dwpr)) == paste(locname, '2', sep=''))
    
    pars <- attr(dat[[cc]], 'coefs')
    
    # Orthogonal polynomials are a complicated transformation.
    # Add some terms to the intercept
    res[1, lindx] <- -pars$alpha[1] / sqrt(pars$norm2[3])
    res[1, quadx] <- -mean(XX[[cc]]^2) / sqrt(pars$norm2[4]) +
      mean(XX[[cc]]) * sum(dat[[cc]][,1]*XX[[cc]]^2) / sqrt(pars$norm2[4]) / sqrt(pars$norm2[3])
    
    # Raw linear coefficient is a linear combination of the orthogonal polynomial linear and quadratic coefficients
    res[lindx, lindx] <- 1 / sqrt(pars$norm2[3])
    res[lindx, quadx] <- -sum(dat[[cc]][,1]*XX[[cc]]^2) / sqrt(pars$norm2[4]) / sqrt(pars$norm2[3])
    
    # Raw quadratic coefficient is just scaled from the orthogonal polynomial quadratic coefficient
    res[quadx, quadx] <- 1 / sqrt(pars$norm2[4])
    
    # lincomb[1, quadx] <- -(mean(XX[[cc]]) * lincomb[lindx, quadx] + mean(XX[[cc]]^2) * lincomb[quadx, quadx])
  }
  
  # back-tranfsorm the coefficients that were only centered and scaled:
  for (cc in lin.coefs) {
    locname <- paste("poly(", cc, ", 1)", sep='')
    lindx <- which(names(coef(obj$dwpr)) == locname)
    
    pars <- attr(dat[[cc]], 'coefs')
    
    # Orthogonal polynomials are a complicated transformation.
    # Add some terms to the intercept
    res[1, lindx] <- -pars$alpha[1] / sqrt(pars$norm2[3])
    
    # Raw linear coefficient is a linear combination of the orthogonal polynomial linear and quadratic coefficients
    res[lindx, lindx] <- 1 / sqrt(pars$norm2[3])
  }
  
  # return the result
  res
}



# load the data from both models
load('~/Desktop/brown-obj.rdata')
obj.brown <- obj
XX.brown <- rbind(obj.brown$coxdata, obj.brown$quaddata)
X.brown <- model.matrix(~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1), data=XX.brown)
dat <- data.frame(XX.brown$logMS_Hist)
dat <- within(dat, {
  logMS_Hist <- poly(XX.brown$logMS_Hist, 1);
  S1_93_11 <- poly(XX.brown$S1_93_11, 2);
  SLOPE <- poly(XX.brown$SLOPE, 1);
  offset <- poly(XX.brown$offset, 1);
})
lincomb.brown <- lincomb(obj.brown, dat)


# load the data from both models
load('~/Desktop/brook-obj.rdata')
obj.brook <- obj
XX.brook <- rbind(obj.brook$coxdata, obj.brook$quaddata)
X.brook <- model.matrix(~ poly(logMS_Hist, 1) + poly(S1_93_11, 2) + poly(SLOPE, 1) + poly(offset, 1), data=XX.brook)
dat <- data.frame(XX.brook$logMS_Hist)
dat <- within(dat, {
  logMS_Hist <- poly(XX.brook$logMS_Hist, 1);
  S1_93_11 <- poly(XX.brook$S1_93_11, 2);
  SLOPE <- poly(XX.brook$SLOPE, 1);
  offset <- poly(XX.brook$offset, 1);
})
lincomb.brook <- lincomb(obj.brook, dat)




# Need to back out of the orthogonal polynomial transformations in order to make plots on the scale of measured covariate values
spaced <- function(x) seq(min(x), max(x), len=400)
mm.brook <- cbind(1, spaced(XX.brook$logMS_Hist), spaced(XX.brook$S1_93_11), spaced(XX.brook$S1_93_11)^2,
            spaced(XX.brook$SLOPE), spaced(XX.brook$offset))
mm.brown <- cbind(1, spaced(XX.brown$logMS_Hist), spaced(XX.brown$S1_93_11), spaced(XX.brown$S1_93_11)^2,
            spaced(XX.brown$SLOPE), spaced(XX.brown$offset))


# Get the median effect of each covariate
median.covars.brook <- c(1, median(XX.brook$logMS_Hist), median(XX.brook$S1_93_11),
                         median(XX.brook$S1_93_11^2), median(XX.brook$SLOPE), median(XX.brook$offset))
median.effect.brook <- apply(X.brook, 2, median) * obj.brook$cox$beta


median.covars.brown <- c(1, median(XX.brown$logMS_Hist), median(XX.brown$S1_93_11),
                         median(XX.brown$S1_93_11^2), median(XX.brown$SLOPE), median(XX.brown$offset))
median.effect.brown <- apply(X.brown, 2, median) * obj.brown$cox$beta


# labels and legends may differ by what covariate we're plotting
xcap <- c(NA,
          expression(paste('Natural logrithm of flow (', ft^3, s^{-1}, ')')),
          expression(paste('Summer stream temperature (', degree, 'C)')),
          expression(paste('Summer stream temperature (', degree, 'C)')),
          'Slope')

xvars <- c(NA, 'logMS_Hist', 'S1_93_11', 'S1_93_11', 'SLOPE')
legendloc <- c(NA, 'topleft', 'topright', 'topright', 'topright')


layout(matrix(1:3, 1, 3))

for (j in list(2, c(3,4), 5)) {
  y.brook <- rowSums(mm.brook %*% lincomb.brook[,j] %*% as.matrix(obj.brook$cox$beta[j])) + sum(median.effect.brook[1:5][-j])
  y.brown <- rowSums(mm.brown %*% lincomb.brown[,j] %*% as.matrix(obj.brown$cox$beta[j])) + sum(median.effect.brown[1:5][-j])
  
  
  
  
  # calculate the confidence intervals for the profile plots
  tmp.brook <- mm.brook
  tmp.brook[,-j] <- rep(median.covars.brook[-j], each=nrow(tmp.brook))
  stderr.brook <- diag(tmp.brook %*% lincomb.brook %*% summary(obj.brook$dwpr$fit$internal$glmfit)$cov.unscaled %*%
                         t(tmp.brook %*% lincomb.brook)) #+ exp(obj.brook$cox$log.sigma.prior / 2)
  
  tmp.brown <- mm.brown
  tmp.brown[,-j] <- rep(median.covars.brown[-j], each=nrow(tmp.brown))
  stderr.brown <- diag(tmp.brown %*% lincomb.brown %*% summary(obj.brown$dwpr$fit$internal$glmfit)$cov.unscaled %*%
                         t(tmp.brown %*% lincomb.brown)) #+ exp(obj.brown$cox$log.sigma.prior / 2)
  
  
  stds <- 2
  
  # get the limits of the plotting area
  xx <- range(c(mm.brown[,j[1]], mm.brook[,j[1]]))
  yy <- range(c(y.brook + stds * stderr.brook, y.brook - stds * stderr.brook,
                y.brown - stds * stderr.brown, y.brown - stds * stderr.brown))
  
  
  # make the plotting window
  plot.new()
  plot.window(xlim=xx, ylim=yy)
  
  
  # Begin with confidence intervals so that they lay on the bottom of the finished plot
  polygon(x=c(mm.brook[,j[1]], rev(mm.brook[,j[1]])), y=c(y.brook + stds*stderr.brook, rev(y.brook - stds*stderr.brook)), border=NA, col='grey80')
  polygon(x=c(mm.brown[,j[1]], rev(mm.brown[,j[1]])), y=c(y.brown + stds*stderr.brown, rev(y.brown - stds*stderr.brown)), border=NA, col='grey80')
  
  par(new=TRUE)
  
  # plot a profile of the marginal response, holding other covariates at their median
  plot(mm.brook[,j[1]], y.brook, xlim=xx, ylim=yy, #log='y',
       type='l', bty='n', xlab=xcap[j[1]], ylab='Natural logarithm of intensity', col='blue', lwd=2)#, cex.axis=1.6, cex.lab=1.6)
  par(new=TRUE)
  plot(mm.brown[,j[1]], y.brown, xlim=xx, ylim=yy, #log='y',
       type='l', bty='n', ann=FALSE, xaxt='n', yaxt='n', col='brown', lwd=2)
  
  # draw an axis for the density
  dens <- density(XX.brown[[xvars[j[1]]]])
  # limdens <- max(dens$y) / 0.75
  # axis(4, at=seq(yy[1], yy[2], len=4), labels=round(seq(0, limdens, len=4), 1))
  
  # add a density curve
  dens$y <- dens$y * 0.75 * diff(yy) / diff(range(dens$y))
  dens$y <- dens$y - (min(dens$y) - yy[1])
  lines(dens, lty=2)
  
  legend(x=legendloc[j[1]], legend=c('brook trout', 'brown trout', 'density of points'), #cex=1.6,
         bty='n', border='white', col=c('blue', 'brown', 'black'), lwd=c(2, 2, 1), lty=c(1, 1, 2), bg='white', box.col='white', y.intersp=2.5)

  
}

