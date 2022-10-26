load("/home/z3509569/Desktop/brown-obj.rdata")
obj.brown <- obj
npar <- nrow(obj.brown$cox$hessian)

# use linear combination matrix to transform the coefficients and covariances
beta <- lincomb.brown %*% obj.brown$cox$beta
covmat <- lincomb.brown %*% solve(obj.brown$cox$hessian[(npar-6):(npar-1), (npar-6):(npar-1)]) %*% t(lincomb.brown)


lincomb.brown %*% obj.brown$dwpr$fit$coef
sqrt(diag(lincomb.brown %*% summary(obj.brown$dwpr$fit$internal$glmfit)$cov.unscaled %*% t(lincomb.brown))) %>% round(3)


# plot the inhomogeneous K envelope for the DWPR model
pdf("figures/brown-dwpr-inhom-K.pdf", width=5, height=5)
  brown$Kenvelope <- apply(brown$Kfun, 1, function(z) quantile(z, c(0.05, 0.95)))
  yy <- range(c(brown$Kenvelope, obj.brown$dwpr$inhomogeneousK$est))
  plot(x=obj.brown$dwpr$inhomogeneousK$r, y=brown$Kenvelope[1,], type='l', lty=2, bty='n', main="", xlab='Stream distance (m)', ylab='K', ylim=yy)
  par(new=TRUE)
  plot(x=obj.brown$dwpr$inhomogeneousK$r, y=brown$Kenvelope[2,], type='l', lty=2, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
  par(new=TRUE)
  plot(x=obj.brown$dwpr$inhomogeneousK$r, y=obj.brown$dwpr$inhomogeneousK$est, type='l', lty=1, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
dev.off()




# plot the inhomogeneous K envelope for the cox model
pdf("figures/brown-cox-inhom-K.pdf", width=5, height=5)
  brown.cox$Kenvelope <- apply(brown.cox$Kfun, 1, function(z) quantile(z, c(0.05, 0.95)))
  yy <- range(c(brown.cox$Kenvelope, obj.brown$cox$inhomogeneousK$est))
  plot(x=obj.brown$cox$inhomogeneousK$r, y=brown.cox$Kenvelope[1,], type='l', lty=2, bty='n', main="", xlab='Stream distance (m)', ylab='K', ylim=yy)
  par(new=TRUE)
  plot(x=obj.brown$cox$inhomogeneousK$r, y=brown.cox$Kenvelope[2,], type='l', lty=2, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
  par(new=TRUE)
  plot(x=obj.brown$cox$inhomogeneousK$r, y=obj.brown$cox$inhomogeneousK$est, type='l', lty=1, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
dev.off()




#######################################################
# Plot a histogram of DWPR parametric bootstrap log-likelihood resamples
pdf("figures/brown-dwpr-nll-boot.pdf", width=5, height=5)
  hist(-brown$ll.envelope, breaks=12, xlab="negative log-likelihood", main="")
  lines(x=rep(-obj.brown$dwpr$fit$maxlogpl, 2), y=c(0, 55), lty=2, col='blue', lwd=2)
  text(x=11700, y=45, labels=paste("p=", format(round(2 * (1-mean(-obj.brown$dwpr$fit$maxlogpl >= -brown$ll.envelope)), 2), nsmall=2)))
dev.off()


#######################################################
# Plot a histogram of Cox parametric bootstrap log-likelihood resamples
pdf("figures/brown-cox-nll-boot.pdf", height=5, width=5)
  hist(brown.cox$ll.envelope, breaks=12, xlab="negative log-likelihood", main="")
  lines(x=rep(obj.brown$cox$neg.loglik, 2), y=c(0, 35), lty=2, col='blue', lwd=2)
  text(x=-24100, y=25, labels=paste("p=", round(2 * (1-mean(obj.brown$cox$neg.loglik >= brown.cox$ll.envelope)), 2)))
dev.off()
