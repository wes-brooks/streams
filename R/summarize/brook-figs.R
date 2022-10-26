load("output/brook-obj.rdata")
obj.brook <- obj
npar <- nrow(obj.brook$cox$hessian)

# use linear combination matrix to transform the coefficients and covariances
beta <- lincomb.brook %*% obj.brook$cox$beta
covmat <- lincomb.brook %*% solve(obj.brook$cox$obj$hessian[(npar-6):(npar-1), (npar-6):(npar-1)]) %*% t(lincomb.brook)


lincomb %*% obj.brook$dwpr$fit$coef
sqrt(diag(lincomb.brook %*% summary(obj.brook$dwpr$fit$internal$glmfit)$cov.unscaled %*% t(lincomb.brook))) %>% round(3)



########################################################
# plot the inhomogeneous K envelope for the DWPR model
pdf("figures/brook-dwpr-inhom-K.pdf", width=5, height=5)
  brook$Kenvelope <- apply(brook$Kfun, 1, function(z) quantile(z, c(0.05, 0.95)))
  yy <- range(c(brook$Kenvelope, obj$dwpr$inhomogeneousK$est))
  plot(x=obj.brook$dwpr$inhomogeneousK$r, y=brook$Kenvelope[1,], type='l', lty=2, bty='n', main="", xlab='Stream distance (m)', ylab='K', ylim=yy)
  par(new=TRUE)
  plot(x=obj.brook$dwpr$inhomogeneousK$r, y=brook$Kenvelope[2,], type='l', lty=2, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
  par(new=TRUE)
  plot(x=obj.brook$dwpr$inhomogeneousK$r, y=obj.brook$dwpr$inhomogeneousK$est, type='l', lty=1, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
dev.off()



######################################################
# plot the inhomogeneous K envelope for the Cox model
pdf("figures/brook-cox-inhom-K.pdf", width=5, height=5)
  brook.cox$Kenvelope <- apply(brook.cox$Kfun, 1, function(z) quantile(z, c(0.05, 0.95)))
  # yy <- range(brook.cox$Kenvelope) 
  yy <- range(c(brook.cox$Kenvelope, obj.brook$cox$inhomogeneousK$est))
  plot(x=obj.brook$cox$inhomogeneousK$r, y=brook.cox$Kenvelope[1,], type='l', lty=2, bty='n', main="", xlab='Stream distance (m)', ylab='K', ylim=yy)
  par(new=TRUE)
  plot(x=obj.brook$cox$inhomogeneousK$r, y=brook.cox$Kenvelope[2,], type='l', lty=2, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
  par(new=TRUE)
  plot(x=obj.brook$cox$inhomogeneousK$r, y=obj.brook$cox$inhomogeneousK$est, type='l', lty=1, bty='n', xaxt='n', yaxt='n', ann=FALSE, ylim=yy)
dev.off()


#######################################################
# Plot a histogram of DWPR parametric bootstrap log-likelihood resamples
pdf("figures/brook-dwpr-nll-boot.pdf", width=5, height=5)
  hist(-brook$ll.envelope, breaks=12, xlab="negative log-likelihood", main="")
  lines(x=rep(-obj.brook$dwpr$fit$maxlogpl, 2), y=c(0, 48), lty=2, col='blue', lwd=2)
  text(x=17867, y=42, labels=paste("p=", round(2 * (1-mean(-obj.brook$dwpr$fit$maxlogpl < -brook$ll.envelope)), 2)))
dev.off()


#######################################################
# Plot a histogram of Cox parametric bootstrap log-likelihood resamples
pdf("figures/brook-cox-nll-boot.pdf", height=5, width=5)
  hist(brook.cox$ll.envelope, breaks=12, xlab="negative log-likelihood", main="")
  lines(x=rep(obj.brook$cox$neg.loglik, 2), y=c(0, 35), lty=2, col='blue', lwd=2)
  text(x=-19300, y=25, labels=paste("p=", round(mean(obj.brook$cox$neg.loglik > brook.cox$ll.envelope)), 2))
dev.off()