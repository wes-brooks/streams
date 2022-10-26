############################
# brook trout
coefs <- matrix(NA, 0, 6)
Kfun <- matrix(NA, 513, 0)
ll.envelope <- vector() 

for (i in (1:200)[-193]) {
    load(paste("~/Dropbox/streams/output/brook-cox/brook-cox-results", i, "RData", sep='.'))
    Kfun <- cbind(Kfun, inhom$est)
    ll.envelope <- c(ll.envelope, cox$neg.loglik)
    coefs <- rbind(coefs, cox$beta)
}

brook.cox <- list()
brook.cox$coefs <- coefs
brook.cox$Kfun <- Kfun
brook.cox$ll.envelope <- ll.envelope



############################
# brown trout
coefs <- matrix(NA, 0, 6)
Kfun <- matrix(NA, 513, 0)
ll.envelope <- vector() 

for (i in 1:200) {
    load(paste("brown-cox-results", i, "RData", sep='.'))
    Kfun <- cbind(Kfun, inhom$est)
    ll.envelope <- c(ll.envelope, cox$neg.loglik)
    coefs <- rbind(coefs, cox$beta)
}

brown.cox <- list()
brown.cox$coefs <- coefs
brown.cox$Kfun <- Kfun
brown.cox$ll.envelope <- ll.envelope
