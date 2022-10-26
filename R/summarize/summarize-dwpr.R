############################
# brook trout
coefs <- matrix(NA, 0, 6)
Kfun <- matrix(NA, 513, 0)
ll.envelope <- vector() 

for (i in 1:200) {
    load(paste("~/hdrive/streams/output/brook-dwpr/brook-dwpr-results", i, "RData", sep='.'))
    Kfun <- cbind(Kfun, inhom$est)
    ll.envelope <- c(ll.envelope, loglik)
    coefs <- rbind(coefs, beta)
}

brook <- list()
brook$coefs <- coefs
brook$Kfun <- Kfun
brook$ll.envelope <- ll.envelope



##########################################
# Brown trout
coefs <- matrix(NA, 0, 6)
Kfun <- matrix(NA, 513, 0)
ll.envelope <- vector() 

for (i in 1:200) {
    load(paste("~/hdrive/streams/output/brown-dwpr/brown-dwpr-results", i, "RData", sep='.'))
    Kfun <- cbind(Kfun, inhom$est)
    ll.envelope <- c(ll.envelope, loglik)
    coefs <- rbind(coefs, beta)
}

brown <- list()
brown$coefs <- coefs
brown$Kfun <- Kfun
brown$ll.envelope <- ll.envelope
