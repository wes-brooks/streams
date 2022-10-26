datagen.runifpoisppOnLines <- function (lambda, L, nzw.rids=NULL, correction.factor=NULL, ...) {
  # L <- as.psp(L)
  # mu <- lambda * sum(lengths.psp(L))
  
  # if we have specified that not all lines have positive intensity then exclude the zero-intensity lines
  if (!is.null(nzw.rids)) {
    indx <- L$lines$rid %in% nzw.rids
  } else indx <- rep(TRUE, length(L$lines$Length))
  
  # specify a correction factor because some lines have nonzero intensity only on part of the lines.
  if (is.null(correction.factor)) {
    correction.factor <- data.frame(rid=L$lines$rid, scale=rep(1, length(L$lines$Length)))
  }
  
  scale.indx <- match(L$lines$rid, correction.factor$rid)
  mu <- lambda * sum(L$lines$Length[indx] * correction.factor$scale[scale.indx][indx])
  n <- rpois(rep.int(1, length(mu)), mu)
  if (length(n) > 1) 
      names(n) <- names(lambda)
  df <- datagen.runifpointOnLines(n, L, nzw.rids=nzw.rids, correction.factor=correction.factor, ...)
  return(df)
}