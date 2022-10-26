# this bisquare function is vectorized to accelerate computation
bisquare <- function(dist, bw) {
  outer(dist, bw, FUN=function(d, h) ifelse(d>h, 0, (1 - (d/h)^2)^2))
}