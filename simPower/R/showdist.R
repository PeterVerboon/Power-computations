showdist <-
function(dat, var, plot = TRUE) {
  options(digits = 3)
  b <- dat[order(dat[,var]),]
  cat("mean  : ", mean(b[,var])," \n")
  cat("median: ", median(b[,var])," \n")
  cat("95% coverage interval: ", b[nrow(b)*c(0.05, 0.95),var]," \n")
  cat("99% coverage interval: ", b[nrow(b)*c(0.01, 0.99),var])
  if (plot) {plot(density(b[,var]), main=var, xlab= "")}
}
