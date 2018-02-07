
# returns a data frame of two variables which correlate with a population correlation of rho
# If desired, one of both variables can be fixed to an existing variable by specifying x

getBiCop <- function(n, rho, mar.fun=rnorm, x = NULL, ...) {

    if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(n, ...)}

    if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
  
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  X2 <- X2 <- mar.fun(n,mean(x),sd(x)) 
  X <- cbind(X1,X2)
  
  # induce correlation (does not change X1)
  df <- X %*% C
  
  ## if desired: check results
  ## all.equal(X1,X[,1])
  ## cor(X)
  
  return(df)
}


aa <- getBiCop(n=50,rho=.5,x=a)
cor(aa)

