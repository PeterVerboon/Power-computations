
# Compute power for simple moderation model

require(userfriendlyscience)
require(MASS)
require(lm.beta)

N=400
rho1 <- .7
rho2 <- .7
rep <- 500

out <- matrix(0,nrow=rep,ncol=7)
varxz <- rep(0,rep)
vary <- rep(0,rep)

for (i in 1:rep) {
  
  a <- simDataSet(N, varNames = c("x","z"),  corr = c(rho1,rho2), seed = NULL, silent = TRUE)
  x <- a$x; z <- a$z
  xz <- x*z
  varxz[i] <- var(xz)


e <-  rnorm(N,0,3)   # the amount of error determines the ES

y <- sqrt(.5)*x + sqrt(.3)*z + sqrt(.2)*xz + e

vary[i] <- var(y)

# res <- lm( y ~ x + z + xz) 
# res <- lm.beta(res)
# sig <- summary(res)$coefficients[4,5] < .05
#   
# 
# out[i,] <- c(res$coef[-1], res$standardized.coefficients[-1], sig)
}


apply(out,2,mean)
mean(varxz)
mean(vary)


### ES as function of variance of error

sqrt(.5)/sqrt(10)
sqrt(.3)/sqrt(10)
sqrt(.2)/sqrt(10)


## expected variance of product xz

1 + cov((x**2),(z**2)) - (cov(x,z)**2)

1 + (2*rho1**2) - (cov(x,z)**2)

