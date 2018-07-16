
# Compute power for simple moderation model

require(userfriendlyscience)
require(MASS)
require(lm.beta)

N=500
rho1 <- .5
rep <- 500
e <- sqrt(3)

out <- matrix(0,nrow=rep,ncol=8)
varxz <- rep(0,rep)
vary <- rep(0,rep)

for (i in 1:rep) {
  
  a <- simDataSet(N, varNames = c("x","z"),  means=c(0,0), sds = c(1, 1), specifiedCorrelations = list(c('x','z',rho1)), ranges = NULL, silent=TRUE, seed = NULL)
  x <- a$x; z <- a$z
  xz <- x*z
  varxz[i] <- var(xz)

  error <-  rnorm(N,0,sqrt(e))   # the amount of error determines the ES

   y <- sqrt(.5)*x + sqrt(.3)*z + sqrt(.2)*xz + error

   vary[i] <- var(y)

   res <- lm( y ~ x + z + xz) 
   res <- lm.beta(res)
   sig <- summary(res)$coefficients[4,5] < .05
   rsq <- summary(res)$r.squared
 
  out[i,] <- c(rsq, res$coef[-1], res$standardized.coefficients[-1], sig)
}


apply(out,2,mean)
mean(varxz)
mean(vary)

# R squared

(mean(vary) - e)/(mean(vary) )


### ES as function of variance of error

sqrt(.5)/sqrt(10)
sqrt(.3)/sqrt(10)
sqrt(.2)/sqrt(10)

## with correlated predictors

sqrt(.5) * (1 / sqrt(mean(vary)))
sqrt(.3) * (1 / sqrt(mean(vary)))
sqrt(.2) * (sqrt(mean(varxz)) / sqrt(mean(vary)))



## expected variance of product xz

1 + cov((x**2),(z**2)) - (cov(x,z)**2)

1 + (2*rho1**2) - (cov(x,z)**2)


str(res)



           