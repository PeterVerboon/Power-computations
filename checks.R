
ve <- 9    # 1, 3, 9
r <- .3    # 0, 0.3, 0.5, 0.8



vy <- (.5) + (.3) + (1+r**2)*(.2) + ve

b1 <- sqrt(.5)/sqrt(vy)
b2 <- sqrt(.3)/sqrt(vy)
b3 <- sqrt(.2)/sqrt(vy)

rsq <- ((.5) + (.3) + (1+r**2)*(.2))/vy


## extra tests


res1 <- res

# R squared

(mean(vary) - mean(err))/(mean(vary) )


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

