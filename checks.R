
ve <- 9    # 1, 3, 9
r <- .3    # 0, 0.3, 0.5, 0.8
ar <- .8    # 0, 0.2, 0.8

b1 <- .5
b2 <- .3
b3 <- .2
b4 <- ar

vy1 <- (b1) + (b2) + (1+r**2)*(b3) + 2*sqrt(b1)*sqrt(b2)*r                                # with interaction, no AR
vy <- (b1) + (b2) + (1+r**2)*(b3) + 2*sqrt(b1)*sqrt(b2)*r + (ar)*sqrt(vy1)                # with interaction, with AR
vy <- (b1) + (b2) + (1+r**2)*(b3) + 2*sqrt(b1)*sqrt(b2)*r + (ar)*sqrt(vy1) + (ve)      # with interaction, with AR, with error


be1 <- (b1)/sqrt(vy)
be2 <- (b2)/sqrt(vy)
be3 <- (b3)*(sqrt(1+r**2)/sqrt(vy))
be4 <- (b4)*(sqrt(vy1)/sqrt(vy))


# R squared

(vy - ve)/vy

(mean(vary) - mean(err))/(mean(vary) )


### ES as function of variance of error

sqrt(.5)/sqrt(2.5)
sqrt(.3)/sqrt(2.5)
sqrt(.2)/sqrt(2.5)

## with correlated predictors

sqrt(.5) * (1 / sqrt(mean(vary)))
sqrt(.3) * (1 / sqrt(mean(vary)))
sqrt(.2) * (sqrt(mean(varxz)) / sqrt(mean(vary)))



## expected variance of product xz

1 + cov((x**2),(z**2)) - (cov(x,z)**2)

1 + (2*rho1**2) - (cov(x,z)**2)


str(res)

