
# Compute power for simple moderation model

require(userfriendlyscience)
require(MASS)
require(lm.beta)

options(digits = 3)
 

simPower.moderation <- function(samSize = c(50,100,150,200,250,300), 
                                errlevel = c(1,3,9), 
                                rholevel = c(0,.3,.5,.8),
                                betas = c(.5, .3, .2),
                                rep = 1000) 
      {   
    numrow <- length(errlevel)*length(rholevel)*length(samSize)
    result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= 8))
    colnames(result) <- c("N","e","rho", "rsq", "b1","b2","b3","power")
    result[,"e"] <- rep(sort(rep(errlevel,length(rholevel))),length(samSize))
    result[,"rho"] <- rep(rholevel,length(errlevel)*length(samSize))
    result[,"N"] <- sort(rep(samSize,length(errlevel)*length(rholevel)))
    
    b1 <- betas[1]
    b2 <- betas[2]
    b3 <- betas[3]
    
     for (N in samSize)
     {
       for (e in errlevel) 
         { 
         for (rho in rholevel)
           {
           sigma = matrix(c(1,rho,rho,1), nrow=2, ncol=2)
           
           out <- matrix(0,nrow=rep,ncol=5)
           
           for (i in 1:rep) 
             {
  
      a <- as.data.frame(mvrnorm(n = N, Sigma=sigma, mu=c(0,0)))
      colnames(a)  <- c("x","z")
      a$xz <- a$x*a$z

      error <-  rnorm(N,0,sqrt(e))   # the amount of error determines the ES

      a$y <- sqrt(b1)*a$x + sqrt(b2)*a$z + sqrt(b3)*a$xz + error

      res <- lm( y ~ x + z + xz, data=a) 
      res <- lm.beta(res)
      sig <- summary(res)$coefficients[4,5] < .05
      rsq <- summary(res)$r.squared
 
      out[i,] <- c(rsq,  res$standardized.coefficients[-1], sig)
  
            }   # loop over replications
           
           
           result[((result$e == e) & (result$rho == rho) & result$N == N),c(4:8)] <- apply(out,2,mean)
           print(c(N, e,rho))
           
          }     # loop over correlations between predictors
         }      # loop over error levels
       }        # loop over sample sizes
    
         return(result)

} # end function



res <- simPower.moderation(samSize = c(100,150,200,250,300),  errlevel = c(1,3,9), rholevel = c(0.3), betas = c(.5,.3,.2), rep=1000) 


res[(res$rho == 0.3),]

  
save(res, file= "result_moderatie.Rdata")
  
require(pander)
pander(res) 

require(ggplot2)


p <- ggplot(data=res, aes(y=power, x=N, group=e, colour=e)) + geom_point() + geom_line()
p <- p + geom_hline(yintercept=0.80, linetype="dashed", color = "red")


  
p
  
   ## extra tests
  
  


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



           