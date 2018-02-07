
options(digits=3);


require(MASS);
require('userfriendlyscience');
# require(lme4);
require(nlme)
require(DataCombine)



# START FUNCTION
# Simulate Power for simple ESM design
# based on effect size (ES), number of beeps (nbeep) , number of subjects (nsubj), Sd of random effect across subjects (sdsub) and autocorrelation (ar)
# alpha level can be corrected for number of tests (ntest) by Bonferroni correctie 
# P. Verboon, january, 2017


simPower.ESM.autocor <- function(nbeep, nsubj, ES=.2, maxiter=1000, ar = 0.5, sdsub = 0, ntest=1) {   
  
  alpha <- 0.05/ntest
  
  ntot <- nbeep*nsubj
  res <- matrix(data=0,nrow=maxiter, ncol=5)
  colnames(res) <- c("regression coefficient","autocorrelation"," Power Beta","Power AR1", "StdDev Random")

  beeps <- rep(seq(1:nbeep),nsubj)
 
    
  for (i in 1:maxiter) {  
    
    yvar <- arima.sim(list(order = c(1,0,0), ar = ar), n = ntot, sd = 1)         # An autocorrelation simulation
    yvar <- scale(yvar)                                                          # standardize yvar
    dat0 <- data.frame(getBiCop(n=ntot,rho=ES,x=yvar) )                          # add variable with given correlation
    subjnr <- sort(rep(seq(1:nsubj),nbeep))
    dat1 <- cbind(subjnr, beeps,dat0)
    names(dat1) <- c("subjnr","beepnr","y","x1")

    dat2 <- slide(dat1, Var = "y", GroupVar = c("subjnr"), NewVar = "ylag1", slideBy = -1, reminder = FALSE)  # with more days this is only correct for first subjects 
    ## dat2[(dat2[,"beepnr"]==1),"yvar2"]  <-- NA                                       # add NA at first beep of the day
    
      
  
    # add random effect across subjects
    
    subjnr <- seq(1:nsubj)
    esub <- rnorm(nsubj,0,sdsub)
    dat3 <- data.frame(cbind(subjnr,esub))
    
    dat <- merge(dat2, dat3, by = "subjnr")
    dat$y <- dat$y + dat$esub  
    
    fit <-  lme(y ~  x1 + ylag1 , random = ~1 |subjnr, data=dat,na.action = (na.omit));
    
    res[i,1] <- fixef(fit)[2]
    res[i,2] <- fixef(fit)[3]
    res[i,3] <- summary(fit)$tTable[,"p-value"][2]  < alpha
    res[i,4] <- summary(fit)$tTable[,"p-value"][3]  < alpha
    res[i,5] <- as.numeric(VarCorr(fit)[1,2] )
    
  }
  
  res <- res[order(res[,1]),]
  return(res)
  
}  # END FUNCTION





# Call to function:

res <- simPower.ESM.autocor (nbeep = 10, nsubj = 50, ES = .20, maxiter=400, ar =.50, sdsub = 1, ntest=3)


# show mean coefficient and power:

apply(res,2,mean)



# Show 95% coverage interval :

l <- length(res)*.5
CI <- res[c((0.025*l),(0.975*l))]
names(CI)  <- c("Lower limit", "    Upper limit")
CI

# Plot density of coefficients from simulation:

plot(density(res[,1]))
plot(density(res[,2]))

