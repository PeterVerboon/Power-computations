
require(MASS);
require('userfriendlyscience');
require(lme4);
require(nlme)

options(digits=3);


# Generate data according tot autocorrelation model

require(DataCombine)

dat <- slide(dat, Var = "e", GroupVar = "subjnr", NewVar = "yvar2", slideBy = -1)
dat$yvar2 <- dat$yvar2 + dat$e + dat$esub
dat$yvar2[is.na(dat$yvar2)] <- dat$e[is.na(dat$yvar2)]


# START FUNCTION
# Simulate Power for simple ESM design
# based on effect size (rho), number of beeps (nbeep) and number of subjects (nsubj)
# alpha level can be corrected for number of tests (ntest) by Bonferroni correctie 
# P. Verboon, october, 2016


simPower.ESM <- function(nbeep, nsubj, rho=.2, maxiter=1000, ntest=1) {   
  
  
  alpha <- 0.05/ntest
  
  ntot <- nbeep*nsubj
  res <- matrix(data=0,nrow=maxiter, ncol=2)
  colnames(res) <- c("regression coefficient","Power")

  beeps <- rep(seq(1:nbeep),nsubj)
  cor=matrix(c( 1,  rho , rho, 1 ), ncol=2, byrow=TRUE);
  
    
  for (i in 1:maxiter) {  
    
    dat2 <- data.frame(dat2 <- mvrnorm(n=ntot,mu = c(0,0), Sigma=cor)) ; 
    names(dat2) <- c("y","x1")
    subjnr <- sort(rep(seq(1:nsubj),nbeep))
    dat1 <- data.frame(cbind(subjnr, beeps,dat2))
    
  
    # add random effect across subjects
    
    subjnr <- seq(1:50)
    esub <- rnorm(50,0,1)
    dat2 <- data.frame(cbind(subjnr,esub ))
    
    dat <- merge(dat1, dat2, by = "subjnr")
    dat$y <- dat$y + dat$esub  
    
    fit <-  lme(y ~  x1 , random = ~1 |subjnr, data=dat,na.action = (na.omit));
    
    res[i,1] <- fixef(fit)[2]
    res[i,2] <- summary(fit)$tTable[,"p-value"][2]  < alpha
    
  }
  
  res <- res[order(res[,1]),]
  return(res)
  
}  # END FUNCTION




# Call to function:

res <- simPower.ESM(nbeep = 10, nsubj = 65, maxiter=500, rho =.15, ntest=3)


# show mean coefficient and power:

apply(res,2,mean)



# Show 95% coverage interval :

l <- length(res)*.5
CI <- res[c((0.025*l),(0.975*l))]
names(CI)  <- c("Lower limit", "    Upper limit")
CI

# Plot density of coefficients from simulation:

plot(density(res[,1]))

