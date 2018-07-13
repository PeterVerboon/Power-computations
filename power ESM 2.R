
require(MASS);
require(userfriendlyscience);
require(lme4)
require(lmerTest);
require(dplyr)

options(digits=3);




# START FUNCTION
# Simulate Power for simple ESM design
# effect size (rho) refers to the correlation between two interval variables 
#
# number of beeps (nbeep), number of days (nday) and number of subjects (nsubj)
# level of autocorrelation (ar) can be added
# random variance across subjects (sd.subj) and days (sd.day) can be adapted
# percentage of missings for beeps (pMisBeep) and days (pMisDay) can be added
# alpha level can be corrected for number of tests (ntest) by Bonferroni correctie 
# anayses are with random intercept across subjects only (randomDay = FALSE) or across subjects and days
#
# functions used: getBiCop
#
# author: P. Verboon, 
# date: october, 2016 adapted july, 2018


simPower.ESM <- function(nbeep, nday, nsubj, rho, ar, sd.subj = 1, sd.day = 0, 
                         maxiter=1000, ntest=1, randomDay = FALSE, estAR=TRUE, pMisDay = 0, pMisBeep = 0) 
  {   
  
  
  alpha <- 0.05/ntest   # correction for multiple tests
  
  ntot <- nbeep*nsubj*nday
  res <- matrix(data=0,nrow=maxiter, ncol=5)
  colnames(res) <- c("regression coefficient","Power", "auto_correlation", "sd_Subject", "sd_Day")

  beepnr <- rep(seq(1:nbeep),(nsubj*nday))

  # Initiate the Progress bar
  pb <- txtProgressBar(min = 0, max = maxiter, style = 3)
  
  # loop over replications
    
  for (i in 1:maxiter) {  
    
    # construct variable y with amount (ar) of auto-regression
    
    y1 <- arima.sim(list(ar = ar), n = ntot, sd=1)
    
    # construct second variable correlated (rho) with y
    
    dat1 <- as.data.frame(getBiCop(n=ntot, rho=rho, x=y1))
    
    # build dataframe
    
    colnames(dat1) <- c("y","x1")
    subjnr <- sort(rep(seq(1:nsubj),(nbeep*nday)))
    daynr <- rep(sort((rep(seq(1:nday),nbeep))),nsubj)
    dat1 <- data.frame(cbind(subjnr, daynr, beepnr,dat1))

   
    # add random effect across subjects
    
    subjnr <- seq(1:nsubj)
    esub <- rnorm(nsubj,0,sd.subj)
    dat3 <- data.frame(cbind(subjnr,esub ))
    dat <- merge(dat1, dat3, by = "subjnr")
    dat$y <- dat$y + dat$esub  

    # add random effect and missings across days
    
    daynr <- seq(1:nday)
    eday <- rnorm(nday,0,sd.day)
    misday <- rep(0, nday)
    misday[runif(nday,0,1) < pMisDay ] <- NA
    dat2 <- data.frame(cbind(daynr,eday,misday ))
    dat <- merge(dat, dat2, by = "daynr")
    dat$y <- dat$y + dat$eday + dat$misday
 

    # add missings at beep level
    
    misbeep <- rep(0, dim(dat)[1])
    misbeep[runif(dim(dat)[1],0,1) < pMisBeep ] <- NA
    dat$y <- dat$y + misbeep
    
    
    # make lagged variable for auto-correlation term in model
    
    dat <- arrange(dat,subjnr,daynr,beepnr)
    dat <- LagESM(dat, subjnr="subjnr",daynr="daynr",beepnr="beepnr", lagn=1, varnames= "y")
    
    # analysis
    
    if (estAR == FALSE)  {
      if (randomDay == TRUE) 
      {fit <-  lmer(y ~  x1  + ( 1 |subjnr) + (1 | daynr), data=dat) }
      else
      {fit <- lmer(y ~  x1  + ( 1 |subjnr), data=dat)}
    }
    
    if (estAR == TRUE)  {
      if (randomDay == TRUE) 
      {fit <-  lmer(y ~  x1 + yL1 + ( 1 |subjnr) + (1 | daynr), data=dat) }
      else
      {fit <- lmer(y ~  x1 + yL1 + ( 1 |subjnr), data=dat)}
    }
    
    
    res[i,1] <- summary(fit)$coefficients[2,1]
    res[i,2] <- summary(fit)$coefficients[2,5]  < alpha
    if (estAR == TRUE) {res[i,3] <- summary(fit)$coefficients[3,1] }
    res[i,c(4,5)]<- sqrt( unlist(summary(fit)$varcor))
    
    setTxtProgressBar(pb, i)
  }
  
  cat("\n")
  
  res <- res[order(res[,1]),]
  return(res)
  
}  # END FUNCTION




# Call to function:

res <- simPower.ESM(nbeep = 8, 
                    nday = 14, 
                    nsubj = 25, 
                    sd.subj = 1, 
                    sd.day = 0, 
                    maxiter=10, 
                    rho = 0.15, 
                    ar = 0.05,
                    ntest = 4, 
                    randomDay = FALSE, 
                    estAR = TRUE,
                    pMisDay = .40, 
                    pMisBeep = .50)


# show mean coefficient and power:

apply(res,2,mean)



# Show 95% coverage interval :

l <- dim(res)[1]
CI <- res[c((0.025*l),(0.975*l))]
names(CI)  <- c("Lower limit", "    Upper limit")
CI

# Plot density of coefficients from simulation:

plot(density(res[,1]))

