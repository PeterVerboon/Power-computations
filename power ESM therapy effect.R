
require(MASS);
require('userfriendlyscience');
require(lme4)
require(lmerTest);
#require(nlme)

options(digits=3);




# START FUNCTION
# Simulate Power for simple repeated measures design in (single case) therapy research
# based on Cohen's d effect size (ES), number of days before therapy (nbefore), number of days total (nday) and number of subjects (nsubj)
# random variance across subjects (sd.subj) and days (sd.day) can be adapted
# alpha level can be corrected for number of tests (ntest) by Bonferroni correctie 
# P. Verboon, february, 2017


simPower.ESM.th <- function(nbefore=6, nday, nsubj, ES, sd.subj = 1, maxiter=1000, ntest=1) {   
  
  
  alpha <- 0.05/ntest
  
  ntot <- nsubj*nday
  res <- matrix(data=0,nrow=maxiter, ncol=3)
  colnames(res) <- c("regression coefficient","Power", "sd_Subject")


  for (i in 1:maxiter) {  
    
    dat1 <- data.frame(dat1 <- mvrnorm(n=ntot,mu = c(0,0), Sigma=cor)) ; 
    ybefore <- rnorm(n=nbefore*nsubj, sd=1, mean=0)
    yafter  <- rnorm(n=(nday-nbefore)*nsubj, sd=1, mean=ES)
    y <- c(ybefore, yafter)
    
 
    daynr <- (rep(seq(1:nday),nsubj))
    subjnr <- sort(rep(seq(1:nsubj),nday))
    phase <- rep(c(rep(0,nbefore),rep(1,(nday-nbefore))),nsubj)
    
  
     dat1 <- data.frame(cbind(subjnr, daynr,phase,dat1))
     dat2 <- dat1[order(dat1$phase),]
     dat2["y"] <- c(ybefore, yafter)
    
     # res1 <- computeEffectSize_d(dat2$phase, dat2$y) 
     # res1$es
     
    # add random effect across subjects
    
    subjnr <- seq(1:nsubj)
    esub <- rnorm(nsubj,0,sd.subj)
    dat3 <- data.frame(cbind(subjnr,esub ))
    dat <- merge(dat2, dat3, by = "subjnr")
    dat$y <- dat$y + dat$esub  

    
    
    # analysis
    
    fit <-  lmer(y ~  phase + ( 1 |subjnr), data=dat);
 
    res[i,1] <- summary(fit)$coefficients[2,1]
    res[i,2] <- summary(fit)$coefficients[2,5]  < alpha
    res[i,3]<- sqrt( unlist(summary(fit)$varcor))
    
  }
  
  res <- res[order(res[,1]),]
  return(res)
  
}  # END FUNCTION




# Call to function:

res <- simPower.ESM.th(nbefore=6, nday = 30, nsubj = 50, sd.subj = 1,  maxiter=300, ES =.20, ntest=2)


# show mean coefficient and power:

apply(res,2,mean)



# Show 95% coverage interval :

l <- dim(res)[1]
CI <- res[c((0.025*l),(0.975*l))]
names(CI)  <- c("Lower limit", "    Upper limit")
CI

# Plot density of coefficients from simulation:

plot(density(res[,1]))

