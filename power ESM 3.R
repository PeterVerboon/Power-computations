
require(MASS);
require(userfriendlyscience);
require(lme4)
require(lmerTest);   # necessary for p-values
require(dplyr)
require(MuMIn);   # for Rsquare from lmer

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
# date:  july, 2018


simPower.ESM <- function(nbeep = 5,
                         nday = 7,
                         samSize = c(10,20,30,40,50), 
                         rho = 0.3,
                         arlevel = c(0, 0.2, 0.8), 
                         errlevel = c(1,3,9),
                         betas = c(.5, .3, .2),
                         alpha = 0.05,
                         sd.subj = 1, 
                         sd.day = 0, 
                         maxiter=500, 
                         ntest=1, 
                         randomDay = FALSE, 
                         estAR=TRUE, 
                         pMisDay = 0, 
                         pMisBeep = 0) 
  {   
  
  b1 <- betas[1]
  b2 <- betas[2]
  b3 <- betas[3]
  
  numrow <- length(errlevel)*length(arlevel)*length(samSize)
  result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= 10))
  colnames(result) <- c("N","e","ar","x","Power_x","interaction","Power_int", "auto_correlation", "Power_ar", "R-square")
  result[,"e"] <- rep(sort(rep(errlevel,length(arlevel))),length(samSize))
  result[,"ar"] <- rep(arlevel,length(errlevel)*length(samSize))
  result[,"N"] <- sort(rep(samSize,length(errlevel)*length(arlevel)))
  
  cor=matrix(c( 1, rho,rho, 1), ncol=2, byrow=TRUE);
  
  # Initiate the Progress bar
  pb <- txtProgressBar(min = 0, max = maxiter*numrow, style = 3)
  j <- 0
  
  for (nsubj in samSize) {
      for (e in errlevel) {
        for (ar in arlevel) {
  
  ntot <- nbeep*nsubj*nday
  out <- matrix(data=0,nrow=maxiter, ncol=7)
  colnames(out) <- c("x","power_x","interaction","Power_int", "auto_correlation", "Power_ar", "R-square")

  
   # loop over replications
    
   for (i in 1:maxiter) {  
    
    j <- j + 1
     
    # construct data eith given correlation structure
    
    dat1 <- data.frame(dat1 <- mvrnorm(n=ntot,mu = c(0,0), Sigma=cor)) ; 
    
    colnames(dat1) <- c("x","z")
    dat1$xz <- dat1$x*dat1$z
    
    beepnr <- rep(seq(1:nbeep),(nsubj*nday))
    subjnr <- sort(rep(seq(1:nsubj),(nbeep*nday)))
    daynr <- rep(sort((rep(seq(1:nday),nbeep))),nsubj)
    dat1 <- data.frame(cbind(subjnr, daynr, beepnr,dat1))
    
    dat1$y <- sqrt(b1)*dat1$x + sqrt(b2)*dat1$z + sqrt(b3)*dat1$xz  
    
    # construct lagged variabel for data construction, using ESM design
    
    dat1$yL1 <- lag(dat1$y)
    
    error <-  rnorm(ntot,0,sqrt(e))   # the amount of error determines the ES
    
    # construct dependent variable
    
    dat2 <- dat1
    dat2[(is.na(dat2$yL1)),"yL1"] <- 0
    dat2$y2 <- sqrt(1-ar)*dat2$y + sqrt(ar)*dat2$yL1 + error
    
  #   dat2 <- LagESM(dat2, subjnr="subjnr",daynr="daynr",beepnr="beepnr", lagn=1, varnames= "y2")
  #   var(dat2$y2) ; sd(dat2$y2)  
  #   var(dat2$x) ;  sd(dat2$x)
  #   var(dat1$y) ; sd(dat1$y)
  #   var(dat2$yL1, na.rm = TRUE); sd(dat2$yL1,na.rm = TRUE)
  #   var(error)
  #   dat2 <- LagESM(dat2, subjnr="subjnr",daynr="daynr",beepnr="beepnr", lagn=1, varnames= "y2")
  #   summary(lm(y2 ~ x + z + xz + y2L1, data = dat2))
  # #   
  # cor(dat1[,-c(1:3)], use = "complete.obs")
  # cor(dat2[,-c(1:3)], use = "complete.obs")
  # cov(dat2[,-c(1:3)], use = "complete.obs")
    
  
   
    # add random effect across subjects
    
    subjnr <- seq(1:nsubj)
    esub <- rnorm(nsubj,0,sd.subj)
    dat3 <- data.frame(cbind(subjnr,esub ))
    dat <- merge(dat2, dat3, by = "subjnr")
    dat$y2 <- dat$y2 + dat$esub  
    
 #   var(dat$y2) ; sd(dat$y2)
 #   cor(dat[,-c(1:3)], use = "complete.obs")
    
    # add random effect and missings across days
    
    daynr <- seq(1:nday)
    eday <- rnorm(nday,0,sd.day)
    misday <- rep(0, nday)
    misday[runif(nday,0,1) < pMisDay ] <- NA
    dat2 <- data.frame(cbind(daynr,eday,misday ))
    dat <- merge(dat, dat2, by = "daynr")
    dat$y2 <- dat$y2 + dat$eday + dat$misday
 

    # add missings at beep level
    
    misbeep <- rep(0, dim(dat)[1])
    misbeep[runif(dim(dat)[1],0,1) < pMisBeep ] <- NA
    dat$y2 <- dat$y2 + misbeep
    
    
    # make lagged variable for auto-correlation term in model
    
    dat <- arrange(dat,subjnr,daynr,beepnr)
    dat <- LagESM(dat, subjnr="subjnr",daynr="daynr",beepnr="beepnr", lagn=1, varnames= "y2")
    
    # analysis
    
    if (estAR == FALSE)  {
      if (randomDay == TRUE) 
      {fit <-  lmer(y2 ~  x + z + xz +  ( 1 |subjnr) + (1 | daynr), data=dat) }
      else
      {fit <- lmer(y2 ~  x + z + xz +  ( 1 |subjnr), data=dat)}
    }
    
    if (estAR == TRUE)  {
      if (randomDay == TRUE) 
      {fit <-  lmer(y2 ~  x + z + xz +  y2L1 + ( 1 |subjnr) + (1 | daynr), data=dat) }
      else
      {fit <- lmer(y2 ~  x + z + xz + y2L1 + ( 1 |subjnr), data=dat)}
    }
    
    stfit <- lm.beta.lmer(fit)
    
    out[i,1] <- stfit[1]
    out[i,2] <- summary(fit)$coefficients[2,5]  < alpha
    out[i,3] <- stfit[3]
    out[i,4] <- summary(fit)$coefficients[4,5]  < alpha
    if (estAR == TRUE) {
      out[i,5] <- stfit[4]
      out[i,6] <- summary(fit)$coefficients[5,5]  < alpha
      }
    out[i,7]<- r.squaredGLMM(fit)[1]
    
    setTxtProgressBar(pb, j)
    
    
   }    # end loop over replications
  
  result[((result$e == e) & (result$ar == ar) & result$N == nsubj),c(4:10)] <- apply(out,2,mean)
 # print(c(nsubj, e, ar))
  
  
        }  # end loop over ar level 
      }  # end loop over error level 
  }  # end loop over sample sizes
  
  cat("\n")
  return(result)
  
}  # END FUNCTION




# Call to function:

res <- simPower.ESM(nbeep = 5, 
                    nday = 14, 
                    samSize = c(10,20,30,40), 
                    rho = 0.3,
                    arlevel = c(0.3), 
                    errlevel = c(3,9),
                    betas = c(.1, .3, .2),
                    alpha = 0.05,
                    sd.subj = 1, 
                    sd.day = 0, 
                    maxiter=200, 
                    ntest=1, 
                    randomDay = FALSE, 
                    estAR=TRUE, 
                    pMisDay = 0.2, 
                    pMisBeep = 0.1) 



save(res, file= "result_E39AR03_alpha05.Rdata")
 load("result_E39AR03_alpha05.Rdata")

E1 <- res
E39 <- res
res1 <- rbind(E1,E39)

## Plot the results

require(ggplot2)

res1 <- res[res$ar == 0.3,]
res1$effectSize <- ordered(res1$e, levels=c(1,3,9), labels= c("0.31","0.23","0.15"))


p <- ggplot(data=res1, aes(y=Power_x, x=N, group=effectSize, colour=effectSize)) + geom_point() + geom_line()
p <- p + geom_hline(yintercept=0.80, linetype="dashed", color = "red")
p <- p + geom_hline(yintercept=0.90, linetype="dashed", color = "blue")
p <- p + coord_cartesian(ylim=c(0.1, 1.0)) + scale_y_continuous(breaks=seq(0.10, 1, 0.10))
p <- p + coord_cartesian(xlim=c(9, 41)) + scale_x_continuous(breaks=seq(10, 40, 10))
p <- p + ggtitle("Power for interaction term for alpha = 0.05")

p


