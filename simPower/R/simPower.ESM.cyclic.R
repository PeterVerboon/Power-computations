simPower.ESM.cyclic <-
function(nbeep, nday, nsubj, rho, b0=0, b1, b2, alpha , sd.subj = 1, sd.day = 0, maxiter=1000, ntest=1) {   
  
  
  alpha <- 0.05/ntest
  
  ntot <- nbeep*nsubj*nday
  res <- matrix(data=0,nrow=maxiter, ncol=8)
  colnames(res) <- c("Amplitude","Phase", "Predictor","Power_pred","Power_b1", "Power_C", "Power_S", "Random_Intercept")

  beepnr <- rep(seq(1:nbeep),(nsubj*nday))
   cor=matrix(c( 1,  rho , rho, 1 ), ncol=2, byrow=TRUE);
  
    
  for (i in 1:maxiter) {  
    
    dat1 <- data.frame(dat1 <- mvrnorm(n=ntot,mu = c(0,0), Sigma=cor)) ; 
    names(dat1) <- c("y","x")
    subjnr <- sort(rep(seq(1:nsubj),(nbeep*nday)))
    daynr <- rep(sort((rep(seq(1:nday),nbeep))),nsubj)
    dat1 <- data.frame(cbind(subjnr, daynr, beepnr,dat1))
    
  
    # add random effect across subjects
    
    subjnr <- seq(1:nsubj)
    esub <- rnorm(nsubj,0,sd.subj)
    dat2 <- data.frame(cbind(subjnr,esub ))
    dat <- merge(dat1, dat2, by = "subjnr")
    dat$y <- dat$y + dat$esub  

    # add random effect across days
    
    daynr <- seq(1:nday)
    eday <- rnorm(nday,0,sd.day)
    dat2 <- data.frame(cbind(daynr,eday ))
    dat <- merge(dat, dat2, by = "daynr")
    dat$y <- dat$y + dat$eday 
    
    # add cyclic effect within days
    
    dat$y <- (b0 + b1*cos((2*pi/nbeep)*(dat$beepnr-b2))) + dat$y
    
    
    # analysis
    
    dat$subj <- factor(x=dat$subjnr)
    dat$cvar <- cos((2*pi/nbeep)*dat$beepnr)
    dat$svar <- sin((2*pi/nbeep)*dat$beepnr)
    
    fit <- lmer(y ~ cvar + svar + x + (1 |subj), data = dat, REML = FALSE);         
 
    a0 <- fixef(fit)[1]
    a1 <- fixef(fit)[2]
    a2 <- fixef(fit)[3]
    b3 <- fixef(fit)[4]
    
    par <- cycpar(a1,a2,P = nbeep)   
    
    res[i,1] <- par[1]
    res[i,2] <- par[2]
    res[i,3] <- b3
    res[i,4] <- summary(fit)$coefficients[4,5]  < alpha
    res[i,5] <- par[1] > (b1-1)
    res[i,6] <- summary(fit)$coefficients[2,5]  < alpha
    res[i,7] <- summary(fit)$coefficients[3,5]  < alpha
    res[i,8]<- sqrt( unlist(summary(fit)$varcor))
    
  }
  
  res <- res[order(res[,1]),]
  return(res)
  
}
