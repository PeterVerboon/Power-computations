
require(MASS);
require('userfriendlyscience');
require(lme4)
require(lmerTest);
#require(nlme)

options(digits=5);




# START FUNCTION
# Simulate Power for  ESM design with cyclic terms
# based on effect size (rho), number of beeps (nbeep), number of days (nday) and number of subjects (nsubj)
# random variance across subjects (sd.subj) and days (sd.day) can be adapted
# alpha level can be corrected for number of tests (ntest) by Bonferroni correctie 
# P. Verboon, october, 2016


simPower.cyc.ESM <- function(nbeep, nday, nsubj, rho, b0=0, b1, b2, alpha.cond , sd.subj = 1, sd.day = 0, maxiter=1000, ntest=1) {   
  
  
  alpha <- 0.05/ntest
  
  ntot <- nbeep*nsubj*nday
  res <- matrix(data=0,nrow=maxiter, ncol=6)
  colnames(res) <- c("Amplitude","Phase", "Power_b1", "Power_C", "Power_S", "Random_Intercept")

  beepnr <- rep(seq(1:nbeep),(nsubj*nday))
   cor=matrix(c( 1,  rho , rho, 1 ), ncol=2, byrow=TRUE);
  
    
  for (i in 1:maxiter) {  
    
    dat1 <- data.frame(dat1 <- mvrnorm(n=ntot,mu = c(0,0), Sigma=cor)) ; 
    names(dat1) <- c("y","x1")
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
    
    fit <- lmer(y ~ cvar + svar + (1 |subj), data = dat, REML = FALSE);         
 
    a0 <- fixef(fit)[1]
    a1 <- fixef(fit)[2]
    a2 <- fixef(fit)[3]
    
    par <- cycpar(a1,a2,P = nbeep)   
    
    res[i,1] <- par[1]
    res[i,2] <- par[2]
    res[i,3] <- par[1] > alpha.cond
    res[i,4] <- summary(fit)$coefficients[2,5]  < alpha
    res[i,5] <- summary(fit)$coefficients[3,5]  < alpha
    res[i,6]<- sqrt( unlist(summary(fit)$varcor))
    
  }
  
  res <- res[order(res[,1]),]
  return(res)
  
}  # END FUNCTION




# Call to function:

res <- simPower.cyc.ESM(nbeep = 5, nday = 7, nsubj = 50, b1=.3, b2=3, alpha.cond=.05, sd.subj = 1, sd.day = 0, maxiter=10, rho =.00, ntest=1)


# show mean coefficient and power:

apply(res,2,mean)



# Show 95% coverage interval :

l <- dim(res)[1]
CI <- res[c((0.01*l),(0.025*l),(0.05*l),(0.950*l),(0.975*l),(0.990*l)),1]
names(CI)  <- c("1.0% Perc", "2.5% Perc","5.0% Perc","95.0% Perc","97.5% Perc","99.0% Perc")
CI

# Plot density of coefficients from simulation:

plot(density(res[,1]))

# loop over simulations
#  i : 25, 50, 75 (subjects)
#  j : 1, 7 (days)
#  k : 5, 10 (beeps)
#  es:  0.1, 0.2, 0.3
# total = 3 x 2 x 2 x 3 = 36

 # START SIMULATION 

  out1 <- matrix(data=NA, nrow=36, ncol=10)
  out2 <- matrix(data=NA, nrow=36, ncol=6)
  colnames(out1) <- c("nsubj","ndays", "nbeeps", "effectsize", "Amplitude", "Phase","Power_b1","Power_C", "Power_S", "Random_Intercept" )
  colnames(out2) <- c("1.0% Perc", "2.5% Perc","5.0% Perc","95.0% Perc","97.5% Perc","99.0% Perc")
    
    maxiter <- 1000
    rownum <- 0
    b2 <- 3
    for (i in c(25,50,75)) { 
     print(paste0("nsubj: ",i))
         for (j in c(1,7)) { 
        print(paste0("ndays: ",j))
              for (k in c(5,10)) { 
                print(paste0("nbeeps: ",k))
                for (es in c(0.1,0.2,0.3)) { 
                  print(paste0("effectsize: ",es))
                  rownum <- rownum + 1

     out1[rownum,c(1:4)] <- c(i,j,k,es)             
     a <- H0.95perc[rownum,2]
       
    res <- simPower.cyc.ESM(nbeep = k, nday = j, nsubj = i, b1= es, b2=b2, alpha.cond = a, sd.subj = 1, sd.day = 0, maxiter=maxiter, rho =.00, ntest=1)
 
    out1[rownum,c(5:10)] <- apply(res,2,mean)  
    
    out2[rownum,] <- CI <- res[c((0.01*maxiter),(0.025*maxiter),(0.05*maxiter),(0.950*maxiter),(0.975*maxiter),(0.990*maxiter)),1]
       
                 } # (es)
             } # (k)
         }  # (j)
} # end loop (i)

    # dit ook nog draaien under H0 too obtain significance levels 1%, 5% en 10% for b1
    # dit geeft een 12 (condities) x 3 (levels) tabel
    # resultaten heirvan gebruiken in de simulatie: vergelijk elke b1 met level om power te verkrijgen
    
    
   sink("output1 simulation 17okt.txt")
   out1
   sink()
  
   sink("output2 simulation 17okt.txt")
   out2
   sink()
 
     
  plot(res[,1])
  density(res[,1])
  plot(density(res[,1]))
  bias <- out1[,5] - out1[,4]
  
   out.es03 <- out1[out1[,4] == 0.3,]
   