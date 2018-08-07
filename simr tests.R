
require(psych)
require(lme4)
require(lmerTest);   # necessary for p-values
require(dplyr)
require(MuMIn);




## Function for replications

rep.ESMsim <- function(nbeep = 5,
                       nday = 7,
                       samSize = c(10,20), 
                       rho = 0.1,
                       arlevel = c(0, 0.2, 0.8), 
                       errlevel = c(1,3,9),
                       betas = c(.5, .3, .2),
                       alpha = 0.05,
                       sd.subj = 1, 
                       sd.day = 0, 
                       maxiter=100, 
                       ntest=1, 
                       pMisDay = 0, 
                       pMisBeep = 0)
   {


numrow <- length(errlevel)*length(arlevel)*length(samSize)
result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= 10))
colnames(result) <- c("N","e","ar","x","Power_x","SD(y)","SD(x)", "auto_correlation", "Power_ar", "R-square")
result[,"e"] <- rep(sort(rep(errlevel,length(arlevel))),length(samSize))
result[,"ar"] <- rep(arlevel,length(errlevel)*length(samSize))
result[,"N"] <- sort(rep(samSize,length(errlevel)*length(arlevel)))


for (nsubj in samSize) {
  
  cat("nsubj:",nsubj)
  cat("\n")
  
  for (e in errlevel) {

    cat("error:",e)
    cat("\n")
    
    for (ar in arlevel) {

      cat("ar:",ar)
      cat("\n")
      
      ntot <- nbeep*nsubj*nday
      out <- matrix(data=0,nrow=maxiter, ncol=7)
      colnames(out) <- c("x","power_x","sd_y","null", "auto_correlation", "Power_ar", "R-square")
  
          for (i in 1:maxiter) {  
            
          #em <- matrix(c(1,1),ncol=2, nrow=nsubj, byrow = TRUE)
            
          dat <- sim.multi(n.obs=nsubj,
                           nvar=2,
                           nfact=2,
                           loading=c(1,1),
                           ntrials= nbeep*nday,
                           days=nday,
                           sigma=sd.subj, 
                           sin.i = 0, 
                           cos.i = 0,
                           phi.i = rho,
                           sigma.i = 1,
                           AR1 = ar, 
                           plot = F)
 
  colnames(dat) <- c("F1","F2","y","x","time","id")
  beepnr <- rep(seq(1:nbeep),(nsubj*nday))
  subjnr <- sort(rep(seq(1:nsubj),(nbeep*nday)))
  daynr <- rep(sort((rep(seq(1:nday),nbeep))),nsubj)
  dat1 <- data.frame(cbind(subjnr, daynr, beepnr,dat))
  
  dat1$y <- dat1$y + rnorm(ntot,0,sqrt(e)) 
 
  dat <- LagESM(dat1, subjnr="subjnr",daynr="daynr",beepnr="beepnr", lagn=1, varnames= "y")
  
  # add missings at beep level
  
  misbeep <- rep(0, dim(dat)[1])
  misbeep[runif(dim(dat)[1],0,1) < pMisBeep ] <- NA
  dat$y <- dat$y + misbeep
  
       
  fit <- lmer(y ~    x  + yL1 + ( 1 |subjnr), data=dat)
  
  stfit <- lm.beta.lmer(fit)
 
  out[i,1] <- stfit[1]
  out[i,2] <- summary(fit)$coefficients[2,5]  < alpha
  out[i,3] <- sd(dat$y,na.rm =TRUE)
  out[i,4] <- sd(dat$x, na.rm =TRUE)
    out[i,5] <- stfit[2]
    out[i,6] <- summary(fit)$coefficients[3,5]  < alpha
  out[i,7]<- r.squaredGLMM(fit)[1]
  
  
          } # replications
      
      
  result[((result$e == e) & (result$ar == ar) & result$N == nsubj),c(4:10)] <- apply(out,2,mean)
        
 
 
               }   # ar loop
           }  # end loop over error level 
   }   # nsubj loop

      cat("\n")
      return(result)

} # end function


res <- rep.ESMsim(nbeep=10, 
                  nday=7, 
                  samSize = c(10,25,40),
                  errlevel = c(1,3,9), 
                  arlevel = c(0, 0.2, 0.8), 
                  rho = .10,
                  sd.subj = 1,
                  pMisBeep = .40,
                  maxiter=200)

save(res, file= "sim_E139AR028_al05rho10.Rdata")
load("sim_E39AR28_alpha05.Rdata")

## Plot the results

require(ggplot2)

res <- res1
res2 <- res
res1 <- res
res2$condition <- NA
res1 <- rbind(res1,res)

res1$condition <- factor(1,levels=c("error = 1, ar = 0.0",
                                    "error = 3, ar = 0.0",
                                    "error = 9, ar = 0.0",
                                    "error = 1, ar = 0.2",
                                    "error = 3, ar = 0.2",
                                    "error = 9, ar = 0.2",
                                    "error = 1, ar = 0.8",
                                    "error = 3, ar = 0.8",
                                    "error = 9, ar = 0.8"))

res1[(res1$ar == .0) & (res1$e == 1),"condition"] <- "error = 1, ar = 0.0"
res1[(res1$ar == .0) & (res1$e == 3),"condition"] <- "error = 3, ar = 0.0"
res1[(res1$ar == .0) & (res1$e == 9),"condition"] <- "error = 9, ar = 0.0"
res1[(res1$ar == .2) & (res1$e == 1),"condition"] <- "error = 1, ar = 0.2"
res1[(res1$ar == .2) & (res1$e == 3),"condition"] <- "error = 3, ar = 0.2"
res1[(res1$ar == .2) & (res1$e == 9),"condition"] <- "error = 9, ar = 0.2"
res1[(res1$ar == .8) & (res1$e == 1),"condition"] <- "error = 1, ar = 0.8"
res1[(res1$ar == .8) & (res1$e == 3),"condition"] <- "error = 3, ar = 0.8"
res1[(res1$ar == .8) & (res1$e == 9),"condition"] <- "error = 9, ar = 0.8"



p <- ggplot(data=res1, aes(y=Power_x, x=N, group=condition, colour=condition)) + geom_point() + geom_line()
p <- p + geom_hline(yintercept=0.80, linetype="dashed", color = "red")
p <- p + geom_hline(yintercept=0.90, linetype="dashed", color = "blue")
p <- p + coord_cartesian(ylim=c(0.1, 1.0)) + scale_y_continuous(breaks=seq(0.10, 1, 0.10))
#p <- p + coord_cartesian(xlim=c(9, 41)) + scale_x_continuous(breaks=seq(10, 40, 10))
p <- p + ggtitle("Power ESM predictor for rho=0.10 and alpha = 0.05, pmis=.30 ")

p


acf(dat[,6])$acf[2]
acf(dat[,7])$acf[2]

summary(fit)
cor(dat, use="complete.obs")
apply(dat,2,mean)
apply(dat,2,sd)
statsBy(dat, group = "id",cors=T)



## simulate ML data

rho1 <- .3
rho2 <- .4
rho3 <- .6
eta <- c(.2,.5)

rwg=matrix(c( 1, rho1,rho2,
              rho1, 1, rho1,
              rho2, rho1, 1), ncol=3, byrow=TRUE);
rbg=matrix(c( 1, rho1,rho3,
              rho1, 1, rho1,
              rho3, rho1, 1), ncol=3, byrow=TRUE);

eta <- rep(.5,3)

a <- sim.multilevel(ncases=100, nvar=3, ngroups = 5, rwg=rwg, rbg=rbg, eta=eta)

cor(a$wg)
cor(a$bg)
a$xy
