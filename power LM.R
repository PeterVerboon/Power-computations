
setwd("~/Onderzoek/Project Ellin")

setwd("~/Documents/Open Universiteit/Onderzoek/Project Ellin ")
getwd()

options(digits=3);


require(MASS)
require(lme4);
require(nlme);
library(reshape2)


# Fixed parameters in simulation 
npred=5
predNamesm = paste0("m", 1:npred)
predNamesw = paste0("w", 1:npred)
predNamesy = paste0("y", 1:npred)

## error correlation matrix ##

ID <- matrix(0,ncol=5,nrow=5)
diag(ID) <- 1

# contrasten 

T0vsT14 <- c(-.5, 1/8, 1/8, 1/8, 1/8)
T01vsT24 <- c(-.25, -.25, 1/3, 1/3, 1/3)
T0vsT34 <- c(-.5, 0, 0, .25, .25)


# expected means in repeated measures ##
pattern <-    c(0,  0.73, 0.93 , 0.93  , 0.93)     # treatment group, DV
pattern_m <-  c(0,  0.50, 0.50 , 0.50  , 0.50)     # treatment group, mediator 1
pattern_w <-  c(0,  0.00, 0.50 , 0.50  , 0.50)     # treatment group, mediator 2
pattern0 <-   c(0,  0.23, 0.23 , 0.23  , 0.23)     # control group, DV (waiting list effect)

# pattern of proportion of missing values (attrition)
mislev <- c(.10,.20,.30,.40)*1


#### START FUNCTION 

simPower.Moderated.Mediation = function (n=50, maxiter=50, slope1 =.20, slope2 = .20, e = .50, rho =.30) {

  # n       = number observation in treatment (is equal to number in control), sample size is 2n
  # maxiter = number of simulations
  # slope   = effect of mediator 1 and 2 on DV
  # e       = error: unexplained variance in DV (1 -R2)
  # rho     = auto correlation between the time points
  # flmer   = function LMER is used (TRUE) or else LME (FALSE)
  
  
  ## sample size per condition, so total sample is 2*n ##  
  nsub = 2*n
  
  ## indicator for intervention (effect coding)
  CGT <- c(rep(1,n),rep(0,n)) - .5
  
  ## autocorrelation (AR(1) matrix (rho is auto correlation)
  rho2 <- rho**2
  rho3 <- rho**3
  rho4 <- rho**4;
  
  cor=matrix(c(  1,  rho , rho2, rho3, rho4,
                 rho,    1 , rho , rho2, rho3,
                 rho2, rho ,    1, rho , rho2,
                 rho3, rho2, rho ,   1 , rho ,
                 rho4, rho3, rho2, rho ,   1), ncol=5, byrow=TRUE);
  
  respow <- res <- matrix(nrow=maxiter, ncol=8)
  colnames(res) <- c("a1", "a2", "b1", "b2","a1b1", "a2b2" ,"cp", "c")
  colnames(respow) <- c("Power(a1)","Power(a2)","Power(b1)","Power(b2)","Power(a1b1)","Power(a2b2)","Power(cp)","Power(c)")
  
  
# preparation for data construction (total variance is 1)

  e1 <- sqrt(e);   # e is percentage error in variance of DV
  b2 <- slope1;    # mediation effect mediator 1
  b3 <- slope2;    # mediation effect mediator 2
  b1 <- sqrt (1 - b2**2 - b3**2 - e);  # time effect
 

  ## START LOOP ##
      
for (i in 1:maxiter)
{
  
## data intervention condition ##

error <- mvrnorm(n=n,mu = c(0,0,0,0,0),Sigma=ID) 
iq <- rnorm(n,0,1)    ## moderator 
iq <- ((iq > 0)*1) - .5      ## cut off point, effect coding
no_effect <- mvrnorm(n=n,mu = c(0,0,0,0,0),Sigma=cor)
datam1 <- mvrnorm(n=n,mu = pattern_m,Sigma=cor)  
datam1[(iq == -.5),] <- no_effect[(iq == -.5),]   # add moderation effect for low IQ

dataw1 <- mvrnorm(n=n,mu = pattern_w,Sigma=cor)  
datay1 <- b1*mvrnorm(n=n,mu = pattern,Sigma=cor)  + b2* datam1 + b3*dataw1 + e1*error
colnames(datam1) <- predNamesm;
colnames(dataw1) <- predNamesw;
colnames(datay1) <- predNamesy;

## data control condition  ##

error <-  mvrnorm(n=n, mu = c(0,0,0,0,0),Sigma=ID)  
datam2 <- mvrnorm(n=n, mu = c(0,0,0,0,0),Sigma=cor) 
dataw2 <- mvrnorm(n=n, mu = c(0,0,0,0,0),Sigma=cor)  
datay2 <- b1*mvrnorm(n=n,mu = pattern0,Sigma=cor)  + b2* datam2 + + b3*dataw2 + e1*error
colnames(datam2) <- predNamesm;
colnames(dataw2) <- predNamesw;
colnames(datay2) <- predNamesy;


## missing pattern : T2 - T5  ##
mislev1 <- mislev*n
df <- nsub - 2*mislev1 - 2; # degrees of freedom 
df3 <- c((5*n - 1), (5*n-1))


  if (sum(mislev1) > 0) 
  {
datam1[1:(mislev1[1]),2] <- dataw1[1:(mislev1[1]),2] <- NA
datam2[1:(mislev1[1]),2] <- dataw2[1:(mislev1[1]),2] <- NA

datam1[1:(mislev1[2]),3] <- dataw1[1:(mislev1[2]),3] <- NA
datam2[1:(mislev1[2]),3] <- dataw2[1:(mislev1[2]),3] <- NA

datam1[1:(mislev1[3]),4] <- dataw1[1:(mislev1[3]),4] <- NA
datam2[1:(mislev1[3]),4] <- dataw2[1:(mislev1[3]),4] <- NA

datam1[1:(mislev1[4]),5] <- dataw1[1:(mislev1[4]),5] <- NA
datam2[1:(mislev1[4]),5] <- dataw2[1:(mislev1[4]),5] <- NA;
  }



## data under H1 ##
datam <- data.frame(rbind(datam1,datam2))
dataw <- data.frame(rbind(dataw1,dataw2))
datay <- data.frame(rbind(datay1,datay2))
datam[,"CGT"] <- CGT  
datam[,"IQ"]<- c(iq,iq)
datam[,"Subj"] <- c(1:nsub)
dlm <- melt(datam, id= c("Subj","CGT","IQ"), measure.vars=c("m1","m2","m3","m4","m5"), variable.name="time", value.name="m")
dlw <- melt(dataw,  measure.vars=c("w1","w2","w3","w4","w5"), variable.name="time", value.name="w")
dly <- melt(datay,  measure.vars=c("y1","y2","y3","y4","y5"), variable.name="time", value.name="y")
w <- dlw[,"w"]
y <- dly[,"y"]



## analyses ##

dl <- cbind(dlm,w,y)
levels(dl$time) <- c(-1/6,-1/6,-1/6,1/4,1/4)
result4<-lmer(y ~ 1 +  time*CGT  + (1 |Subj), data=dl,na.action = (na.omit)); # total effect: c 
res[i,8]<-fixef(result4)[4]     # c 
respow[i,8]<- (2*(1 - pt(coef(summary(result4))[4,"t value"], df=df[2])) < .05)*1


result1<-lmer(y ~  1 +  time*CGT + m + w + (1  |Subj), data=dl,na.action = (na.omit)); # mediation: b1+b2+cp path
 res[i,3:4]<- fixef(result1)[4:5]    # b1 and b2
 res[i,7] <- fixef(result1)[6]      # cp direct effect of intervention on T3 (interaction intervention x T3)
 respow[i,3:4]<- (2*(1 - pt(coef(summary(result1))[4:5,"t value"], df=df3)) < .05)*1
 respow[i,7]<- (2*(1 - pt(coef(summary(result1))[6,"t value"], df=df[4])) < .05)*1
 
dl <- cbind(dlm,w,y)
levels(dl$time) <- c(-.5,1/8,1/8,1/8,1/8)
result2<-lmer(m ~ 1 +   IQ*CGT*time  + (1 |Subj), data=dl,na.action = (na.omit)); # mediation: a1 path 
 res[i,1]<-fixef(result2)[8]     # a1 for high IQ
 respow[i,1]<- (2*(1 - pt(coef(summary(result2))[8,"t value"], df=df[1])) < .05)*1
 res[i,5] <- fixef(result2)[8]*fixef(result1)[4]     # a1*b1

dl <- cbind(dlm,w,y)
levels(dl$time) <- c(-1/4,-1/4,1/6,1/6,1/6)
 result3<-lmer(w ~ 1 +  time*CGT  + (1 |Subj), data=dl,na.action = (na.omit)); # mediation: a2 path 
 res[i,2]<-fixef(result3)[4]     # a2
 respow[i,2]<- (2*(1 - pt(coef(summary(result3))[4,"t value"], df=df[2])) < .05)*1
 res[i,6] <- fixef(result3)[4]*fixef(result1)[5]  ;   # a2*b2;



}
## end iteration ##

# compute significance levels for indirect effects


respow[,5:6] <- res[,5:6]/apply(res[,5:6],2,sd)
dfloop <- maxiter - 1
respow[,5] <- ((2*(1 - pt(respow[,5], df=dfloop))) < .05)*1
respow[,6] <- ((2*(1 - pt(respow[,6], df=dfloop))) < .05)*1

out <- apply(respow,2,mean)  # power estimates
out2 <- apply(res,2,mean)    # parameter estimates

output <- list(power=out, estimates=out2, raw=res, rawpower= respow)

}
# END FUNCTION



a <- simPower.Moderated.Mediation(n=50,maxiter=25,e=.50,slope1=.20,slope2=.20,rho=.3)


a[1]
a[2]
a[3]
