


# intervention-control design
# 4 measurements (pre (T0), post (T1), Follow-up1 (T2), Follow_up2 (T3) )
# 2 moderators
# effect at T3 is small-medium
# effect moderators is small

# Difference scores are computes from T0 for the three following waves

require(lavaan);
require('userfriendlyscience');
require(MASS)

options(digits=3);

rho <-  .20
rho2 <- .00
rho3 <- .50
rho4 <- .30

corMat = matrix(c(  1,  rho , rho2, rho3, rho4,
                  rho,    1 , rho , rho, rho,
                 rho2, rho ,    1,  rho , rho,
                 rho3, rho, rho ,    1 , rho ,
                 rho4, rho, rho,  rho ,   1), ncol=5, byrow=TRUE);



 
simPower.Mediation <- function(n=100, ESt1=.2, ESmod=.2, bpath1 = .4,bpath2 = .3,bpath3 = .2, maxiter=100) {   
  
  res <- matrix(data=0,nrow=maxiter, ncol=10)
  
  colnames(res) <- c("cond-M","p-value cond","moderation effect","p-value mod","effect T1","p-value 1","effect T2", "p-value 2","effect T3", "p-value 3")
  
  for (i in 1:maxiter) {  
    
  data <- simDataSet(n, varNames=c('mod1', 'M', 'y1', 'y2', 'y3'),
                     means = c(0,0,  0, 0, 0), sds = c(1,1, 1, 1, 1), silent=TRUE, seed = NULL,
                     correlations = c(.30, .60),
                     specifiedCorrelations =
                       list(c('mod1', 'y1', .1),
                            c('mod1', 'y2', .1),
                            c('mod1', 'y3', .1),
                            c('mod1', 'M', .2),
                            c('M', 'y1', bpath1),
                            c('M', 'y2', bpath2),
                            c('M', 'y3', bpath3))
                      );

  
  data$cond <- c(rep(-.5,(n/2)),rep(.5,(n/2)))
  data$int1 <- data$cond * data$mod1

  data$M[1:(n/2)] <- data$M[1:(n/2)] + ESt1  # add effect to first group in M
  data$M <- data$M + ESmod*data$int1         # construct interaction effect to mediator

  result<-sem(model, data=data, fixed.x=F, missing="FIML");
  
  res[i,1] <- parameterEstimates(result)[1,5]
  res[i,2] <- parameterEstimates(result)[1,8]
  res[i,3] <- parameterEstimates(result)[3,5]
  res[i,4] <- parameterEstimates(result)[3,8]
  res[i,5] <- parameterEstimates(result)[4,5]
  res[i,6] <- parameterEstimates(result)[4,8]
  res[i,7] <- parameterEstimates(result)[5,5]
  res[i,8] <- parameterEstimates(result)[5,8]
  res[i,9] <- parameterEstimates(result)[6,5]
  res[i,10] <- parameterEstimates(result)[6,8]
  
  
  }
  
  return(res)
  
  }  # END FUNCTION
  


res <- simPower.Mediation(n=250, maxiter=500, ESt1 =.5, bpath1 = .4, bpath2 = .3, bpath3 = .2, ESmod=.5)


b <- res[,c(2,4,6,8,10)] < .05

apply(res[,c(1,3,5,7,9)],2,mean)  # bias : mean of estimates
apply(b,2,mean)    # power: count number of significant effects




summary(data)
associationMatrix(data)
var(data)

# ES small = .2,  medium = .5



# Independent groups t-test
meanDiff(data$y1 ~ data$cond, paired = FALSE, r.prepost = NULL, var.equal = "test", plot = TRUE, digits = 2, envir = parent.frame())





################################ LAVAAN MODEL  ###########

model<-"
! regressions 
M ~ a1*cond
M ~ a2*mod1
M ~ a3*int1

y1~b1*M
y2~b2*M
y3~b3*M


! residuals, variances and covariances

y1 ~~  e1* y2
y1 ~~  e2* y3
y2 ~~  e2* y3

! means
# y1~0
# y2~0

";


###################################

result<-sem(model, data=data, fixed.x=F, missing="FIML");
summary(result, fit.measures=TRUE,standardized=TRUE)




res[i,1] <- parameterEstimates(result)[1,5]
res[i,2] <- parameterEstimates(result)[1,8]
res[i,3] <- parameterEstimates(result)[3,5]
res[i,4] <- parameterEstimates(result)[3,8]
res[i,5] <- parameterEstimates(result)[4,5]
res[i,6] <- parameterEstimates(result)[4,8]
res[i,7] <- parameterEstimates(result)[6,5]
res[i,8] <- parameterEstimates(result)[6,8]

#
parameterEstimates(result)
parameterEstimates(result)[c(1,3,4,5,6),5]
