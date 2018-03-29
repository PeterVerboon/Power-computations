


## Power computation for intervention-control design
## 5 measurements (pre (T0), post (T1), Follow-up1 (T2), Follow_up2 (T3), Follow_up3 (T4) )
## 1 moderator
## 1 mediator
## Effect sizes and sample size can be adjusted
## Power is computed at requested level
##
##
##  Author:   P. Verboon
##  Date:     2017, July
##  Adapted:  2018, March
##############################################################################

# Difference scores are computed from T0 for the three following waves

require(lavaan);
require('userfriendlyscience');
require(MASS)

options(digits=3);


## Function definition for power simulation

 
simPower.Mediation <- function(n=100, ESt1=.2, ESmod=.2, bpath1 = .4,bpath2 = .3,bpath3 = .2, bpath4 = .1, maxiter=100) {   
  
  res <- matrix(data=0,nrow=maxiter, ncol=12)
  
  colnames(res) <- c("cond-M","p-value cond","moderation effect","p-value mod",
                     "effect T1","p-value 1","effect T2", "p-value 2",
                     "effect T3", "p-value 3","effect T4", "p-value 4")
  
  for (i in 1:maxiter) {  
    
  print(i)
  data <- simDataSet(n, varNames=c('mod1', 'Med', 'y1', 'y2', 'y3','y4'),
                     means = c(0,0, 0, 0, 0,0), sds = c(1,1, 1, 1, 1,1), silent=TRUE, seed = NULL,
                     correlations = c(.30, .60),
                     specifiedCorrelations =
                       list(c('mod1', 'y1', .1),
                            c('mod1', 'y2', .1),
                            c('mod1', 'y3', .1),
                            c('mod1', 'y4', .1),
                            c('mod1', 'Med',  .2),
                            c('Med', 'y1', bpath1),
                            c('Med', 'y2', bpath2),
                            c('Med', 'y3', bpath3),
                            c('Med', 'y4', bpath4))
                      );


  data$cond <- c(rep(-.5,(n/2)),rep(.5,(n/2)))
  data$int1 <- data$cond * data$mod1

  data$Med[1:(n/2)] <- data$Med[1:(n/2)] + ESt1   # add effect to first group in M
  data$Med <- data$Med + ESmod*data$int1          # construct interaction effect to mediator

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
  res[i,11] <- parameterEstimates(result)[7,5]
  res[i,12] <- parameterEstimates(result)[7,8]
  
  
  }
  
  return(res)
  
  }  # END FUNCTION
  

################################ LAVAAN MODEL  ###########

model<-"
! regressions 
Med ~ a1*cond
Med ~ a2*mod1
Med ~ a3*int1

y1~b1*Med
y2~b2*Med
y3~b3*Med
y4~b4*Med

! residuals, variances and covariances

y1 ~~  e1* y2
y1 ~~  e2* y3
y1 ~~  e2* y4
y2 ~~  e2* y3
y2 ~~  e2* y4
y3 ~~  e2* y4

";


######################################################



#### Do the power computations


es <- c(-0.5, 0.5, 0.4, 0.3, 0.2, 0.1)         # vector with the effect sizes

res <- simPower.Mediation(n=150, maxiter=500, ESt1 =.5, ESmod=.5, 
                          bpath1 = .4, bpath2 = .3, bpath3 = .2, bpath4 = .2)



b <- res[,c(2,4,6,8,10,12)] < .05



es - apply(res[,c(1,3,5,7,9,11)],2,mean)       # bias : mean of estimates

apply(b,2,mean)                                # power: count number of significant effects












