
# setwd("~/Documents/Open Universiteit/Onderzoek/Project Ank ")
# getwd()



# intervention-control design
# 3 measurements (pre, post1, post2)
# 2 moderators
# effect at T3 is small-medium
# effect moderators is small

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



 
simPower <- function(n=100, ESt1=.2, ESt2=.2, ESmod=.2, maxiter=100) {   
  
  res <- matrix(data=0,nrow=maxiter, ncol=8)
  
  colnames(res) <- c("main effect T2","p-value","moderation effect T2","p-value","main effect T1","p-value","moderation effect T1", "p-value")
  
  for (i in 1:maxiter) {  
    
  data <- simDataSet(n, varNames=c('mod1', 'y0', 'y1', 'y2'),
                     correlations = c(.10, .25),means = c(0,  0, 0, 0), sds = c(1, 1, 1, 1), silent=TRUE, seed = NULL);

  data$cond <- c(rep(1,(n/2)),rep(0,(n/2)))
  data$int1 <- 2*(data$cond - .5) * data$mod1

  data$y2[1:(n/2)] <- data$y2[1:(n/2)] + ESt2
  data$y1[1:(n/2)] <- data$y1[1:(n/2)] + ESt1
  data$y2 <- data$y2 + ESmod*data$int1
  data$y1 <- data$y1 + ESmod*data$int1
  
  result<-sem(model, data=data, fixed.x=F, missing="FIML");
  
  res[i,1] <- parameterEstimates(result)[1,5]
  res[i,2] <- parameterEstimates(result)[1,8]
  res[i,3] <- parameterEstimates(result)[3,5]
  res[i,4] <- parameterEstimates(result)[3,8]
  res[i,5] <- parameterEstimates(result)[4,5]
  res[i,6] <- parameterEstimates(result)[4,8]
  res[i,7] <- parameterEstimates(result)[6,5]
  res[i,8] <- parameterEstimates(result)[6,8]
  
  
  }
  
  return(res)
  
  }  # END FUNCTION
  


res <- simPower(n=200, maxiter=200, ESt1 =.2, ESt2=.2, ESmod=.2)


b <- res[,c(2,4,6,8)] < .05

apply(res,2,mean)  # bias : mean of estimates
apply(b,2,mean)    # power: count number of significant effects




summary(data)
associationMatrix(data)
var(data)

# ES small = .2,  medium = .5



# Independent groups t-test
meanDiff(data$y2 ~ data$cond, paired = FALSE, r.prepost = NULL, var.equal = "test", plot = TRUE, digits = 2, envir = parent.frame())





################################ LAVAAN MODEL  ###########

model<-"
! regressions 
y2~a1*cond
y2~a2*mod1
y2~a3*int1

y1~b1*cond
y1~b2*mod1
y1~b3*int1
y1~b4*y0

! residuals, variances and covariances

y1 ~~  e1* y2
#y0 ~~  e2* y2

! means
 y1~0
 y2~0

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
parameterEstimates(result)[4,5]
