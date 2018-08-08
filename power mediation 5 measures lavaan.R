

##############################################################################
## Power computation for intervention-control design
## 5 measurements (y1 to y4 )
## 1 moderator for the condition --> mediator effect
## 1 mediator
## Effect sizes and sample size can be adjusted
## 
## Result: power and bias for each parameter in the model
##
##
##  Author:   P. Verboon
##  Date:     2018, July
##############################################################################


require(lavaan);
require('userfriendlyscience');
require(MASS)

options(digits=3);



#######   Define Function BuildModel   ###############################################

buildSimModel <- function (EScond = .2, 
                           ESmod = .2,
                           ESint = .2,
                           bpath = c(.4,.3,.2,.1)) 
  {
  
  
  modela1 <- paste0("mediator", " ~ " ,EScond,"*","condition" ,  collapse = " \n ") 
  modela2 <- paste0("mediator", " ~ " ,ESmod,"*","moderator" ,  collapse = " \n ") 
  modela3 <- paste0("mediator", " ~ " ,ESint,"*","interaction" ,  collapse = " \n ") 
  
  modelb1 <- paste0("y1", " ~ " ,bpath[1],"*","mediator" ,  collapse = " \n ") 
  modelb2 <- paste0("y2", " ~ " ,bpath[2],"*","mediator" ,  collapse = " \n ") 
  modelb3 <- paste0("y3", " ~ " ,bpath[3],"*","mediator" ,  collapse = " \n ") 
  modelb4 <- paste0("y4", " ~ " ,bpath[4],"*","mediator" ,  collapse = " \n ") 
  
  model <- paste0(modela1," \n ",modela2," \n ",modela3," \n ",
                  modelb1," \n ", modelb2," \n ",modelb3," \n ", modelb4," \n ")
                
return(model)
  
}  # end function



########### Build LAVAAN model to analyze  ###########

model<-"

mediator ~ a1*condition
mediator ~ a2*moderator
mediator ~ a3*interaction

y1~b1*mediator
y2~b2*mediator
y3~b3*mediator
y4~b4*mediator

";


### Define simulation function ###

simPower.Mediation <- function(n=100, 
                               EScond = .2, 
                               ESmod = .2,
                               ESint = .2,
                               bpath = c(.4,.3,.2,.1), 
                               rho = c(0,0),
                               error = 1,
                               alpha = 0.05,
                               maxiter = 100) 
{   

 ES <- c(EScond,ESmod,ESint,bpath) 
  
# generate lavaan model with specifications
model0 <- buildSimModel(EScond = EScond,
                        ESmod = ESmod,
                        ESint = ESint,
                        bpath = bpath)

res <- matrix(data=0,nrow=maxiter, ncol=14)

colnames(res) <- c("condition","p-value cond",
                   "moderation","p-value mod",
                   "interaction","p-value int",
                   "effect y1","p-value y1",
                   "effect y2", "p-value y2",
                   "effect y3","p-value y3",
                   "effect y4", "p-value y4")

for (i in 1:maxiter) {  
  
print(i)
  
# simulate data according to model specifications
data <- simulateData(model0, sample.nobs = n)

# simulate random data 
  data2 <- simDataSet(n=n, varNames=c('mediator', 'y1', 'y2', 'y3','y4',"condition",'moderator', "interaction"),
                    means = 0, sds = error, 
                    ranges = list(c(condition = c(-3.0, 3.0))),
                    silent=TRUE, seed = NULL, 
                    correlations = c(rho[1], rho[2]))

# add random error
  data3 <- data + data2
  data3$condition <- cut(data3$condition, breaks = 2)

# analyse data using lavaan
  result <- sem(model, data3)

  res[i,1] <- parameterEstimates(result)[1,5]
  res[i,2] <- parameterEstimates(result)[1,8]
  res[i,3] <- parameterEstimates(result)[2,5]
  res[i,4] <- parameterEstimates(result)[2,8]
  res[i,5] <- parameterEstimates(result)[3,5]
  res[i,6] <- parameterEstimates(result)[3,8]
  res[i,7] <- parameterEstimates(result)[4,5]
  res[i,8] <- parameterEstimates(result)[4,8]
  res[i,9] <- parameterEstimates(result)[5,5]
  res[i,10] <- parameterEstimates(result)[5,8]
  res[i,11] <- parameterEstimates(result)[6,5]
  res[i,12] <- parameterEstimates(result)[6,8]
  res[i,13] <- parameterEstimates(result)[7,5]
  res[i,14] <- parameterEstimates(result)[7,8]
  
 }

  sigs <- res[,c(2,4,6,8,10,12,14)] < alpha
  power <- apply(sigs,2,mean)                                # power: count number of significant effects
  bias <- ES - apply(res[,c(1,3,5,7,9,11,13)],2,mean)                   # bias : mean of estimates
  
  output <- list(power = power, bias = bias, raw = res) 
  
  return(output)
  
  
  }  # END FUNCTION
  




#### Do the power computations



res <- simPower.Mediation(n=300, 
                          EScond = .4, 
                          ESmod = .3,
                          ESint = .2,
                          bpath = c(.4,.3,.2,.1), 
                          rho = c(0.1,0.1),
                          error = 1,
                          alpha = 0.05,
                          maxiter = 200) 

res$power
res$bias
apply(res$raw,2,sd)










