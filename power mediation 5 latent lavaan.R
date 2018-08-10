

##############################################################################
## Power computation for intervention-control design
## 5 measurements (y1 to y4 )
## 1 moderator for the condition --> mediator effect
## 1 mediator
## latent growth model
## Effect sizes and sample size can be adjusted
## 
## Result: power and bias for each parameter in the model
##
##
##  Author:   P. Verboon
##  Date:     2018, August
##############################################################################


require(lavaan);
require('userfriendlyscience');
require(MASS)
require(dplyr)
require(ggplot2)
require(viridis)

options(digits=3);



#######   Define Function BuildModel   ###############################################

buildSimModel <- function (EScond = .2, 
                           ESmod = .2,
                           ESint = .2,
                           bpath = c(.4,.3),
                           model = 0) 
  {
  
  modela1 <- paste0("mediator", " ~ " ,EScond,"*","condition" ,  collapse = " \n ") 
  modela2 <- paste0("mediator", " ~ " ,ESmod,"*","moderator" ,  collapse = " \n ") 
  modela3 <- paste0("mediator", " ~ " ,ESint,"*","interaction" ,  collapse = " \n ") 
  
  modelb1 <- paste0("li", " ~ " ,bpath[1],"*","mediator" ,  collapse = " \n ") 
  modelb2 <- paste0("ls", " ~ " ,bpath[2],"*","mediator" ,  collapse = " \n ") 

  modelc1 <- ("li =~ 1*y1 + 1*y2 + 1*y3 + 1*y4")
  modelc2 <- ("ls =~ 3*y1 + 2*y2 + 1*y3 + 0*y4")
  
  modelind1 <- ifelse(model == 0, " ", paste0("ind1 := ",EScond,"*",bpath[1],collapse = " \n "))
  modelind2 <- ifelse(model == 0, " ", paste0("ind2 := ",EScond,"*",bpath[2],collapse = " \n "))
  
  model <- paste0(modela1," \n ",modela2," \n ",modela3," \n ",
                  modelb1," \n ", modelb2," \n ",
                  modelc1," \n ", modelc2," \n ",
                  modelind1," \n ", modelind2 )
return(model)
  
}  # end function




### Define simulation function ###

simPower.Mediation <- function(n=200, 
                               EScond = .2, 
                               ESmod = .2,
                               ESint = .2,
                               bpath = c(.4,.3), 
                               rho = c(0,0),
                               error = 1,
                               alpha = 0.05,
                               maxiter = 100) 
{   

 ES <- c(EScond,ESmod,ESint,bpath,bpath*EScond) 
  
# generate lavaan model with specifications
model0 <- buildSimModel(EScond = EScond,
                        ESmod = ESmod,
                        ESint = ESint,
                        bpath = bpath,
                        model = 0)

# generate lavaan model used in analysis
model <- buildSimModel(EScond = "a1",
                        ESmod = "a2",
                        ESint = "a3",
                        bpath = c("b1","b2"),
                        model = 1)

res <- matrix(data=0,nrow=maxiter, ncol=18)

colnames(res) <- c("condition","pvalue cond",
                   "moderation","pvalue mod",
                   "interaction","pvalue int",
                   "effect li","pvalue li",
                   "effect ls", "pvalue ls",
                   "indirect i","pvalue ind_li",
                   "indirect s","pvalue ind_ls",
                   "cfi","tli","rmsea","srmr")
  
# Initiate the Progress bar
pb <- txtProgressBar(min = 0, max = maxiter, style = 3)

for (i in 1:maxiter) {  

  
# simulate data according to model specifications
data <- simulateData(model0, sample.nobs = n)

# simulate random data 
   data2 <- simDataSet(n=n, varNames=c('mediator', 'y1', 'y2', 'y3','y4',"condition",'moderator', "interaction"),
                     means = 0, sds = error, 
                     ranges = list(c(condition = c(-3.0, 3.0))),
                     silent=TRUE, seed = NULL, 
                     correlations = c(rho[1], rho[2]))

# add random error
  data3 <- data  + data2
  data3$condition <- as.numeric(cut(data3$condition, breaks = 2))

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
  res[i,11] <- filter(parameterEstimates(result), lhs %in% c("ind1"))[,5]
  res[i,12] <- filter(parameterEstimates(result), lhs %in% c("ind1"))[,8]
  res[i,13] <- filter(parameterEstimates(result), lhs %in% c("ind2"))[,5]
  res[i,14] <- filter(parameterEstimates(result), lhs %in% c("ind2"))[,8]
  res[i,c(15:18)] <- fitmeasures(result)[c("cfi","tli","rmsea", "srmr")]
  
 
  setTxtProgressBar(pb, i)
 }

  sigs <- res[,c(2,4,6,8,10,12,14)] < alpha
  power <- apply(sigs,2,mean)                                # power: count number of significant effects
  bias <- ES - apply(res[,c(1,3,5,7,9,11,13)],2,mean)                   # bias : mean of estimates
  
  output <- list(power = power, bias = bias, raw = res) 
  
  return(output)
  
  
  }  # END FUNCTION
  



##################################################

## Do the power computations for fixed sample size and fixed effect size

res <- simPower.Mediation(n=200, 
                          EScond = .5, 
                          ESmod = .3,
                          ESint = .2,
                          bpath = c(.4,.3), 
                          rho = c(0.1,0.1),
                          error = 1,
                          alpha = 0.05,
                          maxiter = 1000) 


### inspect results

showdist <- function(dat, var, plot = TRUE) {
   options(digits = 3)
   b <- res[order(res[,var]),]
   cat("mean  : ", mean(b[,var])," \n")
   cat("median: ", median(b[,var])," \n")
   cat("95% coverage interval: ", b[maxiter*c(0.05, 0.95),var]," \n")
   cat("99% coverage interval: ", b[maxiter*c(0.01, 0.99),var])
   if (plot) {plot(density(b[,var]), main=var, xlab= "")}
}

colnames(res$raw)

showdist(dat = res$raw, var = "indirect s", plot = TRUE)

res$power
res$bias
apply(res$raw,2,sd)
apply(res$raw,2,mean)




##
## loop over sample sizes and effect sizes

 samSizes <- seq(from=100, to=600, by=50)
 esSizes <- c(.15,.30,.50)
 out <- matrix(0,nrow=(length(samSizes)*length(esSizes)),ncol=9)
 colnames(out) <- c("ES","N","condition","moderation", "interaction","effect_li", "effect_ls", "indirect_i","indirect_s" )
 out[,1] <- rep(esSizes,each = length(samSizes))
 out[,2] <- rep(samSizes,length(esSizes))
 
 for (es in esSizes) {
     for (n in samSizes) {
       
  cat("\n","ES =",es, " ##  N =",n, "\n")
       
  res <- simPower.Mediation(n=n, 
                            EScond = es, 
                            ESmod = .3,
                            ESint = .2,
                            bpath = c(.5,.2), 
                            rho = c(0.0,0.0),
                            error = 1,
                            alpha = 0.05,
                            maxiter = 500) 
  
  out[out[,"N"] == n & out[,"ES"] == es,-c(1,2)] <- res$power
     }
}

 save(out, file="modmedlat_alpah05.Rdata")
 
 ## plot the results
 
 out1 <- data.frame(out)
 out1$EffectSize <- as.factor(out1$ES)
 
 p <- ggplot(data=out1, aes(y=condition, x=N,  colour=EffectSize)) +
   geom_point(size=2) + geom_line(size=1) +
   geom_hline(yintercept=0.80, linetype="dashed", color = "red") +
   geom_hline(yintercept=0.90, linetype="dashed", color = "blue") +
   scale_y_continuous(breaks=seq(0.10, 1, 0.10)) + scale_x_continuous(breaks=seq(100, 600, 50) ) +
   scale_color_viridis(discrete=TRUE) +
   theme_bw(base_size = 14) +
   ggtitle("Power of moderated mediation latent growth model, alpha=0.05") +
   theme(plot.title = element_text(size=10, hjust=0)) 
 p

 # save the plot as pdf
 ggsave(plot=p,filename="Modmedlatent_alpha05_cond.pdf", width=7, height=5)


