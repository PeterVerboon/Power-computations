
# require(lavaan);
# require('userfriendlyscience');
# require(MASS)
# require(dplyr)
# require(ggplot2)
# require(viridis)




#' Simulation for power of moderated mediation repeated measures (mmrm) model
#'
#' This function allows you to estimate the power for a given sample size 
#' of a repeated measures model with k measurements, one mediation effect and one interaction effect.
#' @param n sample size
#' @param EScond effect size condition
#' @param ESmod effect size moderator
#' @param ESint effect size interaction term
#' @param bpath vector of regression coefficients from mediator to latent intercept and slope respectively
#' @param ndepend number of dependent variables
#' @param rho vector with minimum and maximum value of randomly selected corrrelation between the variables
#' @param error standard deviations of the random error added to the data 
#' @param alpha alpha level
#' @param maxiter number of iterations
#' @keywords SEM latent growth mediation
#' @export 
#' @import MASS
#' @import lavaan
#' @return List with the following elements
#' @return power: estimate for all effects
#' @return bias: estimate for all effects
#' @return raw: the raw results of the simulation
#' @return input: the input paramters used in the simulation
#' @examples
#' simPwr.mmrm()
      simPwr.mmrm <- function(n=100, 
                               EScond = .2, 
                               ESmod = .2,
                               ESint = .2,
                               bpath = c(.4,.3,.2,.1),
                               ndepend = 4,
                               rho = c(0,0),
                               error = 1,
                               alpha = 0.05,
                               maxiter = 100) 
{   

  input <- c(EScond,ESmod,ESint,bpath,ndepend, rho, error,alpha,n, maxiter)
  blabel <- "b1"
  for (i in 2:length(bpath)){
    blabel <- paste0(blabel,",","b",i)
  }
  blabel <- (unlist(strsplit(blabel, ",", fixed = TRUE)))
  
  names(input) <- c("EScondition","ESmoderator","ESinteraction",blabel,"nDependent","rho","error","alpha","sampleSize","maxIter")
  
  ES <- c(EScond,ESmod,ESint,bpath,bpath[1]*EScond) 
  
# generate lavaan model with specifications
model0 <- buildSimModel(EScond = EScond,
                         ESmod = ESmod,
                         ESint = ESint,
                         bpath = bpath,
                         ndepend = ndepend,
                         model = 0)

# generate lavaan model used in analysis
model <- buildSimModel(EScond = "a1",
                        ESmod = "a2",
                        ESint = "a3",
                        bpath = blabel,
                        ndepend = ndepend,
                        model = 0)

res <- matrix(data=0,nrow=maxiter, ncol=20)

colnames(res) <- c("condition","power cond",
                   "moderation","power mod",
                   "interaction","power int",
                   "effect y1","power y1",
                   "effect y2", "power y2",
                   "effect y3","power y3",
                   "effect y4", "power y4",
                   "indirect", "power_ind",
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
  data3 <- data + data2
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
  res[i,11] <- parameterEstimates(result)[6,5]
  res[i,12] <- parameterEstimates(result)[6,8]
  res[i,13] <- parameterEstimates(result)[7,5]
  res[i,14] <- parameterEstimates(result)[7,8]
  res[i,15] <- filter(parameterEstimates(result), lhs %in% c("ind1"))[,5]
  res[i,16] <- filter(parameterEstimates(result), lhs %in% c("ind1"))[,8]
  res[i,c(17:20)] <- fitmeasures(result)[c("cfi","tli","rmsea", "srmr")]
  
  setTxtProgressBar(pb, i)
 }

  sigs <- res[,c(2,4,6,8,10,12,14,16)] < alpha
  power <- apply(sigs,2,mean)                                # power: count number of significant effects
  bias <- ES - apply(res[,c(1,3,5,7,9,11,13,15)],2,mean)                   # bias : mean of estimates
  
  output <- list(power = power, bias = bias, raw = res, input = input) 
  
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
                          maxiter = 100) 


### inspect results

showdist <- function(dat, var, plot = TRUE) {
  options(digits = 3)
  b <- dat[order(dat[,var]),]
  cat("mean  : ", mean(b[,var])," \n")
  cat("median: ", median(b[,var])," \n")
  cat("95% coverage interval: ", b[nrow(b)*c(0.05, 0.95),var]," \n")
  cat("99% coverage interval: ", b[nrow(b)*c(0.01, 0.99),var])
  if (plot) {plot(density(b[,var]), main=var, xlab= "")}
}

colnames(res$raw)

showdist(dat = res$raw, var = "condition", plot = TRUE)

res$power
res$bias
apply(res$raw,2,sd)


##
## loop over sample sizes and effect sizes

samSizes <- seq(from=100, to=600, by=50)
esSizes <- c(.15,.30,.50)
out <- matrix(0,nrow=(length(samSizes)*length(esSizes)),ncol=10)
colnames(out) <- c("ES","N","condition","moderation","interaction",
                   "effect_y1","effect_y2","effect_y3","effect_y3","effect_ind")
out[,1] <- rep(esSizes,each = length(samSizes))
out[,2] <- rep(samSizes,length(esSizes))

for (es in esSizes) {
  for (n in samSizes) {
    
    cat("\n","ES =",es, " ##  N =",n, "\n")
    
    res <- simPower.Mediation(n=n, 
                              EScond = es, 
                              ESmod = .3,
                              ESint = .2,
                              bpath = c(.4,.3,.2,.1), 
                              rho = c(0.0,0.0),
                              error = 1,
                              alpha = 0.05,
                              maxiter = 500) 
    
    out[out[,"N"] == n & out[,"ES"] == es,-c(1,2)] <- res$power
  }
}

save(out, file="modmed_alpah05.Rdata")

## plot the results

out1 <- data.frame(out)
out1$EffectSize <- as.factor(out1$ES)

p <- ggplot(data=out1, aes(y=effect_y1, x=N,  colour=EffectSize)) +
  geom_point(size=2) + geom_line(size=1) +
  geom_hline(yintercept=0.80, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.90, linetype="dashed", color = "blue") +
  scale_y_continuous(breaks=seq(0.10, 1, 0.10)) + scale_x_continuous(breaks=seq(100, 600, 50) ) +
  scale_color_viridis(discrete=TRUE) +
  theme_bw(base_size = 14) +
  ggtitle("Power of moderated mediation lavaan model, alpha=0.05") +
  theme(plot.title = element_text(size=10, hjust=0)) 
p

# save the plot as pdf
ggsave(plot=p,filename="Modmed_alpha05_b1.pdf", width=7, height=5)




require(semPlot)
result <- sem(model, data)
semPaths(result, layout = "tree2")


