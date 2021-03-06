#' Simulation for power of two lavaan models
#'
#' This function allows you to estimate the power for a given sample size 
#' of a latent growth model with k measurements, one mediation effect and one interaction effect.
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
#' simPwr.growth()
simPwr.growth <- function(n=200, 
                         EScond = .2, 
                         ESmod = .2,
                         ESint = .2,
                         bpath = c(0.4,0.1), 
                         ndepend = 5,
                         rho = c(0.1,0.3),
                         error = 1,
                         alpha = 0.05,
                         maxiter = 1000) 
{   
  
  input <- c(EScond,ESmod,ESint,bpath,error,alpha, n, maxiter)
  blabel <- "b1"
  for (i in 2:length(bpath)){
    blabel <- paste0(blabel,",","b",i)
  }
  blabel <- (unlist(strsplit(blabel, ",", fixed = TRUE)))
  
  names(input) <- c("EScondition","ESmoderator","ESinteraction",blabel,"error","alpha","sampleSize", "maxIter")
  
  ES <- c(EScond,ESmod,ESint,bpath,bpath*EScond) 
  
  
  # generate lavaan model with specifications
  model0 <- buildSimModel(EScond = EScond,
                          ESmod = ESmod,
                          ESint = ESint,
                          bpath = bpath,
                          ndepend = ndepend,
                          model = "growth")
  
  # generate lavaan model used in analysis
  model <- buildSimModel(EScond = "a1",
                         ESmod = "a2",
                         ESint = "a3",
                         bpath = c("b1","b2"),
                         ndepend = ndepend,
                         model = "growth")
  
  res <- matrix(data=0,nrow=maxiter, ncol=18)
  
  colnames(res) <- c("condition","sig cond",
                     "moderation","sig mod",
                     "interaction","sig int",
                     "effect li","sig li",
                     "effect ls", "sig ls",
                     "indirect i","sig ind_li",
                     "indirect s","sig ind_ls",
                     "cfi","tli","rmsea","srmr")
  
  # Initiate the Progress bar
  pb <- txtProgressBar(min = 0, max = maxiter, style = 3)
  set.seed(1234)
  
  for (i in 1:maxiter) {  
    
    # simulate data according to model specifications
    data <- lavaan::simulateData(model0, sample.nobs = n, model.type = "growth")
    k <- dim(data)[2]
    
    # simulate random data 
    sigma <- matrix(runif(k*k, min = rho[1], max = rho[2]), nrow = k)
    diag(sigma) <- error
    data2 <-  MASS::mvrnorm(n = n, mu = rep(0,k), Sigma = sigma)
    
    # add random error
    data3 <- data  + data2
    data3$condition <- as.numeric(cut(data3$condition, breaks = 2))
    
    # analyse data using lavaan
    result <- lavaan::sem(model, data3)
    
    res[i,1] <- lavaan::parameterEstimates(result)[1,5]
    res[i,2] <- lavaan::parameterEstimates(result)[1,8]
    res[i,3] <- lavaan::parameterEstimates(result)[2,5]
    res[i,4] <- lavaan::parameterEstimates(result)[2,8]
    res[i,5] <- lavaan::parameterEstimates(result)[3,5]
    res[i,6] <- lavaan::parameterEstimates(result)[3,8]
    res[i,7] <- lavaan::parameterEstimates(result)[4,5]
    res[i,8] <- lavaan::parameterEstimates(result)[4,8]
    res[i,9] <- lavaan::parameterEstimates(result)[5,5]
    res[i,10] <- lavaan::parameterEstimates(result)[5,8]
    res[i,11] <- lavaan::parameterEstimates(result)[(lavaan::parameterEstimates(result)[,"lhs"] == "ind1"),5]
    res[i,12] <- lavaan::parameterEstimates(result)[(lavaan::parameterEstimates(result)[,"lhs"] == "ind1"),8]
    res[i,13] <- lavaan::parameterEstimates(result)[(lavaan::parameterEstimates(result)[,"lhs"] == "ind2"),5]
    res[i,14] <- lavaan::parameterEstimates(result)[(lavaan::parameterEstimates(result)[,"lhs"] == "ind2"),8]
    res[i,c(15:18)] <- lavaan::fitmeasures(result)[c("cfi","tli","rmsea", "srmr")]
    
    setTxtProgressBar(pb, i)
  }
  
  sigs <- res[,c(2,4,6,8,10,12,14)] < alpha
  power <- apply(sigs,2,mean)                                           # power: count number of significant effects
  bias <- ES - apply(res[,c(1,3,5,7,9,11,13)],2,mean)                   # bias : mean of estimates
  
  output <- list(power = power, bias = bias, raw = res, input = input) 
  
  class(output) <- "simPwr.growth"
  return(output)
  
} 





#' Print method for simPwr.growth
#'
#' This function allows you to print and plot the result of simPwr.growth.
#' @param x simPwr.growth object
#' @param var the effect that is printed and optionally plotted. Use dimnames(x$raw)[[2]] to see which names are available.
#' @param plot whether a plot is shown (plot = TRUE)
#' @export 
#' @keywords SEM latent growth mediation
#' @examples print(x, var= "indirect li")
print.simPwr.growth <- function(x, var, plot = TRUE) {
  dat <- x$raw
  b <- dat[order(dat[,var]),]
  cat("mean  : ", mean(b[,var])," \n")
  cat("median: ", median(b[,var])," \n")
  cat("95% coverage interval: ", b[nrow(b)*c(0.05, 0.95),var]," \n")
  cat("99% coverage interval: ", b[nrow(b)*c(0.01, 0.99),var])
  if (plot) {plot(density(b[,var]), main=var, xlab= "")}
}


## simPwr.growth(n=100, rho = c(.1,.4), bpath = c(.4,.3, .1))
