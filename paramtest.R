
# load all packages used in this vignette
library('paramtest')
library('pwr')
library('ggplot2')
library('knitr')
library('nlme')
library('lavaan')
library('dplyr')
library(MASS)
library(lm.beta)

options(digits=2, scipen=999)


lm_test_interaction <- function(simNum=100, N=100, b1=.5, b2=.3, b3=.2, b0=0, x1m=0, x1sd=1,
                                x2m=0, x2sd=1, rho=0) {
  
  require(MASS)
  require(lm.beta)
  

  sigma = matrix(c(x1sd,rho,rho,x2sd), nrow=2, ncol=2)
  dat <- as.data.frame(mvrnorm(n = N, Sigma=sigma, mu=c(x1m,x2m)))
  names(dat) <- c("x1", "x2")
  
  modelVar <- b1^2 + b2^2 + b3^2
  if (modelVar < 1) { yvar <- sqrt(1 - modelVar)    
     } else  { yvar <- 0; print("warning: model variance equals or exceeds 1") }
  
  dat$y <- rnorm(N, b0 + b1*dat$x1 + b2*dat$x2 + b3*dat$x1*dat$x2, yvar)
  model <- lm(y ~ x1 * x2, data=dat)
  
  # pull output from model (two main effects and interaction)
  est_x1 <- coef(summary(model))['x1', 'Estimate']
  p_x1 <- coef(summary(model))['x1', 'Pr(>|t|)']
  sig_x1 <- p_x1 < .05
  est_x2 <- coef(summary(model))['x2', 'Estimate']
  p_x2 <- coef(summary(model))['x2', 'Pr(>|t|)']
  sig_x2 <- p_x2 < .05
  est_int <- coef(summary(model))['x1:x2', 'Estimate']
  p_int <- coef(summary(model))['x1:x2', 'Pr(>|t|)']
  sig_int <- p_int < .05
  rsq <- summary(model)$r.squared
  rsq.adj <- summary(model)$adj.r.squared
  
  res <- c(est_x1, p_x1, sig_x1, 
           est_x2, p_x2, sig_x2,
           est_int, p_int, sig_int,
           rsq, rsq.adj, sd(dat$y),
           lm.beta(model)$standardized.coefficients[2:3])
  names(res) <- c("est_x1", "p_x1", "sig_x1", 
                  "est_x2", "p_x2", "sig_x2",
                  "est_int", "p_int", "sig_int",
                  "rsq", "rsq.adj","sdy",
                  "std_x1", "std_x2")
  return(res)
}

a <- lm_test_interaction(simNum=100, N=100, b1=.221, b2=.171, b3=.146, b0=0, x1m=0, x1sd=1,
                    x2m=0, x2sd=1, rho=.3)

# varying N at 200 and 300; setting coefficient of x1 = .15, coefficient of
# x2 = 0, and coefficien of interaction = .3

power_lm_int <- grid_search(lm_test_interaction, params=list(N=c(100,200, 300, 400)),
                            n.iter=100, output='data.frame', b1=.221, b2=.171, b3=.146, rho = .3, parallel='snow', ncpus=4)
results(power_lm_int) %>%
  group_by(N.test) %>%
  summarise(
    power_x1=mean(sig_x1),
    power_x2=mean(sig_x2),
    power_int=mean(sig_int),
    rsq = mean(rsq),
    rsq.adj = mean(rsq.adj))


print(apply(results(power_lm_int),2,mean), digits=1)


# test with grid_search

power_lm_int <- grid_search(regrPwrSim, params=list(n=c(100,200)),
                            predictors= 3,
                            n.iter= 10,
                            betas= c(.1,.2,.3), 
                            cor= c(.1,.2,.3), 
                            interactions = NULL,
                            output= 'data.frame', 
                            parallel= 'snow', ncpus= 4)

results(power_lm_int) %>%
  group_by(n.test) %>%
  summarise(
    power_x1=mean(sig_x1),
    power_x2=mean(sig_x2),
    power_int=mean(sig_int),
    rsq = mean(rsq),
    rsq.adj = mean(rsq.adj))


print(apply(results(power_lm_int),2,mean), digits=1)
