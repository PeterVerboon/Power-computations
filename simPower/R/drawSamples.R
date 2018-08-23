drawSamples <-
function(preds = NULL,
                        Sigma = NULL,
                        coef = NULL,
                        predictNames = NULL,
                        dependName = "y",
                        repetitions = 100, samples=100,
                        interactionTerms = NULL) {
  
  
  output <- matrix(data=0, nrow=samples, ncol=(preds + length(interactionTerms) + 1))
  mu0 = rep(0, preds)
  
  for (i in 1:samples) {
    dat <- data.frame( mvrnorm(n = repetitions, mu = mu0, Sigma = Sigma))
    
    colnames(dat) <- c(predictNames)
    
    regressionFormula <- paste(dependName, "~",
                               paste0(predictNames, collapse=" + "));
    
    
    if (!is.null(interactionTerms)) {
      for (interaction in 1:length(interactionTerms)) {
        
        names(interactionTerms)[interaction] <- paste0(interactionTerms[[interaction]][1],interactionTerms[[interaction]][2])
        dat[, names(interactionTerms)[interaction]] <-
          dat[, interactionTerms[[interaction]][1]] *
          dat[, interactionTerms[[interaction]][2]];
      }
      regressionFormula <- paste(regressionFormula, " + ",
                                 paste(names(interactionTerms), collapse = " + "));
    }
    
    modelVar <- sum(coef**2)
    
    if (modelVar < 1) { yvar <- sqrt(1 - modelVar)    
    } else  { yvar <- 0; print("warning: model variance equals or exceeds 1") }
    
    datm <- as.matrix(dat)
    mu <- datm %*% coef
    
    dat$y <- rnorm(repetitions, mu, yvar)
    
    regressionOutcome <- lm(formula(regressionFormula), dat);
    
    output[i,-1] <- tail(summary(regressionOutcome)$coefficients[,4], -1);
    output[i, 1] <- summary(regressionOutcome)$r.squared
  }
  
  colnames(output) <- c("rsq",predictNames,names(interactionTerms))
  return(output);
  
}
