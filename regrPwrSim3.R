 require(pwr)
 require(MASS)
 #require(mosaic)
 
a <- regrPwrSim(100, predictors=4, betas = c(.3,.4,.2,.1), cor=c(.1,.1,.1,.2,.3,.4))

a <- regrPwrSim(n=100, cor=matrix(c( 1, .2, .3, .5,
                                    .2,  1, .2, .3,
                                    .3, .2,  1, .4,
                                    .5, .3, .4,  1), ncol=4, byrow=TRUE),
                     betas = c(.3,.4,.2,.1,.3,.2),
                     interactions = list(c("x1","x2"), c("x1","x3")),
                     maxiter = 100);

### Power simulations for regression analyses using roughly equally
### correlated predictors and one dependent variable

regrPwrSim <- function(n, predictors = NULL, cor = c(.3, .5), betas= NULL,
                       maxiter=100, sig.level=.05, interactions = NULL,
                       digits=2) {
  
  res <- list(intermediate = list(cor = cor),output = list());
  
  if (is.vector(cor) ) {
    if (is.null(predictors)) {
      stop("If you don't supply a correlation matrix in the 'cor' argument, ",
           "you have to specify the number of predictors in the 'predictors' ",
           "argument!");
    }
    if (length(cor) > (predictors*(predictors-1)/2)) {
      warning("You specified too many correlations. Only the first will be used")
      cor <- cor[1:(predictors*(predictors-1)/2)];
    }
    if (length(cor) < (predictors*(predictors-1)/2)) {
      warning("You specified too few correlations. Zero's will be added");
      cor <- c(cor, rep(0, ((predictors*(predictors-1)/2)-length(cor))))
    }
    M <- matrix(1, predictors, predictors)
    M[lower.tri(M, diag = F)] <-  cor 
    M <-  M + t(M) - 1
    res$intermediate$cor <- M
  }
  
  else {
    predictors <- (ncol(cor));
    res$intermediate$cor <- cor
  }

  if (length(betas) != (predictors + length(interactions)) ) {
      stop("The nUmber of beta's that you specified does not equal the sum of the 
             predictors and the interactions terms")
  }
  
  res$input <- as.list(environment())
  
  predictorNames = paste0("x", 1:predictors)
    
  res$intermediate$regressionAnalyses <- drawSamples(preds = predictors,
                                                     Sigma = res$intermediate$cor,
                                                     predictNames = predictorNames,
                                                     dependName = "y",
                                                     repetitions = n, 
                                                     coef = betas,
                                                     samples=maxiter,
                                                     interactionTerms = interactions)
  
  
  res$intermediate$power <- apply(res$intermediate$regressionAnalyses,2,
                                         function(x) { sum(x < sig.level); })/maxiter ;
  
  res$intermediate$CI <- apply(res$intermediate$regressionAnalyses,2,
                                  function(x) { 
                                   binom.test(x=sum(x < sig.level), n=maxiter)$conf.int
                                   }) ;
  
  res$intermediate$rsq <- mean((res$intermediate$regressionAnalyses)[,1])
  
  res$output$regres <- data.frame(regressionPower = res$intermediate$power[-1])
  res$output$rsq <- (res$intermediate$rsq )
  res$output$CI <- data.frame(res$intermediate$CI[,-1] )
  
  
  class(res) <- 'regrPwrSim';
  
  return(res);
  
}

print.regrPwrSim <- function(x, digits=x$input$digits, ...) {
  cat(" Ran ",x$input$maxiter , " regression analyses, each with a sample size of ",
      x$input$n, " and ", x$input$predictors, " predictors,  \n and ",
      length(a$input$interactions)," interaction terms. \n\n",
      " You provided this correlation between predictors matrix:\n\n",sep="");
  print(x$intermediate$cor);
  cat("\n These are regression coefficients you provided:\n\n", sep="");
  print(x$input$betas);
  cat("\n These are the power estimates:\n\n", sep="");
  print(x$output$regres, digits=digits);
  cat("\n These are 95% confidence intervals of the power estimates:\n\n", sep="");
  print(x$output$CI, digits=digits);
}




drawSamples <- function(preds = NULL,
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
  
}   # End drawSamples

