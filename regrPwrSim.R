# require(pwr)
# require(MASS)
# require(mosaic)
regrPwrSim(100, predictors=3)

regrPwrSim(100, cor=matrix(c( 1, .2, .3, .5,
                             .2,  1, .2, .3,
                             .3, .2,  1, .4,
                             .5, .3, .4,  1), ncol=4, byrow=TRUE),
           predictorNames = c('abuseSeverity',
                              'autoMutilation',
                              'suicideAttempts'),
           dependentName = 'treatmentSatisfaction');

### Power simulations for regression analyses using roughly equally
### correlated predictors and one dependent variable

regrPwrSim <- function(n, predictors=NULL, cor = c(.3, .5), 
                       predictorNames = paste0("predictor_", 1:predictors),
                       dependentName = "dependent", 
                       samples=100, sig.level=.05, interactions = NULL,
                       digits=2) {
  
  res <- list(intermediate = list(cor = cor),
              output = list());
  
  if (is.vector(cor) && (length(cor) == 2)) {
    if (is.null(predictors)) {
      stop("If you don't supply a correlation matrix in the 'cor' argument, ",
           "you have to specify the number of predictors in the 'predictors' ",
           "argument!");
    }
    res$intermediate$cor <- data.frame(matrix(cor[1], nrow=predictors,
                                              ncol=predictors));
    res$intermediate$cor[, predictors+1] <- cor[2];
    res$intermediate$cor[predictors+1, ] <- cor[2];
    diag(res$intermediate$cor) <- 1;
  }
  else if (!is.matrix(cor) && (ncol(cor) == nrow(cor)) &&
             (ncol(cor)== predictors)) {
    stop("The 'cor' argument has to be either a vector with two elements, ",
         "specifying the correlation between the predictors with each other ",
         "and the dependent variable, or a matrix representing the ",
         "correlations.");
  }
  else {
    predictors <- (ncol(cor) - 1);
  }

  res$input <- as.list(environment());
    
  colnames(res$intermediate$cor) <- rownames(res$intermediate$cor) <- c(predictorNames, dependentName);
  
  res$intermediate$predictorCorrelations <-
    res$intermediate$cor[-nrow(res$intermediate$cor), -ncol(res$intermediate$cor)];
  
  res$intermediate$meanPredictorsCor <-
    mean(res$intermediate$predictorCorrelations[lower.tri(res$intermediate$predictorCorrelations)]);
  
  res$intermediate$regressionAnalyses <-
    do(samples) * function(mu = rep(0, (predictors+1)),
                           Sigma = res$intermediate$cor,
                           predictNames = predictorNames,
                           dependName = dependentName,
                           repetitions = n,
                           interactionTerms = interactions) {
      dat <- data.frame(mvrnorm(n = repetitions, mu = mu, Sigma = Sigma));
      names(dat) <- c(predictNames, dependName);
      regressionFormula <- paste(dependName, "~",
                                  paste0(predictNames, collapse=" + "));
      if (!is.null(interactionTerms)) {
        for (interaction in 1:length(interactionTerms)) {
          if (!(min(interactionTerms[[interaction]]) < 1 ||
                  !(min(interactionTerms[[interaction]]) > predictors))) {
            stop("The lowest number in the list provided as argument 'interactions' ",
                 "is lower than 1 or higher than the number of predictors!");
          }
          if (!(max(interactionTerms[[interaction]]) > predictors)) {
            stop("The highest number in the list provided as argument 'interactions' ",
                 "is higher than the number of predictors!");
          }
          dat[, names(interactionTerms[[interaction]])] <-
            dat[, interactionTerms[[interaction]][1]] *
            dat[, interactionTerms[[interaction]][2]];
        }
        regressionFormula <- paste(regressionFormula,
                                   paste(names(interactionTerms), collapse = " + "));
      }
      regressionOutcome <- lm(formula(regressionFormula), dat);
      res <- tail(summary(regressionOutcome)$coefficients[,4], -1);
      return(res);
    };
  
  res$intermediate$significant <- apply(res$intermediate$regressionAnalyses,
                                        2, function(x) { sum(x < sig.level); }) / samples;
  
  res$intermediate$bivariatePower <- sapply(res$intermediate$cor[-nrow(res$intermediate$cor), ncol(res$intermediate$cor)],
                                            function(x, n) {
                                              return(pwr.r.test(r=x, n=n, sig.level=sig.level)$power);
                                            }, n=n);
  
  res$intermediate$bivariatePower.adjusted <- sapply(res$intermediate$cor[-nrow(res$intermediate$cor), ncol(res$intermediate$cor)],
                                            function(x, n) {
                                              return(pwr.r.test(r=x, n=n, sig.level=sig.level/predictors)$power);
                                            }, n=n);
  
  
  res$output$dat <- data.frame(biVarCorrelation = res$intermediate$cor[-nrow(res$intermediate$cor), ncol(res$intermediate$cor)],
                               regressionPower = res$intermediate$significant,
                               correlationPower = res$intermediate$bivariatePower,
                               corPower.adjusted = res$intermediate$bivariatePower.adjusted);
  
  row.names(res$output$dat) <- c(predictorNames);
  
  class(res) <- 'regrPwrSim';
  
  return(res);
  
}

print.regrPwrSim <- function(x, digits=x$input$digits, ...) {
  cat("Ran ", x$input$samples, " regression analyses, each with a sample size of ",
      x$input$n, " and ", x$input$predictors, " predictors, with an average ",
      "correlation between predictors of ",
      x$intermediate$meanPredictorsCor, " and an average correlation between ",
      "predictors and the criterion of ", mean(x$output$dat$biVarCorrelation),
      ". You provided this correlation matrix:\n\n",
      sep="");
  print(x$intermediate$cor);
  cat("\nThese are the bivariate correlations, achieved power, and the power ",
      "for the bivariate correlations:\n\n", sep="");
  print(x$output$dat, digits=digits);
}
