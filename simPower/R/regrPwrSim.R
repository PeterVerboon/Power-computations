regrPwrSim <-
function(n = 100, 
                       predictors = NULL, 
                       cor = NULL, 
                       betas = NULL,
                       interactions = NULL,
                       sig.level=.05,
                       maxiter=100,
                       monteCarlo = FALSE,
                       digits=2) {
  
  res <- list(intermediate = list(cor = cor),output = list());
  
  if (is.null(cor)) { 
    print("You did not specify correlations.They are randomly selected from from -0.50 to 0.50")
    cor <- runif((predictors*(predictors-1)/2),-0.50,0.50)
  }
  
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
    if (!isSymmetric(cor)) {
      stop("The specified correlation matrix is not symmetric") 
    }
    if (sum(eigen(cor, only.values=TRUE)$values > 0) < ncol(cor)) {
      warning("The specified correlation matrix is not positive definite.
              It is smoothed to be positive defintite") 
      cor <- cor.smooth(cor)
    }
    predictors <- (ncol(cor));
    res$intermediate$cor <- cor
  }
  
  
  if (is.null(betas)) {
    warning("No regression coefficients were specified. 
            They are randomly selected from Normal distribution, mean = 0.00, sd = 0.25")
    betas <- rnorm((predictors + length(interactions)),0,.25)
  }
  if (length(betas) != (predictors + length(interactions)) ) {
      stop("The number of beta's that you specified does not equal the sum of the 
             predictors and the interactions terms")
  }
  
   
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
  
  if (monteCarlo == FALSE) {
    
     res$input <- as.list(environment())
  
     res$intermediate$CI <- apply(res$intermediate$regressionAnalyses,2,
                                  function(x) { 
                                   binom.test(x=sum(x < sig.level), n=maxiter)$conf.int
                                   }) ;
  
     res$intermediate$rsq <- mean((res$intermediate$regressionAnalyses)[,1])
  
     res$output$regres <- data.frame(regressionPower = res$intermediate$power[-1])
     res$output$rsq <- (res$intermediate$rsq )
     res$output$CI <- data.frame(res$intermediate$CI[,-1] )
  
     class(res) <- 'regrPwrSim';
  }
  
  if (monteCarlo == TRUE) {
     res <- unlist(res$intermediate$power[-1])
  }
  
  return(res)
}
