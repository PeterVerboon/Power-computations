

### Power simulations for regression analyses using correlated predictors,
### interaction effects (optional) and one dependent variable.

 #' Simulation for power moderation for a given sample size
 #'
 #' This function allows you to estimate the power in an linear model with interaction effects.
 #' @param n The sample size
 #' @param predictors A number of predictor variables
 #' @param cor A vector or matrix containing the correlation bewtween the predictors, if NULL random correlations will be used
 #' @param betas A vector regression coeffcients (effect sizes)
 #' @param interactions A list of vectors containing a pair of variables
 #' @param sig.level Type I error level for power analysis
 #' @param maxiter The number of independent replications (iterations, samples)
 #' @param monteCarlo Logical indicating whether function is used in a Monte Carlo experiment (TRUE)
 #' @keywords power regression interaction moderation
 #' @export
 #' @import MASS
 #' @examples SimPwr.regr(n=100, predictors = 3, betas = c(.4,.3.,.2), maxiter= 1000)
regrPwrSim <- function(n = 100,
                       predictors = NULL,
                       cor = NULL,
                       betas = NULL,
                       interactions = NULL,
                       sig.level=.05,
                       maxiter=1000,
                       monteCarlo = FALSE) {

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
  res$intermediate$CI <- apply(res$intermediate$regressionAnalyses,2,
                               function(x) {
                                 binom.test(x=sum(x < sig.level), n=maxiter)$conf.int
                               }) ;
  rownames(res$intermediate$CI) <- c(" low_CI ", "high_CI ")


  if (monteCarlo == FALSE) {

     res$input <- as.list(environment())

     res$intermediate$rsq <- mean((res$intermediate$regressionAnalyses)[,1])

     res$output$regres <- data.frame(regressionPower = res$intermediate$power[-1])
     res$output$rsq <- (res$intermediate$rsq )
     res$output$CI <- data.frame(res$intermediate$CI[,-1] )

     class(res) <- 'simPwr.regr';

  }

  if (monteCarlo == TRUE) {
     res <- rbind(t(res$intermediate$power[-1]),res$intermediate$CI[,-1])
  }

  return(res)
}




#' Print method for function regrPwrSim
#'
#' This function prints the results of a power analyses of the linear model
#' @param x simPwr.regr object
#' @examples print(x)
print.simPwr.regr <- function(x, digits=x$input$digits, ...) {
  cat(" Ran ",x$input$maxiter , " regression analyses, each with a sample size of ",
      x$input$n, " and ", x$input$predictors, " predictors,  \n and ",
      length(x$input$interactions)," interaction terms. \n\n",
      " You provided this correlation between predictors matrix:\n\n",sep="");
  print(x$intermediate$cor);
  cat("\n These are regression coefficients you provided:\n\n", sep="");
  print(x$input$betas);
  cat("\n These are the power estimates:\n\n", sep="");
  print(x$output$regres, digits=digits);
  cat("\n These are 95% confidence intervals of the power estimates:\n\n", sep="");
  print(x$output$CI, digits=digits);
}





