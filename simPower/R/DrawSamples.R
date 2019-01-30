#' Simulation for power of moderation model
#'
#' This function allows you to estimate the power in an linear model with interaction effects, used by simPwr.regr().
#' @param preds The number of predictor variables
#' @param Sigma The correlation matrix of the predictor variables
#' @param coef The regression coefficients of respectively the predictors and (if present) the interaction terms
#' @param predictNames Names of the predictors (usually left NULL, by default called "x1", "x2" etc.)
#' @param dependName Name of the dependent variable (usually left NULL, by default called "y")
#' @param interactionTerms List of vectors of length two, that specifies the two-way interaction terms, e.g. interactionTerms = list(c("x1","x2"), c("x1","x3"))
#' @param repetitions The number of records (sample size) in a single data frame
#' @param samples The number of samples drawn from the specified model (replications)
#' @details Instead of interaction terms also quadratic terms can be specified similar to interaction terms, e.g. interactionTerms = list(c("x1","x1"))
#' @keywords power regression interaction moderation
#' @return A matrix with in the rows representing the samples, and columns the R-squared values and the p-values, respectively 
#' @export
#' @import MASS
#' @examples drawSamples(preds = 3, coef = c(.4,.3.,2), repetitions = 100, samples=100)
drawSamples <- function(preds = NULL,
                        Sigma = NULL,
                        coef = NULL,
                        predictNames = NULL,
                        dependName = NULL,
                        interactionTerms = NULL,
                        repetitions = 100,
                        samples=1000)   {
  
  if (is.null(preds)) stop("The number of predictors must be specified")
  if (is.null(Sigma)) { Sigma = matrix(0, preds, preds) ; diag(Sigma) <- 1 }
  if (is.null(dependName)) dependName <- "y"
  if (is.null(predictNames)) predictNames = paste0("x", 1:preds)
  
  output <- matrix(data=0, nrow=samples, ncol=(preds + length(interactionTerms) + 1))
  mu0 = rep(0, preds)
  
  
  for (i in 1:samples) {
    dat <- data.frame( mvrnorm(n = repetitions, mu = mu0, Sigma = Sigma))
    
    colnames(dat) <- c(predictNames)
    
    regressionFormula <- paste(dependName, "~", paste0(predictNames, collapse=" + "));
    
    
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

