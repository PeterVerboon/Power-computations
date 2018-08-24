 #require(pwr)
 require(MASS)
 #require(psych)

 # Function drawSamples
 
 #' Simulation for power of moderation model 
 #'
 #' This function allows you to estimate the power in an linear model with interaction effects.
 #' @param preds The number of predictor variables
 #' @param Sigma The correlation matrix of the predictor variables
 #' @param coef The regression coefficients of respectively the predictors and (if present) the interaction terms
 #' @param predictNames Names of the predictors (usually left NULL, by default called "x1", "x2" etc.)
 #' @param dependName Name of the dependent variable (usually left NULL, by default called "y")
 #' @param interactionTerms List of vectors of length two, that specifies the two-way interaction terms, e.g. interactionTerms = list(c("x1","x2"), c("x1","x3")) 
 #' @param repetitions The number of records (sample size) in single data frame
 #' @param samples The number of samples drawn from the specified model 
 #' @details Instead of interaction terms also quadratic terms can be specified similar to interaction terms, e.g. interactionTerms = list(c("x1","x1"))
 #' @keywords power regression interaction moderation
 #' @value A matrix with in the rows representing the samples, and columns the R-squared values and the p-values, respectively 
 #' @export
 #' @import MASS
 #' @examples
 #' drawSamples(preds = 3, coef = c(.4,.3.,2), repetitions = 100, samples=100)
 drawSamples <- function(preds = NULL,
                         Sigma = NULL,
                         coef = NULL,
                         predictNames = NULL,
                         dependName = NULL,
                         interactionTerms = NULL,
                         repetitions = 100, samples=100)   {
   
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
   
 }   # End drawSamples
 
 
### Power simulations for regression analyses using correlated predictors,
### interaction effects (optional) and one dependent variable.

 #' Simulation for power moderation function for one sample size
 #'
 #' This function allows you to estimate the power in an linear model with interaction effects.
 #' @param predictors A number of predictor variables
 #' @param cor A vector or matrix containing the correlation bewtween the predictors, if NULL random correlations will be used
 #' @param betas A vector regression coeffcients (effect sizes) 
 #' @param interactions A list of vectors containing a pair of variables
 #' @param sig.level Type I error level for power analysis  
 #' @keywords power regression interaction moderation
 #' @export
 #' @import MASS
 #' @examples
 #' regrPwrSim(n=100, predictors = 3, betas = c(.4,.3.,.2), maxiter= 1000)
regrPwrSim <- function(n = 100, 
                       predictors = NULL, 
                       cor = NULL, 
                       betas = NULL,
                       interactions = NULL,
                       sig.level=.05,
                       maxiter=1000,
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
  
     class(res) <- 'regrPwrSim';
  }
  
  if (monteCarlo == TRUE) {
     res <- rbind(t(res$intermediate$power[-1]),res$intermediate$CI[,-1])
  }
  
  return(res)
}


# test function

# a <- regrPwrSim(200, predictors=3, cor=c(.0,.0,.0), betas=c(.1,.2,.3), maxiter=1000, monteCarlo = FALSE)
# 
# 
# a <- regrPwrSim(n=100, cor=matrix(c( 1, .2, .3, .5,
#                                      .2,  1, .2, .3,
#                                      .3, .2,  1, .4,
#                                      .5, .3, .4,  1), ncol=4, byrow=TRUE),
#                 betas = c(.3,.4,.2,.1,.3,.2),
#                 interactions = list(c("x1","x2"), c("x1","x3")),
#                 maxiter = 100,
#                 monteCarlo = FALSE);


# print method

print.regrPwrSim <- function(x, digits=x$input$digits, ...) {
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





