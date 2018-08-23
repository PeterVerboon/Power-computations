\name{drawSamples}
\alias{drawSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
This functions draws samples using a linear regression design with optional interaction terms. 
}
\usage{
drawSamples(preds = NULL, Sigma = NULL, coef = NULL, predictNames = NULL, 
            interactionTerms = NULL,dependName = "y", 
            repetitions = 100, samples = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{preds}{
The number of predictors in the model
}
  \item{Sigma}{
The correlation matrix of the predictor variables
}
  \item{coef}{
The regression coefficients of respectively the predictors and (if present) the interaction terms
}
  \item{predictNames}{
Names of the predictors (usually left NULL, by default called "x1", "x2" etc.)
}
  \item{dependName}{
Name of the dependent variable (usually left NULL, by default called "y")
}
  \item{repetitions}{
The number of records (sample size) in single data frame
}
  \item{samples}{
The number of samples drawn from the specified model
}
  \item{interactionTerms}{
List of vectors of length two, that specifies the two-way interaction terms, 
e.g. interactionTerms = list(c("x1","x2"), c("x1","x3")) 
}
}
\details{
Instead of interaction terms also quadratic terms can be specified similar to interaction terms, 
e.g. interactionTerms = list(c("x1","x1"))
}
\value{
A matrix with in the rows representing the samples, and in the first column the R-squared values from each analysis 
and in the other columns the p-values corresponding to the regression terms in the model.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (preds = NULL, Sigma = NULL, coef = NULL, predictNames = NULL, 
    dependName = "y", repetitions = 100, samples = 100, interactionTerms = NULL) 
{
    output <- matrix(data = 0, nrow = samples, ncol = (preds + 
        length(interactionTerms) + 1))
    mu0 = rep(0, preds)
    for (i in 1:samples) {
        dat <- data.frame(mvrnorm(n = repetitions, mu = mu0, 
            Sigma = Sigma))
        colnames(dat) <- c(predictNames)
        regressionFormula <- paste(dependName, "~", paste0(predictNames, 
            collapse = " + "))
        if (!is.null(interactionTerms)) {
            for (interaction in 1:length(interactionTerms)) {
                names(interactionTerms)[interaction] <- paste0(interactionTerms[[interaction]][1], 
                  interactionTerms[[interaction]][2])
                dat[, names(interactionTerms)[interaction]] <- dat[, 
                  interactionTerms[[interaction]][1]] * dat[, 
                  interactionTerms[[interaction]][2]]
            }
            regressionFormula <- paste(regressionFormula, " + ", 
                paste(names(interactionTerms), collapse = " + "))
        }
        modelVar <- sum(coef^2)
        if (modelVar < 1) {
            yvar <- sqrt(1 - modelVar)
        }
        else {
            yvar <- 0
            print("warning: model variance equals or exceeds 1")
        }
        datm <- as.matrix(dat)
        mu <- datm \%*\% coef
        dat$y <- rnorm(repetitions, mu, yvar)
        regressionOutcome <- lm(formula(regressionFormula), dat)
        output[i, -1] <- tail(summary(regressionOutcome)$coefficients[, 
            4], -1)
        output[i, 1] <- summary(regressionOutcome)$r.squared
    }
    colnames(output) <- c("rsq", predictNames, names(interactionTerms))
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
