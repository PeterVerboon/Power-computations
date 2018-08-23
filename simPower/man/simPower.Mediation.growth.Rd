\name{simPower.Mediation.growth}
\alias{simPower.Mediation.growth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simPower.Mediation.growth(n = 200, EScond = 0.2, ESmod = 0.2, ESint = 0.2, bpath = c(0.4, 0.3), rho = c(0, 0), error = 1, alpha = 0.05, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{EScond}{
%%     ~~Describe \code{EScond} here~~
}
  \item{ESmod}{
%%     ~~Describe \code{ESmod} here~~
}
  \item{ESint}{
%%     ~~Describe \code{ESint} here~~
}
  \item{bpath}{
%%     ~~Describe \code{bpath} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{error}{
%%     ~~Describe \code{error} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
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
function (n = 200, EScond = 0.2, ESmod = 0.2, ESint = 0.2, bpath = c(0.4, 
    0.3), rho = c(0, 0), error = 1, alpha = 0.05, maxiter = 100) 
{
    input <- c(EScond, ESmod, ESint, bpath, error, alpha, n, 
        maxiter)
    blabel <- "b1"
    for (i in 2:length(bpath)) {
        blabel <- paste0(blabel, ",", "b", i)
    }
    blabel <- (unlist(strsplit(blabel, ",", fixed = TRUE)))
    names(input) <- c("EScondition", "ESmoderator", "ESinteraction", 
        blabel, "error", "alpha", "sampleSize", "maxIter")
    ES <- c(EScond, ESmod, ESint, bpath, bpath * EScond)
    model0 <- buildSimModel(EScond = EScond, ESmod = ESmod, ESint = ESint, 
        bpath = bpath, model = "growth")
    model <- buildSimModel(EScond = "a1", ESmod = "a2", ESint = "a3", 
        bpath = c("b1", "b2"), model = "growth")
    res <- matrix(data = 0, nrow = maxiter, ncol = 18)
    colnames(res) <- c("condition", "power cond", "moderation", 
        "power mod", "interaction", "power int", "effect li", 
        "power li", "effect ls", "power ls", "indirect i", "power ind_li", 
        "indirect s", "power ind_ls", "cfi", "tli", "rmsea", 
        "srmr")
    pb <- txtProgressBar(min = 0, max = maxiter, style = 3)
    for (i in 1:maxiter) {
        data <- simulateData(model0, sample.nobs = n)
        data2 <- simDataSet(n = n, varNames = c("mediator", "y1", 
            "y2", "y3", "y4", "condition", "moderator", "interaction"), 
            means = 0, sds = error, ranges = list(c(condition = c(-3, 
                3))), silent = TRUE, seed = NULL, correlations = c(rho[1], 
                rho[2]))
        data3 <- data + data2
        data3$condition <- as.numeric(cut(data3$condition, breaks = 2))
        result <- sem(model, data3)
        res[i, 1] <- parameterEstimates(result)[1, 5]
        res[i, 2] <- parameterEstimates(result)[1, 8]
        res[i, 3] <- parameterEstimates(result)[2, 5]
        res[i, 4] <- parameterEstimates(result)[2, 8]
        res[i, 5] <- parameterEstimates(result)[3, 5]
        res[i, 6] <- parameterEstimates(result)[3, 8]
        res[i, 7] <- parameterEstimates(result)[4, 5]
        res[i, 8] <- parameterEstimates(result)[4, 8]
        res[i, 9] <- parameterEstimates(result)[5, 5]
        res[i, 10] <- parameterEstimates(result)[5, 8]
        res[i, 11] <- filter(parameterEstimates(result), lhs \%in\% 
            c("ind1"))[, 5]
        res[i, 12] <- filter(parameterEstimates(result), lhs \%in\% 
            c("ind1"))[, 8]
        res[i, 13] <- filter(parameterEstimates(result), lhs \%in\% 
            c("ind2"))[, 5]
        res[i, 14] <- filter(parameterEstimates(result), lhs \%in\% 
            c("ind2"))[, 8]
        res[i, c(15:18)] <- fitmeasures(result)[c("cfi", "tli", 
            "rmsea", "srmr")]
        setTxtProgressBar(pb, i)
    }
    sigs <- res[, c(2, 4, 6, 8, 10, 12, 14)] < alpha
    power <- apply(sigs, 2, mean)
    bias <- ES - apply(res[, c(1, 3, 5, 7, 9, 11, 13)], 2, mean)
    output <- list(power = power, bias = bias, raw = res, input = input)
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
