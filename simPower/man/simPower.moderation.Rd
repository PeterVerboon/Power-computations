\name{simPower.moderation}
\alias{simPower.moderation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simPower.moderation(samSize = c(50, 100, 150, 200, 250, 300), alphalevel = c(0.5, 0.1), n.predictors = 3, l.interactions = NULL, cor.pred = NULL, bpar = c(0.5, 0.3, 0.2, 0.2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samSize}{
%%     ~~Describe \code{samSize} here~~
}
  \item{alphalevel}{
%%     ~~Describe \code{alphalevel} here~~
}
  \item{n.predictors}{
%%     ~~Describe \code{n.predictors} here~~
}
  \item{l.interactions}{
%%     ~~Describe \code{l.interactions} here~~
}
  \item{cor.pred}{
%%     ~~Describe \code{cor.pred} here~~
}
  \item{bpar}{
%%     ~~Describe \code{bpar} here~~
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
function (samSize = c(50, 100, 150, 200, 250, 300), alphalevel = c(0.5, 
    0.1), n.predictors = 3, l.interactions = NULL, cor.pred = NULL, 
    bpar = c(0.5, 0.3, 0.2, 0.2)) 
{
    numrow <- length(alphalevel) * length(samSize)
    result <- as.data.frame(matrix(data = 0, nrow = numrow, ncol = (2 + 
        n.predictors + length(l.interactions))))
    colnames(result)[c(1, 2)] <- c("N", "alpha")
    result[, "alpha"] <- rep(sort(rep(alphalevel)), length(samSize))
    result[, "N"] <- sort(rep(samSize, length(alphalevel)))
    cor.pred <- runif((n.predictors * (n.predictors - 1)/2), 
        -0.3, 0.3)
    for (N in samSize) {
        for (a in alphalevel) {
            out <- regrPwrSim(n = N, predictors = n.predictors, 
                cor = cor.pred, betas = bpar, interactions = l.interactions, 
                sig.level = a, maxiter = 500, monteCarlo = TRUE)
            result[((result$alpha == a) & result$N == N), -c(1:2)] <- out
            print(c(N, a))
        }
    }
    colnames(result) <- c("N", "alpha", names(out))
    class(result) <- "simPower.moderation"
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
