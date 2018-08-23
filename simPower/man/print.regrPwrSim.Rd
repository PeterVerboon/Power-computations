\name{print.regrPwrSim}
\alias{print.regrPwrSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
print.regrPwrSim(x, digits = x$input$digits, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, digits = x$input$digits, ...) 
{
    cat(" Ran ", x$input$maxiter, " regression analyses, each with a sample size of ", 
        x$input$n, " and ", x$input$predictors, " predictors,  \n and ", 
        length(a$input$interactions), " interaction terms. \n\n", 
        " You provided this correlation between predictors matrix:\n\n", 
        sep = "")
    print(x$intermediate$cor)
    cat("\n These are regression coefficients you provided:\n\n", 
        sep = "")
    print(x$input$betas)
    cat("\n These are the power estimates:\n\n", sep = "")
    print(x$output$regres, digits = digits)
    cat("\n These are 95\% confidence intervals of the power estimates:\n\n", 
        sep = "")
    print(x$output$CI, digits = digits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
