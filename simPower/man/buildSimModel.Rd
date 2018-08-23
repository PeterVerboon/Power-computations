\name{buildSimModel}
\alias{buildSimModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
buildSimModel(EScond = 0.2, ESmod = 0.2, ESint = 0.2, bpath = c(0.4, 0.3), model = "growth")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (EScond = 0.2, ESmod = 0.2, ESint = 0.2, bpath = c(0.4, 
    0.3), model = "growth") 
{
    modela1 <- paste0("mediator", " ~ ", EScond, "*", "condition", 
        collapse = " \n ")
    modela2 <- paste0("mediator", " ~ ", ESmod, "*", "moderator", 
        collapse = " \n ")
    modela3 <- paste0("mediator", " ~ ", ESint, "*", "interaction", 
        collapse = " \n ")
    if (model == "growth") {
        modelb1 <- paste0("li", " ~ ", bpath[1], "*", "mediator", 
            collapse = " \n ")
        modelb2 <- paste0("ls", " ~ ", bpath[2], "*", "mediator", 
            collapse = " \n ")
        modelb <- paste0(modelb1, "\n ", modelb2)
        modelc1 <- ("li =~ 1*y1 + 1*y2 + 1*y3 + 1*y4")
        modelc2 <- ("ls =~ 3*y1 + 2*y2 + 1*y3 + 0*y4")
    }
    else {
        modelb <- paste0("y1", " ~ ", bpath[1], "*", "mediator", 
            collapse = " \n ")
        for (i in 2:length(bpath)) {
            value <- paste0("y", i, " ~ ", bpath[i], "*", "mediator", 
                collapse = " \n ")
            modelb <- paste0(modelb, " ;  ", value, " ; ", collapse = " \n ")
        }
        modelc1 <- "  "
        modelc2 <- "  "
    }
    modelind1 <- paste0("ind1 := ", EScond, "*", bpath[1], collapse = " \n ")
    modelind2 <- ifelse(model != "growth", "  ", paste0("ind2 := ", 
        EScond, "*", bpath[2], collapse = " \n "))
    model <- paste0(modela1, " \n ", modela2, " \n ", modela3, 
        " \n ", modelb, " \n ", modelc1, " \n ", modelc2, " \n ", 
        modelind1, " \n ", modelind2)
    return(model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
