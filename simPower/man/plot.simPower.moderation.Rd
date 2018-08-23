\name{plot.simPower.moderation}
\alias{plot.simPower.moderation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
The function plots the power curves after a simulation. If more than one alpha levels are used in the simulation,
the alpha level can be chosen. 
}
\usage{
plot.simPower.moderation(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class simPower.moderation 
}
  \item{\dots}{
The value of alpha (p value) can be specified here.
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
Peter Verboon
}
\note{
It needs the packages ggplot2, viridis and reshape
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
function (x, ...) 
{
    for (i in 1:(length(x) - 2)) {
        dat <- data.frame(N = x$N, y = x[[i + 2]], alpha = x$alpha)
        p <- ggplot(data = dat, aes(y = y, x = N, colour = alpha)) + 
            geom_point(size = 2) + geom_line(size = 1) + geom_hline(yintercept = 0.8, 
            linetype = "dashed", color = "red") + geom_hline(yintercept = 0.9, 
            linetype = "dashed", color = "blue") + scale_y_continuous(breaks = seq(0.1, 
            1, 0.1), limits = c(0.1, 1)) + scale_x_continuous(breaks = seq(50, 
            500, 50)) + scale_color_viridis(discrete = TRUE) + 
            theme_bw(base_size = 14) + ggtitle(paste0("Power of: ", 
            i, "-th predictor")) + theme(plot.title = element_text(size = 10, 
            hjust = 0))
        print(p)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
