
#' Simulation for power of latent growth model
#'
#' This function computes the power curves
#' of a latent growth model with k measurements, one mediation effect and one interaction effect.
#' @param samSizes vector of sample size
#' @param esSizes vector of effect sizes of the condition effect
#' @param ESmod effect size moderator
#' @param ESint effect size interaction term
#' @param bpath vector of regression coefficients from mediator to latent intercept and slope respectively
#' @param ndepend number of dependent variables
#' @param rho vector with minimum and maximum value of randomly selected corrrelation between the variables
#' @param error standard deviations of the random error added to the data 
#' @param alpha alpha level
#' @param niter number of iterations
#' @keywords SEM latent growth mediation
#' @export
#' @import MASS
#' @import lavaan
#' @return List with the following elements
#' @return power: estimate for all effects
#' @return bias: estimate for all effects
#' @return raw: the raw results of the simulation
#' @return input: the input paramters used in the simulation
#' @examples
#' simPwr.growth()
simPwr.lgm <- function(samSizes, esSizes, ndepend = 4, ESmod = .3, ESint = .2, 
                       bpath = c(0.5,-0.1), rho = c(0.0, 0.1), error = 1, alpha = .05, maxiter = 1000 ) {

out <- matrix(0,nrow=(length(samSizes)*length(esSizes)),ncol=9)
colnames(out) <- c("ES","N","condition","moderation", "interaction","effect_li", "effect_ls", "indirect_i","indirect_s" )
out[,1] <- rep(esSizes,each = length(samSizes))
out[,2] <- rep(samSizes,length(esSizes))

for (es in esSizes) {
    for (n in samSizes) {
    
    res <- simPwr.growth(n=n,
                         EScond = es,
                         ESmod = ESmod,
                         ESint = Esint,
                         bpath = bpath,
                         ndepend = ndepend,
                         rho = rho,
                         error = error,
                         alpha = alpha,
                         maxiter = maxiter) 
    
    out[out[,"N"] == n & out[,"ES"] == es,-c(1,2)] <- res$power
  }
}

class(out) <- "simPwr.lgm"
return(out)

}









## plot the results

out1 <- data.frame(out)
out1$EffectSize <- as.factor(out1$ES)

p <- ggplot(data=out1, aes(y=indirect_i, x=N,  colour=EffectSize)) +
  geom_point(size=2) + geom_line(size=1) +
  geom_hline(yintercept=0.80, linetype="dashed", color = "red") +
  geom_hline(yintercept=0.90, linetype="dashed", color = "blue") +
  scale_y_continuous(breaks=seq(0.10, 1, 0.10)) + scale_x_continuous(breaks=seq(100, 600, 50) ) +
  scale_color_viridis(discrete=TRUE) +
  theme_bw(base_size = 14) +
  ggtitle("Power of moderated mediation latent growth model, alpha=0.05") +
  theme(plot.title = element_text(size=10, hjust=0)) 
p

# save the plot as pdf
ggsave(plot=p,filename="Modmedlatent_alpha05_cond.pdf", width=7, height=5)
