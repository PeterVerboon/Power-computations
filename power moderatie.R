
# Compute power for simple moderation model

require(userfriendlyscience)
require(MASS)
require(lm.beta)
require(plyr)

options(digits = 3)
 

simPower.moderation <- function(samSize = c(50,100,150,200,250,300), 
                                errlevel = c(1,3,9), 
                                rholevel = c(0,.3,.5,.8),
                                bpar = c(.5, .3, .2),
                                alpha = 0.05,
                                rep = 1000) 
      {   
    numrow <- length(errlevel)*length(rholevel)*length(samSize)
    result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= 8))
    colnames(result) <- c("N","e","rho", "rsq", "b1","b2","b3","power")
    result[,"e"] <- rep(sort(rep(errlevel,length(rholevel))),length(samSize))
    result[,"rho"] <- rep(rholevel,length(errlevel)*length(samSize))
    result[,"N"] <- sort(rep(samSize,length(errlevel)*length(rholevel)))
    
    # Make relative effects sum to one
    b1 <- bpar[1]/sum(bpar)
    b2 <- bpar[2]/sum(bpar)
    b3 <- bpar[3]/sum(bpar)
   
     for (N in samSize)
     {
       for (e in errlevel) 
         { 
         for (rho in rholevel)
           {
           sigma = matrix(c(1,rho,rho,1), nrow=2, ncol=2)
           
           out <- matrix(0,nrow=rep,ncol=5)
           
           for (i in 1:rep) 
             {
  
      a <- as.data.frame(mvrnorm(n = N, Sigma=sigma, mu=c(0,0)))
      colnames(a)  <- c("x","z")
      a$xz <- a$x*a$z

      error <-  rnorm(N,0,sqrt(e))   # the amount of error determines the ES

      a$y <- sqrt(b1)*a$x + sqrt(b2)*a$z + sqrt(b3)*a$xz + error

      res <- lm( y ~ x + z + xz, data=a) 
      res <- lm.beta(res)
      sig <- summary(res)$coefficients[4,5] < alpha
      rsq <- summary(res)$r.squared
 
      out[i,] <- c(rsq,  res$standardized.coefficients[-1], sig)
  
            }   # loop over replications
           
           
           result[((result$e == e) & (result$rho == rho) & result$N == N),c(4:8)] <- apply(out,2,mean)
           print(c(N, e,rho))
           
          }     # loop over correlations between predictors
         }      # loop over error levels
       }        # loop over sample sizes
    
         return(result)

} # end function



# test the function

res <- simPower.moderation(samSize = c(50,100,150,200,250,300,350,400,450,500),  
                           errlevel = c(1,3,9), 
                           rholevel = c(0.3), 
                           betas = c(.5,.3,.2),
                           alpha = 0.05, 
                           rep=1000) 

  
 save(res, file= "result_moderation_alpha05.Rdata")
load("result_moderation_alpha01.Rdata")



## Print tabel of results

require(pander)

pander(res) 


## Plot the results

require(ggplot2)
library(viridis)


res$effectSize <- ordered(res$e, levels=c(1,3,9), labels= c("0.32","0.22","0.14"))
 
 
p <- ggplot(data=res, aes(y=power, x=N,  colour=effectSize)) +
        geom_point(size=2) + geom_line(size=1) +
        geom_hline(yintercept=0.80, linetype="dashed", color = "red") +
        geom_hline(yintercept=0.90, linetype="dashed", color = "blue") +
        scale_y_continuous(breaks=seq(0.10, 1, 0.10)) + scale_x_continuous(breaks=seq(50, 500, 50) ) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size = 14) +
        ggtitle("Power of interaction term, alpha=0.01") +
        theme(plot.title = element_text(size=10, hjust=0)) 
p
  

ggsave(plot = p,filename="Power moderation alpha05.pdf", width=7, height=5)
       

           
