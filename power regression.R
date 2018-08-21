
# Compute power for simple moderation model

require(userfriendlyscience)
require(MASS)
require(lm.beta)
require(plyr)

options(digits = 3)
 

simPower.moderation <- function(samSize = c(50,100,150,200,250,300), 
                                alphalevel = c(0.5, 0.1), 
                                n.predictors = 3,
                                l.interactions = NULL,
                                cor.pred = NULL,
                                bpar = c(.5, .3, .2, .2)) 
      {   
    numrow <- length(alphalevel)*length(samSize)
    result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= (2 + n.predictors + length(l.interactions))))
    colnames(result)[c(1,2)] <- c("N","alpha" )
    result[,"alpha"] <- rep(sort(rep(alphalevel)),length(samSize))
    result[,"N"] <- sort(rep(samSize,length(alphalevel)))
    
    cor.pred <- runif((n.predictors*(n.predictors-1)/2),-0.30,0.30)
   
     for (N in samSize)
     {
       for (a in alphalevel) 
         { 
 
          out <- regrPwrSim(n = N, 
                            predictors = n.predictors, 
                            cor = cor.pred, 
                            betas = bpar, 
                            interactions = l.interactions,
                            sig.level = a, 
                            maxiter = 500, 
                            monteCarlo = TRUE)
             
          result[((result$alpha == a) &  result$N == N),-c(1:2)] <- out
          print(c(N, a))
           
         }      # loop over alpha levels
       }        # loop over sample sizes
    
       colnames(result) <- c("N","alpha",names(out))
       
       class(result) <- "simPower.moderation";
     
       return(result)

} # end function



# test the function

res <- simPower.moderation(samSize = c(50,100,150,200,250,300,350,400),  
                           alphalevel = c(0.05, 0.01), 
                           n.predictors = 3,
                           l.interactions = list(c("x1","x2")),
                           cor.pred = NULL,
                           bpar = c(.5, .3, .2, .2)) 

  
  save(res, file= "result_moderation_alpha05.Rdata")
load("result_moderation_alpha01.Rdata")



## Print tabel of results

require(pander)

pander(res) 


## Plot the results

require(ggplot2)
library(viridis)


res$alpha <- ordered(res$alpha, levels=c(0.05, 0.01))
 
plot.simPower.moderation <- function(x, ...) {
  for (i in 1:(length(x)-2)){
        dat <- data.frame(N=x$N,y=x[[i+2]],alpha=x$alpha)
p <- ggplot(data=dat, aes(y=y, x=N,  colour=alpha)) +
        geom_point(size=2) + geom_line(size=1) +
        geom_hline(yintercept=0.80, linetype="dashed", color = "red") +
        geom_hline(yintercept=0.90, linetype="dashed", color = "blue") +
        scale_y_continuous(breaks=seq(0.10, 1, 0.10), limits = c(0.1,1)) + 
        scale_x_continuous(breaks=seq(50, 500, 50) ) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size = 14) +
        ggtitle(paste0("Power of: ",i, "-th predictor")) +
        theme(plot.title = element_text(size=10, hjust=0)) 
print(p)
}
} 

#ggsave(plot = p,filename="Power moderation alpha05.pdf", width=7, height=5)
       
   plot(x=res)        
