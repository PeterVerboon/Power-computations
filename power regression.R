
# Function simPower.moderation() calls repeatedly (for vector of sample sizes and alpha levels) 
# the function regrPwrSim().
# Function regrPwrSim() uses the function drawSamples().


# Compute power for  moderation model
#

require(userfriendlyscience)
require(MASS)  # mvrnorm
require(lm.beta)
require(plyr)

options(digits = 3)
 

simPower.moderation <- function(samSize = c(50,100,150,200,250,300), 
                                alphalevel = c(0.5, 0.1), 
                                n.predictors = 3,
                                l.interactions = NULL,
                                cor.pred = NULL,
                                bpar = c(.5, .3, .2, .2),
                                niter = 100) 
      {  
  
    output <- list();
    output$input <- as.list(environment());
  
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
                            maxiter = niter, 
                            monteCarlo = TRUE)
             
          result[((result$alpha == a) &  result$N == N),-c(1:2)] <- out
          print(c(N, a))
           
         }      # loop over alpha levels
       }        # loop over sample sizes
    
       colnames(result) <- c("N","alpha",names(out))
       
       output$result <- result
       
       class(output) <- "simPower.moderation";
       
       return(output)

} # end function



# test the function

res <- simPower.moderation(samSize = c(50,100,150,200),  
                           alphalevel = c(0.05, 0.01), 
                           n.predictors = 3,
                           l.interactions = list(c("x1","x1")),
                           cor.pred = NULL,
                           bpar = c(.5, .3, .2, .2),
                           niter = 100) 

  
saveRDS(res, file= "result_moderation_alpha05.Rds")
load("result_moderation_alpha01.Rds")




## Plot the results

require(ggplot2)
library(viridis)
require(broom)


res$alpha <- ordered(res$alpha, levels=c(0.05, 0.01))
 
library(reshape)


plot.simPower.moderation <- function(x, pval=0.05, ...) {
        
        dat <- res$result
        dat <- dat[dat[,"alpha"]== pval,]
        dat[,"alpha"] == pval
        mdat <- melt(dat, id=c("N","alpha"))
        levels(mdat$variable) <-  paste0(levels(mdat$variable)," = ",res$input$bpar)
        mdat <- rename(mdat, c(variable = "effectSize")) 
        
p <- ggplot(data=mdat, aes(y=value, x=N,  colour=effectSize, group = effectSize)) + geom_point(size=2) +
        geom_line(size=.5) + labs( y = "Estimated power") +
        geom_hline(yintercept=0.80, linetype="dashed", color = "red") +
        geom_hline(yintercept=0.90, linetype="dashed", color = "blue") +
        scale_y_continuous(breaks=seq(0.10, 1, 0.10), limits = c(0.1,1)) + 
        scale_x_continuous(breaks=seq(50, 500, 50) ) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size = 12) +
        ggtitle(paste0("alpha = ", pval))
        
print(p)
}


#ggsave(plot = p,filename="Power moderation alpha05.pdf", width=7, height=5)
       
   plot(res)        
plot(res, pval = .01)
