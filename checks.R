

# DEFINE FUNCTION
# Assume x and z ~ N(0,1)
# check expected effect sizes using formulas from power moderation paper
# arlevel is relative effect of auto-regression or other covariate uncorrelated with other predictors
# rholevel is correlation between the predictors (x and z)
# errlevel is error level
# bpar provides the relative effects of x, z and xz

 checkPars  <- function(rholevel = c(0.0, 0.3, 0.8),
                        arlevel = c(0,0.5),
                        errlevel = c(0,1,3,9),
                        bpar = c(.5, .3, .2),
                        relative = TRUE)
   {
   
   numrow <- length(rholevel)*length(arlevel)*length(errlevel)
   result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= 9))
   colnames(result) <- c("AR","error","rho","b1","b2","b3","b4","rsq","f2")
   result[,"error"] <- rep(sort(rep(errlevel,length(rholevel))),length(arlevel))
   result[,"rho"] <- rep(rholevel,length(errlevel)*length(arlevel))
   result[,"AR"] <- sort(rep(arlevel,length(errlevel)*length(rholevel)))
   
   for (ar in arlevel) 
   {
     
     b1 <- bpar[1]
     b2 <- bpar[2]
     b3 <- bpar[3]
     b4 <- ar
     
     # Make optionally relative effects sum to one
     
     if (relative == TRUE) {
        
        b1 <- bpar[1]/(sum(bpar) + ar)
        b2 <- bpar[2]/(sum(bpar) + ar)
        b3 <- bpar[3]/(sum(bpar) + ar)
        b4 <- ar/(sum(bpar) + ar)
     }
     for (e in errlevel) 
     { 
       for (r in rholevel)
       {
    
   # variance of y     
   vary <- b1 + b2 + (1+r**2)*(b3) + 2*sqrt(b1)*sqrt(b2)*r + b4 + e       

   # standardized effect sizes
   be1 <- sqrt(b1)/sqrt(vary)
   be2 <- sqrt(b2)/sqrt(vary)
   be3 <- sqrt(b3)*(sqrt(1+r**2)/sqrt(vary))
   be4 <- sqrt(b4)/sqrt(vary)

   # R squared
   rsq <- (vary - e)/vary 
   f2 = rsq / (1 - rsq)
   
   result[((result$error == e) & (result$rho == r) & result$AR == ar),c(4:9)] <- c(be1,be2,be3,be4,rsq, f2)
   
         }    # end rho level
       }      # end error level
     }        # end ar level
  
  return(result)
   
 } # end function


 
 # test

 require(pander)
 
out <- checkPars(rholevel = c(0.0),
                 arlevel = c(0.0),
                 errlevel = c(6.5),
                 bpar = c(.15,.0, .0),
                 relative = FALSE)



pander(out)
pander(out[,-c(1,7)],caption = "Expected parameter values under various conditions")

## expected variance of product xz

1 + cov((x**2),(z**2)) - (cov(x,z)**2)
1 + (2*rho1**2) - (cov(x,z)**2)


pwr.f2.test(u=1, v=c(200-2, 300-2), f2=(.15^2) / (1 - .15^2))


rsq <- f2/(1+f2)
