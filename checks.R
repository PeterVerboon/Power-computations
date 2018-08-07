

# DEFINE FUNCTION
# Assume x and z ~ N(0,1)
# check expected effect sizes using formulas from power moderation paper

 checkPars  <- function(rholevel = c(0.0, 0.3, 0.8),
                  arlevel = c(0,0.8),
                  errlevel = c(0,1,3,9),
                  bpar = c(.5, .3, .2))
   {
   b1 <- bpar[1]
   b2 <- bpar[2]
   b3 <- bpar[3]
 
   numrow <- length(rholevel)*length(arlevel)*length(errlevel)
   result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= 8))
   colnames(result) <- c("AR","error","rho","b1","b2","b3","b4","rsq")
   result[,"error"] <- rep(sort(rep(errlevel,length(rholevel))),length(arlevel))
   result[,"rho"] <- rep(rholevel,length(errlevel)*length(arlevel))
   result[,"AR"] <- sort(rep(arlevel,length(errlevel)*length(rholevel)))
   

     for (e in errlevel) 
     { 
       for (r in rholevel)
       {
         for (ar in arlevel) 
         {
           
   
   vy0 <- (b1) + (b2) + (1+r**2)*(b3) + 2*sqrt(b1)*sqrt(b2)*r             # with interaction, no AR
   vy2 <- vy0 + (ar)*vy0 + e                               # with interaction, with AR, with error
   
   be1 <- sqrt(b1)/sqrt(vy2)
   be2 <- sqrt(b2)/sqrt(vy2)
   be3 <- sqrt(b3)*(sqrt(1+r**2)/sqrt(vy2))
   be4 <- sqrt(ar)/sqrt(vy2)

   rsq <- (vy2 - e)/vy2   # R squared
   
   result[((result$error == e) & (result$rho == r) & result$AR == ar),c(4:8)] <- c(be1,be2,be3,be4,rsq)
   
         } # end ar level
       } # end rho level
     } # end error level
  
  return(result)
   
 } # end function


 
 # test
 
out <- checkPars(rholevel = c(0.0, 0.3, 0.6),
                 arlevel = c(0.0,0.5),
                 errlevel = c(0,1,3,9),
                 bpar = c(.5, .3, .2))




## expected variance of product xz

1 + cov((x**2),(z**2)) - (cov(x,z)**2)
1 + (2*rho1**2) - (cov(x,z)**2)



