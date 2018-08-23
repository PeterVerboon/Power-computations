simPower.moderation <-
function(samSize = c(50,100,150,200,250,300), 
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

}
