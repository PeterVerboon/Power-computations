
#' Simulation for power moderation function
#'
#' This function allows you to estimate the power in an linear model with interaction effects.
#' @param samSize vector with sample sizes
#' @param alphalevel vector of alpha levels
#' @param n.predictors number of predictors
#' @param l.interactions list of interaction terms, e.g. list(c("x1","x2),c("x1","x3")) 
#' @param cor.pred correlations between predictors
#' @param bpar vector of regression coefficients
#' @param niter number of iterations
#' @keywords power regression interaction moderation
#' @export
#' @import MASS
#' @examples simPower.moderation()
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
    predNames = paste0("x", 1:n.predictors)
    if (!is.null(l.interactions)) {
      for (i in 1:length(l.interactions)) {
        names(l.interactions)[i] <- paste0(l.interactions[[i]][1],l.interactions[[i]][2])
        predNames <- c(predNames, names(l.interactions))
        outNames <- c(predNames, paste0("CI_ low_",predNames),paste0("CI_high_",predNames))
      }
    }

    numrow <- length(alphalevel)*length(samSize)
    numcol= (2 + 3*(npredtot <- n.predictors + length(l.interactions)))
    result <- as.data.frame(matrix(data=0, nrow = numrow, ncol= numcol))
    colnames(result) <- c("N","alpha",outNames )
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

          result[((result$alpha == a) &  result$N == N),c(3:(npredtot+2))] <- out[1,]
          result[((result$alpha == a) &  result$N == N),c((npredtot+3):(2*npredtot+2))] <- out[2,]
          result[((result$alpha == a) &  result$N == N),c((2*npredtot+3):(3*npredtot+2))] <- out[3,]
          print(c(N, a))

         }      # loop over alpha levels
       }        # loop over sample sizes

       colnames(result) <- colnames(result) <- c("N","alpha",outNames )

       output$result <- result

       class(output) <- "simPower.moderation";

       return(output)
} 



#' Plot method for simPower.moderation
#'
#' This function allows you to plot the power of an object of class simPower.moderation
#' @param x an object of class simPower.moderation
#' @param pval type I error
#' @keywords plot powercurve regression interaction moderation
#' @import ggplot2
#' @import viridis
#' @import reshape
#' @export
#' @examples plot(x)
plot.simPower.moderation <- function(x, pval=0.05, ...) {

        dat <- res$result
        datpow <- dat[dat[,"alpha"]== pval,c(1:(length(res$input$bpar)+2))]
        datlow <- dat[dat[,"alpha"]== pval,c(1:2,(length(res$input$bpar)+3):(2*length(res$input$bpar)+2))]
        dathigh <- dat[dat[,"alpha"]== pval,c(1:2,(2*length(res$input$bpar)+3):(3*length(res$input$bpar)+2))]

        mdat <- melt(datpow, id=c("N","alpha"))
        ci_low <- melt(datlow, id=c("N","alpha"))[,-c(1:3)]
        ci_high <- melt(dathigh, id=c("N","alpha"))[,-c(1:3)]
        mdat <- cbind(mdat,ci_low,ci_high)

        levels(mdat$variable) <-  paste0(levels(mdat$variable)," = ",res$input$bpar)
        mdat <- rename(mdat, c(variable = "effectSize"))

p <- ggplot(data=mdat, aes(y=value, x=N,  colour=effectSize, group = effectSize)) + geom_point(size=2) +
        geom_line(size=.5) + labs( y = "Estimated power") +
        geom_errorbar(ymin = ci_low, ymax = ci_high, width=0.2, size=.3) +
        geom_hline(yintercept=0.80, linetype="dashed", color = "red") +
        geom_hline(yintercept=0.90, linetype="dashed", color = "blue") +
        scale_y_continuous(breaks=seq(0.10, 1, 0.10), limits = c(0.1,1)) +
        scale_x_continuous(breaks=seq(50, 500, 50) ) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw(base_size = 12) +
        ggtitle(paste0("alpha = ", pval))

print(p)
}



#' Print method for simPower.moderation
#'
#' This function allows you to plot the power of an object of class simPower.moderation
#' @param x an object of class simPower.moderation
#' @keywords print simPower.moderation 
#' @export
#' @examples print(x)
print.simPower.moderation <- function(x,  ...) {
  cat(" Ran ",x$input$niter , " replications with ", x$input$n.predictors, " predictors,  \n and ",
      length(x$input$l.interactions)," interaction terms. \n\n",
      " You provided these regression coeffcients:\n\n",sep="");
  print(x$input$bpar);
  cat("\n These are the power estimates:\n\n", sep="");
  print(x$result[c(1:(length(res$input$bpar)+2))], digits=3);
}


