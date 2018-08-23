plot.simPower.moderation <-
function(x, ...) {
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
