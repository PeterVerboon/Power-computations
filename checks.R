
ve <- 1
r <- .0    # 0, 0.1, 0.5



vy <- (.5) + (.3) + (1+r**2)*(.2) + ve

b1 <- sqrt(.5)/sqrt(vy)
b2 <- sqrt(.3)/sqrt(vy)
b3 <- sqrt(.2)/sqrt(vy)

rsq <- ((.5) + (.3) + (1+r**2)*(.2))/vy
