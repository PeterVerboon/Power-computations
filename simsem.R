

install.packages("simsem_0.5-14.tar.gz", repos = NULL, type = "source")
install.packages("simsem")
library(simsem)
library(semTools)
library(semPlot)

loading <- matrix(0, 4, 2)
loading[1:4, 1] <- c(1, 1, 1,1)
loading[1:4, 2] <- c(3,2,1,0)
LY <- bind(loading)

facCov <- matrix(NA, 2, 2)
facCovVal <- diag(c(0.2, 0.3))
facCovVal[lower.tri(facCovVal)] <- c(0.1)
facCovVal[upper.tri(facCovVal)] <- c(0.1)
PS <- binds(facCov, facCovVal)
PS

errorCov <- diag(NA, 5)
errorCovVal <- diag(c(0.5, 0.5, 0.6, 0.7, 0.8))
TE <- binds(errorCov, errorCovVal)
TE

path <- matrix(NA, 4, 2)
m <- c(.4,.3,.2,.5)
BE <- bind(path, m)
BE

test <- model(LY=LY, PS=PS, TE=TE, BE=BE, modelType="SEM")

Output <- sim(100, n=200, model = test,lavaanfun = "sem")

getCutoff(Output, 0.05)
plotCutoff(Output, 0.05)
summary(Output)

res <- sim(100,n=100:150, model = model, lavaanfun = "sem")
getCutoff(res, 0.05)
plotCutoff(res, 0.05)
summary(res)

cutoff2 <- c(RMSEA = 0.05, CFI = 0.95, TLI = 0.95, SRMR = 0.06)
getPowerFit(res, cutoff=cutoff2)
plotPowerFit(res, cutoff=cutoff2)

###



Output.NULL <- sim(NULL, n = 100:300, model0)
Output <- sim(NULL, n = 100:300, model0, generate = model,lavaanfun = "sem")

cutoff <- getCutoff(Output, alpha = 0.05, nVal = 250)
plotCutoff(Output, alpha = 0.05)


semPaths(result, layout= "tree2")
