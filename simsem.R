

install.packages("simsem_0.5-14.tar.gz", repos = NULL, type = "source")
install.packages("simsem")
library(simsem)
library(semTools)
library(semPlot)


### model created with buildModel function

popNull <- buildSimModel(EScond = .5,
                        ESmod = .3,
                        ESint = .2,
                        bpath = c(.4,.3,.2,.1),
                        model = 0)

analyzeNull <- buildSimModel(EScond = "a1",
                       ESmod = "a2",
                       ESint = "a3",
                       bpath = c("b1","b2","b3","b4"),
                       model = 0)


cutoff2 <- c(RMSEA = 0.05, CFI = 0.95, TLI = 0.95, SRMR = 0.06)

###  simulate for fixed sample size

simDat <- simulateData(popNull, sample.nobs=100)
fit <- sem(analyzeNull, simDat)

semPaths(fit, layout= "tree2")

summary(fit, standardized=TRUE)
fitmeasures(fit)[c("cfi","rmsea")]

output <- sim(nRep=1000, n=5*nrow(simDat), model = fit, generate = fit)
summary(output)
plotCutoff(output, 0.05)
pValue(fit, output)
cutoff <- getCutoff(output, alpha=0.05)
cutoff2 <- c(RMSEA = 0.05, CFI = 0.95, TLI = 0.95, SRMR = 0.06)

getPowerFit(output, cutoff=cutoff2)

summaryParam(output)


### loop over sampel sizes

output1 <- sim(NULL, n = 100:400, analyzeNull, generate = popNull, lavaanfun = "sem")

summary(output1)
coef(output1)
cutoff <- getCutoff(output1, alpha = 0.05, nVal = 250)
plotCutoff(output1, alpha = 0.05)
plotPowerFit(output1, cutoff=cutoff)
plotPowerFit(output1, cutoff=cutoff2)

summaryParam(output1)


