

install.packages("simsem_0.5-14.tar.gz", repos = NULL, type = "source")
install.packages("simsem")
library(simsem)
library(semTools)
library(semPlot)

model1 <- "
Mediator ~ a*IV
DV ~ b*Mediator + c*IV

indirect := a*b
direct   := c
total    := c + (a*b)
"

model1.population <- "
Mediator ~ .5*IV
DV ~ .6*Mediator + .2*IV

indirect := .5*.6
direct   := .2
total    := .2 + (.5*.6)
"


simDat <- simulateData(model1.population, sample.nobs=100)
fit <- sem(model1, simDat)
summary(fit)

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

getPowerFit(output1, cutoff=cutoff2, nVal = 500)

summaryParam(output)


### loop over sampel sizes

output1 <- sim(NULL, n = 100:400, model1, generate = model1.population, lavaanfun = "sem")

summary(output1)
a <- coef(output1)
apply(a, 2, mean)

cutoff <- getCutoff(output1, alpha = 0.05, nVal = 250)
plotCutoff(output1, alpha = 0.05)
plotPowerFit(output1, cutoff=cutoff)
plotPowerFit(output1, cutoff=cutoff2)
getCIwidth(output1, nVal = 500, assurance = .9)
plotCIwidth(output1, "a")
summaryParam(output1)

summaryPopulation(output1)


paramValue(output1)
