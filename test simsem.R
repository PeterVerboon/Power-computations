
popNull <- "
y2 ~ 0.4*y1
y3 ~ 0.4*y2
y4 ~ 0.4*y3
y5 ~ 0.4*y4
y1 ~~ 1*y1
y2 ~~ 0.64*y2
y3 ~~ 0.64*y3
y4 ~~ 0.64*y4
y5 ~~ 0.64*y5
"

popAlt <- "
y2 ~ 0.4*y1
y3 ~ 0.4*y1
y4 ~ 0.4*y2 + 0.4*y3
y5 ~ 0.4*y4
y1 ~~ 1*y1
y2 ~~ 0.64*y2
y3 ~~ 0.64*y3
y4 ~~ 0.40*y4
y5 ~~ 0.64*y5
"

analyzeNull <- "
y2 ~ y1
y3 ~ y2
y4 ~ y3
y5 ~ y4
"

Output.NULL <- sim(NULL, n = 100:500, analyzeNull, generate = popNull, lavaanfun = "sem")
Output.ALT <- sim(NULL, n = 100:500, analyzeNull, generate = popAlt, lavaanfun = "sem")

cutoff <- getCutoff(Output.NULL, alpha = 0.05, nVal = 250)
plotCutoff(Output.NULL, alpha = 0.05)
getPowerFit(Output.ALT, nullObject = Output.NULL, alpha = 0.05, nVal = 250)
getPowerFit(Output.ALT, cutoff = cutoff, nVal = 250, condCutoff = TRUE)
plotPowerFit(Output.ALT, Output.NULL, alpha = 0.05)

summaryParam(Output.NULL)
summaryParam(Output.ALT)

cutoff2 <- c(RMSEA = 0.05, CFI = 0.95, TLI = 0.95, SRMR = 0.06)
getPowerFit(Output.ALT, cutoff = cutoff2, nVal = 250,  condCutoff = FALSE)
plotPowerFit(Output.ALT, cutoff = cutoff2)
