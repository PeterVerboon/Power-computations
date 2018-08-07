require(simr)

x <- rep(1:20)
g <- c('a', 'b', 'c')

X <- expand.grid(x=x, g=g)

b <- c(2, -0.1) # fixed intercept and slope
V1 <- 1.0 # random intercept variance
V2 <- matrix(c(0.5,0.05,0.05,0.1), 2) # random intercept and slope variance-covariance matrix
s <- 1 # residual variance

model1 <- makeLmer(y ~ x + (1|g), fixef=b, VarCorr=V1, sigma=s, data=X)

####

print(model1)

powerSim(model1, nsim=20)

##############

model1 <- makeLmer(y2 ~  x + z + xz + y2L1 + ( 1 |subjnr), fixef=c(0,betas,ar), VarCorr=V1, sigma=1, data=dat)

model1 <- makeLmer(y2 ~  x + z + xz  + ( 1 |subjnr), fixef=c(0,betas), VarCorr=V1, sigma=1, data=dat)


model2 <- extend(model1, along="subjnr", n=15) 

powerSim(model2)

pc2 <- powerCurve(model2) 
print(pc2)
plot(pc2)

save(pc2, file= "result_MLAmod_alpha05.Rdata")
#load("result_moderation_alpha05.Rdata")
