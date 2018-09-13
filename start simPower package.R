
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)


setwd("~/Documents/Open Universiteit/Onderzoek/Methodologie/Power-computations")

#create("simPower2")

setwd("./simPower2")
document()



### TEST

setwd("..")
install("simPower2")

library(simPower2)

res <- simPwr.growth(n=200,  
                     EScond = .4,
                     ESint = .2,
                     bpath = c(.4,-0.1), 
                     rho = c(0.1,0.1),
                     error = 1,
                     alpha = 0.05,
                     maxiter = 100) 

print(res, var = "effect li", plot = TRUE)
