
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)


setwd("~/Documents/Open Universiteit/Onderzoek/Methodologie/Power-computations")

create("simPower")

setwd("./simPower")
document()

use_package("dplyr") 

### TEST

setwd("..")
install("simPower")

library(simPower)

res <- simPwr.growth(n=200,  
                     EScond = .4,
                     ESint = .2,
                     bpath = c(.4,-0.1), 
                     ndepend = 5,
                     rho = c(0.1,0.1),
                     error = 1,
                     alpha = 0.05,
                     maxiter = 50) 

print(res, var = "effect li", plot = TRUE)
print.simPwr.growth(res, var = "effect li", plot = TRUE)
