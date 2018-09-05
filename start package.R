
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)


setwd("~/Documents/Open Universiteit/Onderzoek/Methodologie/Power-computations")

#create("simPower2")

setwd("./simPower2")
document()

setwd("..")
install("simPower2")



library(simPower2)
