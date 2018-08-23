simPower.Mediation.growth <-
function(n=200, 
                               EScond = .2, 
                               ESmod = .2,
                               ESint = .2,
                               bpath = c(.4,.3), 
                               rho = c(0,0),
                               error = 1,
                               alpha = 0.05,
                               maxiter = 100) 
{   

 input <- c(EScond,ESmod,ESint,bpath,error,alpha, n, maxiter)
 blabel <- "b1"
 for (i in 2:length(bpath)){
     blabel <- paste0(blabel,",","b",i)
 }
 blabel <- (unlist(strsplit(blabel, ",", fixed = TRUE)))
 
 names(input) <- c("EScondition","ESmoderator","ESinteraction",blabel,"error","alpha","sampleSize", "maxIter")
 
 ES <- c(EScond,ESmod,ESint,bpath,bpath*EScond) 
 
  
# generate lavaan model with specifications
model0 <- buildSimModel(EScond = EScond,
                         ESmod = ESmod,
                         ESint = ESint,
                         bpath = bpath,
                         model = "growth")

# generate lavaan model used in analysis
model <- buildSimModel(EScond = "a1",
                        ESmod = "a2",
                        ESint = "a3",
                        bpath = c("b1","b2"),
                        model = "growth")

res <- matrix(data=0,nrow=maxiter, ncol=18)

colnames(res) <- c("condition","power cond",
                   "moderation","power mod",
                   "interaction","power int",
                   "effect li","power li",
                   "effect ls", "power ls",
                   "indirect i","power ind_li",
                   "indirect s","power ind_ls",
                   "cfi","tli","rmsea","srmr")
  
# Initiate the Progress bar
pb <- txtProgressBar(min = 0, max = maxiter, style = 3)

for (i in 1:maxiter) {  

  
# simulate data according to model specifications
data <- simulateData(model0, sample.nobs = n)

# simulate random data 
   data2 <- simDataSet(n=n, varNames=c('mediator', 'y1', 'y2', 'y3','y4',"condition",'moderator', "interaction"),
                     means = 0, sds = error, 
                     ranges = list(c(condition = c(-3.0, 3.0))),
                     silent=TRUE, seed = NULL, 
                     correlations = c(rho[1], rho[2]))

# add random error
  data3 <- data  + data2
  data3$condition <- as.numeric(cut(data3$condition, breaks = 2))

# analyse data using lavaan
  result <- sem(model, data3)

  res[i,1] <- parameterEstimates(result)[1,5]
  res[i,2] <- parameterEstimates(result)[1,8]
  res[i,3] <- parameterEstimates(result)[2,5]
  res[i,4] <- parameterEstimates(result)[2,8]
  res[i,5] <- parameterEstimates(result)[3,5]
  res[i,6] <- parameterEstimates(result)[3,8]
  res[i,7] <- parameterEstimates(result)[4,5]
  res[i,8] <- parameterEstimates(result)[4,8]
  res[i,9] <- parameterEstimates(result)[5,5]
  res[i,10] <- parameterEstimates(result)[5,8]
  res[i,11] <- filter(parameterEstimates(result), lhs %in% c("ind1"))[,5]
  res[i,12] <- filter(parameterEstimates(result), lhs %in% c("ind1"))[,8]
  res[i,13] <- filter(parameterEstimates(result), lhs %in% c("ind2"))[,5]
  res[i,14] <- filter(parameterEstimates(result), lhs %in% c("ind2"))[,8]
  res[i,c(15:18)] <- fitmeasures(result)[c("cfi","tli","rmsea", "srmr")]
  
  setTxtProgressBar(pb, i)
 }

  sigs <- res[,c(2,4,6,8,10,12,14)] < alpha
  power <- apply(sigs,2,mean)                                # power: count number of significant effects
  bias <- ES - apply(res[,c(1,3,5,7,9,11,13)],2,mean)                   # bias : mean of estimates
  
  output <- list(power = power, bias = bias, raw = res, input = input) 
  
  return(output)
  
  
  }
