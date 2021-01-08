#Script that 1) creates artificial data and 2) reads in the artificial data such that they comply with Script_MCsimulation_size_power.R

cancerdatasets = list(Breastdata1,CLLdata1,Colorectaldata1,Endometrialdata1,HLdata1,NHLdata1,Ovariandata1,Thyroiddata1)

####FOR REPRODUCIBILITY: Needs the empirical data (run Script_MakeEmpiricalData.R)
#1. Generate artificial data based on graphical model empirical data
N = 190
cancerlabels = c("Breast","CLL","Colorectal","Endometrial","HL","NHL","Ovarian","Thyroid")
varlabels <-c("FA", "CF", "SL", "EF", "PA", "DY", "AP", "CO", "NV")
for (j in 1:length(cancerlabels)) {
  data <- cancerdatasets[[j]]
  S <- cor(data, use = "pairwise.complete.obs")
  test <- EBICglasso(S, n = nrow(data), gamma = 0.5, returnAllResults = TRUE, threshold = TRUE)
  #select model with lowest BIC
  ind = which(test$ebic==min(test$ebic))
  lambdaval = test$lambda[ind]
  glassoresultData = glasso(S, rho = lambdaval)
  sample <- mvrnorm(n = 190, mu = rep(0,ncol(data)), Sigma = glassoresultData$w)
  colnames(sample) <- varlabels
  write.table(sample, file = paste(cancerlabels[j],"data1.dat",sep = ""))
}
#####


#2. Read artificial data
temp = list.files(pattern="*data1.dat")
cancerdatasets = lapply(temp, read.table)
