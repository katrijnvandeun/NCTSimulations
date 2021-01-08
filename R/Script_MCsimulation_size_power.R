#In this syntax, the actual rejection size and power of the network comparison test is estimated using Monte Carlo simulation.
#A data driven approach is taken to define the population networks. Hence the code is structured with a first part that obtains
#the networks for all possible pairs of cancers, relying on the empirical data. For each pair, a significance test on the difference
#in edge weights is performed and, using generated data with a covariance structure defined by the network analysis of the empirical 
#data, the power (rejection size in case of no significant differences) of the hypothesis test is obtained.
#In case of more than one significant edge, the power of detecting at least one significant difference is calculated as well as
#the power of detecting all significant differences (using compound hypothesis testing procedure).
#Binary, continuous, and a mix of binary and continuous data can be used using the mgm function for estimation of the network.
#However, as the mgm function is slow, the glasso function and generation of data from a multivariate normal distribution are
#recommended.


###FIXED INPUT PART FOR CALLING MGM

estimatorArgs <- list()
estimatorArgs$type <- c(rep("g",9))
estimatorArgs$level <- c(rep(1,9))
estimatorArgs$k <- 2
estimatorArgs$alpha <- 1 #lasso only, ridge tuning par = 0
estimatorArgs$sel <- "EBIC"
allcontinuous <- 1 ####ARE ALL VARIABLES CONTINUOUS? Give value 1 for yes, 0 for no. If yes, NCT used and glasso, not mgm!!! (mgm is too slow)
out = NCT_estimator_MGM(x=Breastdata1,argumentlist=estimatorArgs)#test of NCT_estimator_MGM function


######################EXPLORE DATA################

cancerlabels = c("Breast","CLL","Colorectal","Endometrial","HL","NHL","Ovarian","Thyroid")
SUMMARYTABLE = data.frame(matrix(NA,nrow = 28,ncol = 3))
SUMMARYLIST = list()
counter = 1
for (j in 1:(length(cancerlabels)-1)) {
  for (k in (j+1):length(cancerlabels)) {
    data1 <- cancerdatasets[[j]]
    data2 <- cancerdatasets[[k]]
    if (allcontinuous==1){
      compare <- NCT(data1, data2, gamma = 0.5, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges="all", it = 1000)
    } else {
      compare <- NCT_own(data1, data2, gamma = 0.5, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges="all",estimatorArgs = estimatorArgs,it= 100)
    }
    SUMMARYLIST = c(SUMMARYLIST,compare)
    edgelist <- list()
    k2 = 0
    pvls <- c()
    for (j2 in 1:length(compare$einv.pvals$`p-value`)) {
      if (compare$einv.pvals$`p-value`[j2]<.05) {  #Filter out edges that show a significant difference
        k2 = k2+1
        edgelist[[k2]] <- as.numeric(c(compare$einv.pvals$Var1[j2],compare$einv.pvals$Var2[j2]))
        pvls[k2] <- as.numeric(levels(compare$einv.pvals$`p-value`[j2]))
      }
    }
    SUMMARYTABLE[counter,] <- c(cancerlabels[j],cancerlabels[k],k2)
    counter = counter+1
  }
}
#save the summarytable and the summarylist
write.table(SUMMARYTABLE,"Summarytable_glasso.txt",sep=" ",eol="\n")

####START OF MC PART for one particular pair of cancers
N = 190 #sample size of a single dataset
S = 100 #number of dataset pairs to generate in the MC simulation
counter = 1
powerIUT <- vector(length = 28)
powerUIT <- vector(length = 28)
start_time <- Sys.time() #track computation time
for (j2 in 1:(length(cancerlabels)-1)) {
  for (k2 in (j2+1):length(cancerlabels)) {
    data1 <- cancerdatasets[[j2]]
    data2 <- cancerdatasets[[k2]]
    #data1 <- Breastdata1######
    S1 <- cor(data1, use = "pairwise.complete.obs")
    #data2 <- Ovariandata1#######
    S2 <- cor(data2, use = "pairwise.complete.obs")#
    ###Building of population graphical models from which to sample in the MC simulation
    #1. Data driven approach: obtain edges and their weights from the empirical data
    if (allcontinuous==1){
      test <- EBICglasso(S1, n = nrow(data1), gamma = 0.5, returnAllResults = TRUE, threshold = TRUE)
      #select model with lowest BIC
      ind = which(test$ebic==min(test$ebic))
      lambdaval = test$lambda[ind]
      glassoresultData1 = glasso(S1, rho = lambdaval)
      test <- EBICglasso(S2, n = nrow(data2), gamma = 0.5, returnAllResults = TRUE, threshold = TRUE)
      #select model with lowest BIC
      ind = which(test$ebic==min(test$ebic))
      lambdaval = test$lambda[ind]
      glassoresultData2 = glasso(S2, rho = lambdaval)
      NCT_obs <- NCT(data1, data2, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges="all")
    } else {
      out1 = mgm(data1,type=c(rep("g",10)),level = c(rep(1,10)), lambdaSel = "EBIC",k=2)
      sds1 = rep(1,ncol(data1))  #data are standardized so have variance one
      factorinput1 = out1$interactions$indicator #geeft nodige input voor mgmsampler, namelijk de factors en interactions
      out1$interactions$signs[is.na(out1$interactions$signs)]=0
      interactioninput1 <- list()
      interactioninput1[[1]] <- vector("list",length = nrow(out1$interactions$indicator[[1]]))
      for (j in 1:nrow(out1$interactions$indicator[[1]])) {
        interactioninput1[[1]][[j]] = array(out1$interactions$weightsAgg[[1]][[j]]*out1$interactions$signs[[1]][j], dim=c(1, 1))
      }
      out2 = mgm(data2,type=c(rep("g",9)),level = c(rep(1,9)), lambdaSel = "EBIC",k=2)
      factorinput2 = out2$interactions$indicator #geeft nodige input voor mgmsampler, namelijk de factors en interactions
      out2$interactions$signs[is.na(out2$interactions$signs)]=0
      interactioninput2 <- list()
      interactioninput2[[1]] <- vector("list",length = nrow(out2$interactions$indicator[[1]]))
      for (j in 1:nrow(out2$interactions$indicator[[1]])) {
        interactioninput2[[1]][[j]] = array(out2$interactions$weightsAgg[[1]][[j]]*out2$interactions$signs[[1]][j],dim=c(1, 1))
      }
      sds2 = rep(1,ncol(data2))  ####NOTE: mgm standardizes the data
      NCT_obs <- NCT_own(data1, data2, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges="all",estimatorArgs = estimatorArgs)
    }
    #Note: 3 types of sign. testing: gls structures in the output refer to global strength (ie, sum of abs val of weights),
    # nw is for the large difference in weights, while e refers to differences in individual edge weights with -by default - UNcorrected p-values
    ###PROPOSAL: power analysis for all edges showing a significant difference; 
    #
    #Finding edges with p<.05 such that these can be given in the call to nct_own in the MC simulation
    edgelist <- list()
    k = 0
    for (j in 1:length(NCT_obs$einv.pvals$`p-value`)) {
      if (NCT_obs$einv.pvals$`p-value`[j]<.05) {  #Filter out edges that show a significant difference
        k = k+1
        edgelist[[k]] <- as.numeric(c(NCT_obs$einv.pvals$Var1[j],NCT_obs$einv.pvals$Var2[j]))
      }
    }
    #2. Input MC simulatie: number of simulated data sets and sample size (see lines 60, 61)
    #Account for situations where no significantly different edge weight was found and assess size (control over type I error)
    if (k==0){
      out_MC <- matrix(NA,nrow=S,ncol = 1)
    } else {
      out_MC <- matrix(NA,nrow=S,ncol = k)
    }
    for (s in 1:S) {   #THE MC SIMULATION: obtain p-values for S generated data pairs
      if (allcontinuous!=1){
        sample1 = mgmsampler(factors=factorinput1,interactions = interactioninput1,N=N,type=estimatorArgs$type,
                             level=estimatorArgs$level,sds=sds1,thresholds=rep(0,ncol(data1)))  #met mgm sampler 1data simuleren
        sample2 = mgmsampler(factors=factorinput2,interactions = interactioninput2,N=N,type=estimatorArgs$type,
                             level=estimatorArgs$level,sds=sds2,thresholds=rep(0,ncol(data2)))
        nctsamples = NCT_own(sample1$data, sample2$data, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges=edgelist,estimatorArgs = estimatorArgs)
      } else {
        sample1 <- mvrnorm(n = 190, mu = rep(0,ncol(data1)), Sigma = glassoresultData1$w)
        sample2 <- mvrnorm(n = 190, mu = rep(0,ncol(data2)), Sigma = glassoresultData2$w)
        if (k==0){
          NCTsamples <- NCT(sample1, sample2, binary.data=FALSE, paired=FALSE)
          out_MC[s,] <- NCTsamples$nwinv.pval
        } else {
          NCTsamples <- NCT(sample1, sample2, binary.data=FALSE, paired=FALSE, test.edges=TRUE, edges=edgelist)
          out_MC[s,]=as.numeric(paste(NCTsamples$einv.pvals$`p-value`))
        }
      }
    }
    #
    # MGM: approx. 10 minutes for 1 simulation => 17hrs for 100simulations (= one power analysis!!!)
    # glasso: approx. 45 sec for 10 simulations => for 100 = 9 minutes!!!!
    ###glasso: TOOK TWO HOURS AND HALVE for all pairs of cancers
    #
    #Obtain power from out_MC; either looking for at least one significant edge: at least one pval < .05 in row (=> correction for multiple testing needed) 
    #OR
    #all edges significant (this is all pvals in a row < .05)
    #
    if (k==0){
      powerIUT[counter] = sum(rowSums(out_MC < .05)==1)/S
      powerUIT[counter] = sum(rowSums(out_MC < .05)>0)/S
    } else {
      powerIUT[counter] = sum(rowSums(out_MC < .05)==k)/S
      #power
      powerUIT[counter] = sum(rowSums(out_MC < (.05/k))>0)/S
      #power
    }
    counter = counter+1
  }
}
end_time <- Sys.time()
end_time - start_time

SUMMARYTABLE = cbind(SUMMARYTABLE,powerIUT,powerUIT)
write.table(SUMMARYTABLE,"Summarytable_glasso.txt",sep=" ",eol="\n")
