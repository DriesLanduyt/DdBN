# ------------
# Author:       Dries Landuyt
# Email:        drieslanduyt@gmail.com
# Description:  Code to develop data-driven BBNs and to assess parameter settings
# ------------

#load packages
library(bnlearn)
library(infotheo)
library(caret)
library(utils)

#learn BBN models using a set of possible parameter settings
runsettings <- function (dataset, contin_cols, disc_cols, targetvariable, numberofdraws, undersampling, save_net, num_states, num_restarts, disc_type, score_type ,BBN_type, kfold, selectednodes, arcs){
  
  #---
  # Function loops over different parameter settings and returns a data frame with in each row a particular parameter setting with its corresponding model performance (% Error,CCI and Kappa).  
  #
  # Main attributes:
  #   - dataset: data frame 
  #   - contin_cols: a vector with the names of the continuous variables in the original data frame
  #   - targetvariable: name of the target variable in the original data frame
  #   - numberofdraws: number of replicates for each parameter setting
  #   - undersampling: Whether the data has to be undersampled or not (TRUE or FALSE) before model calibration
  #   - save_net: Whether or not the network has to be saved (TRUE or FALSE)
  #   - kfold: number of folds for cross-validation
  #   - num_restarts: number of random restarts of the algorithm
  #
  # Settings (should be provided as vectors, more options are possible for each setting):
  #   - num_states: number of intervals each continuous variable is discretized in
  #   - dis_type: type of discretization applied (equal frequency discretization (equalfreq) or equal width discretization (equalwidth))
  #   - score_type: maximization score, used by the algorithm (Akaike's Information Criterion (aic) or Bayesian Information Criterion (bic))
  #   - BBN_type: A free structure (BBN) or a basic classification structure (Naive Bayes classifier (NB) or Tree Augmented Naive Bayes (TAN))
  #---
  
  #define parametersettings data frame
  parametersettings <- expand.grid(list(states = num_states,disc = disc_type,scoretype= score_type,BBNtype=BBN_type))
  parametersettings$scoretype <- as.character(parametersettings$scoretype)
  parametersettings$BBNtype <- as.character(parametersettings$BBNtype)
  parametersettings$disc <- as.character(parametersettings$disc)
  
  #initiate results for return 
  meanresult <- parametersettings[rep(1:nrow(parametersettings),numberofdraws),] #For each parameter setting, one model result is listed (mean of the three cross-validations)
  runningmeanCCI <- rep(NA,numberofdraws)
  runningstdCCI  <- rep(NA,numberofdraws)
  
  #preprocess dataset
  data <- dataprep(dataset,targetvariable, contin_cols, disc_cols) 
  
  #initiate progressbar
  pb <- txtProgressBar(0,numberofdraws,style = 3)
  
  #start simulation
  for (i in c(0:(numberofdraws-1))){
    
    #subsample dataset
    dataset <- subsample(data,kfold,undersampling,targetvariable) 
    npar = nrow(parametersettings)
    
    #loop over parameter settings
    for (r in c(0:(nrow(parametersettings)-1))){
    
      batcherror = rep(NA,kfold)
      batchkappa = rep(NA,kfold)
      batchscore = rep(NA,kfold)
      
      for (j in c(1:kfold)){
        
        # Split data in test and training set
        test <- subset(dataset,dataset$T==j)
        ndata <- nrow(test)
        test$T <- 'test'
        learn <- subset(dataset,dataset$T!=j)
        learn$T <- 'train'
        dataset_tot = rbind(learn,test)
      
        # Learn BBN model model 
        output <- learnBBN(dataset_tot,contin_cols, disc_cols, targetvariable, states = parametersettings[r+1,1],disc = parametersettings[r+1,2],scoretype = parametersettings[r+1,3],BBNtype = parametersettings[r+1,4],num_restarts = num_restarts,save = save_net,selectednodes=selectednodes,arcs=arcs)
        
        batcherror[j] = (output$FP+output$FN)/(output$FP+output$FN+output$TP+output$TN)
        batchkappa[j] = ((output$TP+output$TN)/ndata-((output$TP+output$FP)/ndata*(output$TP+output$FN)/ndata+(output$TN+output$FN)/ndata*(output$TN+output$FP)/ndata))/(1-((output$TP+output$FP)/ndata*(output$TP+output$FN)/ndata+(output$TN+output$FN)/ndata*(output$TN+output$FP)/ndata))
        batchscore[j] = output$score
      }
      
      # Save mean results for each cross validation
      meanresult[i*npar+r+1,'Error'] <- mean(batcherror)
      meanresult[i*npar+r+1,'Kappa'] <- mean(batchkappa)
      meanresult[i*npar+r+1,'CCI'] <- 1 - mean(batcherror)
      meanresult[i*npar+r+1,'Score'] <- mean(batchscore)
      
    }
    
    #calculate moving average after each simulation
    runningmeanCCI[i+1] <- mean(na.omit(meanresult$CCI))
    runningstdCCI[i+1]  <- sd(na.omit(meanresult$CCI))
    
    #print progress
    setTxtProgressBar(pb,i+1)
      
  }
  PlotResults(meanresult,runningmeanCCI, runningstdCCI)
  return (meanresult)

}

PlotResults <- function(result,runningmean, runningstd){
  
  #---
  # Plot running average and independent effects of parametersettings on model performance 
  #---
  
  par(oma=c(0,0,3,1),mfrow=c(1,2))
  plot(runningmean,type = 'l', xlab = 'Iterations',ylab = 'Mean (CCI)')
  plot(runningstd,type = 'l', xlab = 'Iterations',ylab = 'SD (CCI)')
  mtext("Stability plots",side=3, line=1, outer=TRUE, cex=1.5, font=2)
  
  par(oma=c(0,0,3,1),mfrow=c(2,2))
  boxplot(result$Kappa~result$BBNtype,main = 'Network type',ylab="Cohen's Kappa")
  boxplot(result$Kappa~result$disc,main = 'Discretisation method',ylab="Cohen's Kappa")
  boxplot(result$Kappa~result$states,main = 'Number of states',ylab="Cohen's Kappa")
  boxplot(result$Kappa~result$scoretype,main = 'Score type',ylab="Cohen's Kappa")
  mtext("Cohen's kappa",side=3, line=1, outer=TRUE, cex=1.5, font=2)
  
  par(oma=c(0,0,3,1),mfrow=c(2,2))
  boxplot(result$CCI~result$BBNtype,main = 'Network type', ylab = 'CCI')
  boxplot(result$CCI~result$disc,main = 'Discretisation method', ylab = 'CCI')
  boxplot(result$CCI~result$states,main = 'Number of states', ylab = 'CCI')
  boxplot(result$CCI~result$scoretype,main = 'Score type', ylab = 'CCI') 
  mtext("CCI",side=3, line=1, outer=TRUE, cex=1.5, font=2)
  
  par(oma=c(0,0,3,1),mfrow=c(2,2))
  boxplot(result$Score~result$BBNtype,main = 'Network type', ylab = 'Score')
  boxplot(result$Score~result$disc,main = 'Discretisation method', ylab = 'Score')
  boxplot(result$Score~result$states,main = 'Number of states', ylab = 'Score')
  boxplot(result$Score~result$scoretype,main = 'Score type', ylab = 'Score') 
  mtext("Score",side=3, line=1, outer=TRUE, cex=1.5, font=2)
}
  
dataprep <- function(dataset,targetvariable,contin_cols,disc_cols){
  
  #----
  #Function selects all relevant variables and deletes incomplete rows 
  #---
  
  #select all relevant variables
  columns = c(contin_cols,disc_cols,targetvariable)
  dataset <- dataset[,columns]
  
  #delete incomplete rows
  dataset <- na.omit(dataset)
  dataset <- as.data.frame(dataset)

  #print number of complete records
  print ('number of complete records:')
  print (nrow(dataset))
  
  #print number of complete records in case stratified
  print ('number of records in the stratified dataset:')
  print (2*length(which(dataset[,targetvariable]==1)))
  
  return (dataset)
}


subsample <- function(dataset,kfold, undersampling, targetvariable){
  
  #---
  #This function creates a stratified (in case undersampling == TRUE) dataset and adds an extra column 'T' which can be used to 
  #split the dataset in a testset (one-third of the records) and a trainset (two-thirds of the records).
  #A stratified dataset with an extra column 'T' ('train','test') is returned.
  #---
  
  #stratify dataset
  if (undersampling == TRUE){
    
    #split dataset in presence and absence records
    present <- which(dataset[,targetvariable]==1)
    dataset_present <- dataset[present,]
    absent <- which(dataset[,targetvariable]==0)
    dataset_absent <- dataset[absent,]
    
    # select subset of absence data
    absence_rows <-sample(1:nrow(dataset_absent),nrow(dataset_present))
    
    # use the selected cases to take subsample
    dataset_absent_sample <- dataset_absent[absence_rows,]
    dataset_present_sample <- dataset_present
    
    # merge in one stratified dataset 
    dataset <- rbind(dataset_absent_sample,dataset_present_sample)
  
  }
  
  # reshuffle dataset 
  dataset_subsample <- dataset[sample(nrow(dataset)),]
  
  # split in sets for 3-fold cross validation
  rows <- nrow(dataset_subsample)
  folds <- caret::createFolds(c(1:rows), k = kfold, list = TRUE, returnTrain = FALSE)
  
  #mark each subset using the 'T' column
  for (i in c(1:length(folds))){
    dataset_subsample[as.data.frame(folds[i])[,1],'T'] <- i
  }
  dataset_subsample[,'T'] <- as.factor(dataset_subsample[,'T'])
  
  
  return (dataset_subsample)
}


learnBBN <- function (dataset,contin_cols,disc_cols,targetvariable, states = 5,scoretype = "bic",num_restarts = 100,disc = "equalfreq",BBNtype = "BBN",save = FALSE,plot = FALSE,selectednodes, arcs){
  
  #----
  #This function learns a BBN, based on following function attributes:
  #     - dataset: dataframe used to learn the BBNs
  #     - contin_cols: a vector with the names of the continuous variables
  #     - targetvariable: the column number of the target variable
  #     - states: number of classes used for each continuous variables, default set to 5 (3,4,5,10,20)
  #     - scoretype: model evaluation score used in search algorithm , default set to 'bic'('aic', 'bic')
  #     - num_restarts: number of restarts to avoid local maxima, default set to 100
  #     - disc: type of discretization used, default set to 'eaqualfreq' ('eaqualfreq','equalwidth)
  #     - BBNtype: string denoting BBNtype, default set to 'BBN' ('TAN','NB' or 'BBN')
  #     - save: if TRUE the learned BBN is saved in a Netica-readable .DSC format
  #     - plot: if TRUE the learned BBN is plotted in R
  #
  #----
  
  #discretize continuous data
  dataset[,contin_cols] <- infotheo::discretize(dataset[,contin_cols], disc, nbins=states) 

  #define target variable
  targetvar <- targetvariable
  explanvar <- c(contin_cols,disc_cols)
  
  #convert dataset to usable factor format
  for (cn in names(dataset)){dataset[,cn] <- as.factor(as.character(dataset[,cn]))}
  
  #split dataset in testset and learnset
  sets <- split(dataset,dataset$T)
  dataset <- sets$train[,-length(names(sets$train))]
  testset <- sets$test[,-length(names(sets$test))]
  
  #learn network structure
  if (BBNtype == "EXPERT"){
    net <- empty.graph(selectednodes)
    arcs(net, ignore.cycles = TRUE) = matrix(arcs,ncol = 2, byrow = TRUE,dimnames = list(c(), c("from", "to")))
    if (plot){plot(net)}
  }
  
  if (BBNtype == "BBN"){
    net <- hc(dataset, score = scoretype, restart = num_restarts, max.iter = Inf)
    if (plot){plot(net)}
  }
  
  if (BBNtype == "NB"){
    net <- naive.bayes(targetvar,explanvar,dataset)
    if (plot){plot(net)}
  }
  
  if (BBNtype == "TAN"){
    net <- tree.bayes(dataset,targetvar,explanvar)
    if (plot){plot(net)}
  }
  
  #fit conditional probabilities
  bbn = bn.fit(net,dataset)
  
  #save learned BBN
  if(save){
    filename = paste(states,'states',BBNtype,'.dsc',sep = "")
    write.dsc(filename, bbn)
  }
  
  #calculate confusion matrix
  pred <- predict(object=bbn, node=targetvar, data=testset)
  CM <- table(pred,testset[,targetvariable])
  output <- list(TN = CM[1], FP=CM[2], FN=CM[3], TP=CM[4], score=score(net,dataset,type = scoretype), comp = nparams(bbn),arcs = narcs(bbn))
  
  #output is a list with obtained true negatives ($TN), false positives ($FP), false negatives ($FN), true positives ($TP), structure score ($score) and complexity score as number of parameters ($comp) 
  return (output)
}
