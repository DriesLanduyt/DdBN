# ------------
# Author:       Dries Landuyt
# Email:        drieslanduyt@gmail.com
# Description:  Template form to run 'Data-driven Bayesian Networks' script (DdBN.R)
# ------------

                                   #------------------------------------#
                                   # TEMPLATE FORM TO RUN DdBN.R SCRIPT #
                                   #------------------------------------#

# Fill in all required settings below and run the entire script by selecting all lines and pressing the run button:
# Several plots will be shown and the performance of all developed models will be saved in the data frame 'output'

#------------------------------#
#      LOAD DdBN SCRIPT        #
#------------------------------#

#Specify directory where the 'Data-driven Bayesian Networks' R script (DdBN.R) is located and where the developed networks will be saved
setwd('<here>') 
source ('DdBN.R')

#------------------------------#
#       SPECIFY DATASET        #
#------------------------------#

#specify working directory where dataset is located
setwd('<here>') 

#Specify dataset file (.csv)
dataset = as.data.frame(read.csv('<here>',header=TRUE)) 

#Specify columnnames of the continuous variables in the dataset, exclude response variable
contin_cols = c('<here>','<here>','<here>','<here>','...') 

#Specify columnnames of the discrete variables in the dataset, exclude response variable
disc_cols = c('<here>','<here>','<here>','<here>','...') 

#Specify the response variable
targetvariable = '<here>' 

#Specify whether undersampling is preferred (TRUE or FALSE)
undersampling = TRUE 

#------------------------------#
#  SPECIFY PARAMETER SETTINGS* #
#------------------------------#

# *either for model development (specify optimal setting) or for assessment of the effect of each setting on model performance (specify all possible parameter settings)

#Specify the number of states, it is possible to specify more than one
numberofstates = c(2)  

#Specify the discretiation type, it is possible to specify more than one ('equalfreq','equalwidth')
discretizationtype = c('equalfreq') 

#Specify the score type, it is possible to specify more than one ('bic','aic')
scoretype = c('aic') 

#Specify the network structure type, it is possible to specify more than one ('BBN','NB','TAN','EXPERT')
BBNtype = c('BBN')

#specify the, based on expert knowledge, selected variables, necessary in case the BBN_type vector includes 'EXPERT' 
selectednodes = c('<here>','<here>','<here>','<here>',...)

#specify the causal relations, defined by experts, necessary in case the BBN_type vector includes 'EXPERT'
arcs = c('<parent causal relation 1>','<child causal relation 1>','<parent causal relation 2>','<child causal relation 2>',...)

#------------------------------#
#    SPECIFY OTHER SETTINGS    #
#------------------------------#

#Specify the number of random restarts of the algorithm, it is possible to specify more than one
restarts = 50

#Specify number of replicates for each parametersetting
numberofdraws = 10 

#Specify whether you want to save the network (.dne file) for visualisation (TRUE or FALSE)
savenetwork = FALSE 

#Specify the number of folds for crossvalidation
kfold = 3 

#------------------------------#
#        RUN ALGORITHM         #
#------------------------------#

output = runsettings(dataset, contin_cols, disc_cols, targetvariable, numberofdraws, undersampling, savenetwork, numberofstates, restarts, discretizationtype, scoretype ,BBNtype, kfold, selectednodes, arcs)
