# complex_networks_main.R
# to be submitted to: Journal of Computers in Biology and Medicine 
# purpose: complex network theory reveals the modular structure of diseases and their inter-relationships
# commenced 4/1/18

#load("complexnets-12thJanuary2018.RData") # load in required data - the contents will change regulary
memory.limit(2210241024*1024) # use more RAM memory (22 GBs)
setwd("C:/R-files/complexnetworks")    # point to where my code lives
source("complex_networks_functions.R")  # load in the functions required for this work. 
source("complex_networks_data.R")  # load in the ppi data sets, GO, OMIM, drug targets etc.
source("complex_networks_buildnets.R")  # create a PPI network, detemine hubs and targets
source("complex_networks_plots.R")  # load in the plotting functions used in paper. 







