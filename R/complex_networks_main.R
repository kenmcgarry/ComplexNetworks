# complex_networks_main.R
# to be submitted to: Journal of Computers in Biology and Medicine 
# Objectives: are hub proteins predisposed to be potential candidates as drug targets for  
# therapeutic interventions?
# commenced 4/1/18
# updated 12/1/19

memory.limit(1010241024*1024) # use more RAM memory (10 GBs)
setwd("C:/common_laptop/R-files/complexnetworks")    # point to where my code lives
#load("complexnets_16thFebruary2018.RData")
source("complex_networks_functions.R")  # load in the functions required for this work. 

# only run these files if you intend building up the data from scratch
source("complex_networks_data.R")  # load in the ppi data sets, GO, OMIM, drug targets etc.
source("complex_networks_buildnets.R")  # create a PPI network, detemine hubs and targets
source("complex_networks_plots.R")  # load in the plotting functions used in paper. 

source("complex_networks_proteins.R")  # examine k-coreness of ppi network

# some plots for paper
bar_plot_gg2(drug_targets,1,"red")  # plot all target proteins
bar_plot_gg2(hubtargetlist,2,"blue")  # plot targets that are also hubs 
plot_power2(gppi)   # graph of degree power law

# build igraph object, ensuring only great connected unit is returned, prune out bad gene names
ppi_net <- build_network(ppi)
gs <- get_gstatistics(ppi_net)
gs$nodes$central <- abs(gs$nodes$central)

# create the matrix used for training classifiers, annotate with GOSLIM terms and then class labels
source("GO_slim_functions.R")  # load in the functions required for this work. 
obo <- "C://common_laptop//R-files//GOslim//goslim_generic.obo" # generic terms by GO consortium

bigdata <- go_slim_annotation()
mm <- give_classlabels_mm(bigdata)



