# complex_networks_main.R
# to be submitted to: Journal of Computers in Biology and Medicine 
# purpose: complex network theory reveals the modular structure of diseases and their inter-relationships
# commenced 4/1/18

memory.limit(2210241024*1024) # use more RAM memory (22 GBs)
setwd("C:/R-files/complexnetworks")    # point to where my code lives
source("complex_networks_functions.R")  # load in the functions required for this work. 
source("complex_networks_plots.R")  # load in the plotting functions used in paper. 

#load("complexnets-11thJanuary2018.RData") # load in required data - the contents will change regulary
# example_from_kolaczyk()

# HINT (High-quality INTeractomes) is a curated compilation of high-quality protein-protein 
# interactions from 8 interactome resources (BioGRID, MINT, iRefWeb, DIP, IntAct, HPRD, MIPS 
# and the PDB). Contains 12,429 unique proteins with 59,128 interactions between them. http://hint.yulab.org/
ppi_hint <- read.csv("c:\\R-files\\proteins\\HINT-2017.csv", header=TRUE,stringsAsFactors = FALSE)
drug_targets <- load_drugtargets()  # loads in and preprocesses the file.

bar_plot_drugtargets()  # barplot of the types of drug targets

# Covert gene names all to uppercase
drug_targets$Gene <- toupper(drug_targets$Gene)
ppi_hint <- mutate_all(ppi_hint, funs(toupper))

# get rid of duplicate A and B columns - NB for some reason duplicated() function cannot find them!!!
ppi_hint <- ppi_hint %>% 
    dplyr::filter(Gene_A != Gene_B)

source("complex_networks_buildnets.R")  # create a PPI network, detemine hubs and targets








