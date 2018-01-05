# complex_networks_main.R
# to be submitted to: Journal of Computers in Biology and Medicine 
# purpose: complex network theory reveals the modular structure of diseases and their inter-relationships
# commenced 4/1/18

library(sand)
library(igraph)
library(dplyr)
library(GO.db)
library(GOstats)
library(org.Sc.sgd.db)
library(org.Hs.eg.db)

memory.limit(2210241024*1024) # use more RAM memory (22 GBs)
setwd("C:/R-files/complexnetworks")    # point to where my code lives
source("complex_networks_functions.R")  # load in the functions required for this work. 
#load("complexnets-4thJanuary2018.RData") # load in required data - the contents will change regulary
# example_from_kolaczyk()

# HINT (High-quality INTeractomes) is a curated compilation of high-quality protein-protein 
# interactions from 8 interactome resources (BioGRID, MINT, iRefWeb, DIP, IntAct, HPRD, MIPS 
# and the PDB). Contains 12,429 unique proteins with 59,128 interactions between them. http://hint.yulab.org/
ppi_hint <- read.csv("c:\\R-files\\proteins\\HINT-2017.csv", header=TRUE,stringsAsFactors = FALSE)
drug_targets <- load_drugtargets()


mp <- barplot(table(drug_targets), axes = FALSE, axisnames = FALSE)
text(mp, par("usr")[3], labels = labels, srt = 45, cex=.9)
axis(2)














