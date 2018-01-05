# complex_networks_main.R
# purpose: complex network theory reveals the modular structure of diseases and their inter-relationships
# commenced 4/1/18

library(sand)
library(igraph)
library(dplyr)
library(GO.db)
library(GOstats)
library(org.Sc.sgd.db)
library(org.Hs.eg.db)

source("complex_networks_functions.R")  # load in the functions required for this work. 

# HINT (High-quality INTeractomes) is a curated compilation of high-quality protein-protein 
# interactions from 8 interactome resources (BioGRID, MINT, iRefWeb, DIP, IntAct, HPRD, MIPS 
# and the PDB). Contains 12,429 unique proteins with 59,128 interactions between them. http://hint.yulab.org/
ppi_hint <- read.csv("c:\\R-files\\proteins\\HINT-2017.csv", header=TRUE,stringsAsFactors = FALSE)

setwd("C:/R-files/complexnetworks")    # point to where my code lives
#load("complexnets-4thJanuary2018.RData") # load in required data - the contents will change regulary
memory.limit(2010241024*1024) # use more RAM memory (20 GBs)

# intial code from Kolaczyk and Csardi book
setseed(42)
data(ppi.CC)
upgrade_graph(ppi.CC)

summary(ppi.CC)
V(ppi.CC)[ICSC == 1]$color <- "yellow"
V(ppi.CC)[ICSC == 0]$color <- "lightblue"

plot(ppi.CC, vertex.size=5, vertex.label=NA,
     main="Network of interacting proteins responsible for cell communications: \n yellow=involved in ICSC, blue=not involved")

clu<- clusters(ppi.CC)
ppi.CC.gc <- induced.subgraph(ppi.CC,clu$membership==which.max(clu$csize))
nn.ave <- sapply(V(ppi.CC.gc), function(x) mean(V(ppi.CC.gc)[nei(x)]$ICSC))

par(mfrow=c(2,1))
hist(nn.ave[V(ppi.CC.gc)$ICSC ==1], col="yellow", ylim=c(0,30),xlab="Proportion neighbours w/ ICSC",
     main="Egos w/ ICSC")
hist(nn.ave[V(ppi.CC.gc)$ICSC ==0], col="lightblue", ylim=c(0,30),xlab="Proportion neighbours without ICSC",
     main="Egos w/out ICSC")

nn.pred <- as.numeric(nn.ave > 0.5)
mean(as.numeric(nn.pred) != V(ppi.CC.gc)$ICSC)

# now do the GO stuff, as per page 140 of Kolaczyk
x <- as.list(org.Sc.sgdGO2ALLORFS)
current.icst <- x[names(x) == "GO:0035556"]
ev.code <- names(current.icst[[1]])
icst.ida <- current.icst[[1]][ev.code == "IDA"]
orig.icsc <- V(ppi.CC.gc)[ICSC == 1]$name
candidates <- intersect(icst.ida, V(ppi.CC.gc)$name)
new.icsc <- setdiff(candidates, orig.icsc)
new.icsc

