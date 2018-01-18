# complex_networks_functions.R
# helper functions and calls in libraries
# started: 5/1/18
library(sand)
library(igraph)
library(GO.db)
#library(GOstats)
library(org.Sc.sgd.db)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(igraph)
library(NCBI2R)
library(rentrez)
library(ggplot2)
library(biomaRt)
library(scales)
library(grid)
library(RColorBrewer)
library(xtable)
library(poweRlaw)
library(ontologySimilarity)
library(ontologyIndex)
library(infotheo)
library(clusterProfiler)
library(linkcomm)

# Makes first letter of string uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# Drug targets of the various drugs - usefully contains protein type (e.g. GPCR) as well.
# from http://drugcentral.org/download
# DrugCentral is a comprehensive drug information resource for FDA drugs and drugs approved outside USA. The 
# resources can be searched using: drug, target, disease, pharmacologic action, terms. 
load_drugtargets <- function(){
  #drug_targets <- file.path('C://R-files//disease','drug.target.interaction.tsv.gz') %>% read.delim(na.strings='',header =TRUE,stringsAsFactors = FALSE) 
  drug_targets <- read.csv(file="C://R-files//disease//drug.target.interaction.tsv", header=TRUE, sep="\t",stringsAsFactors = FALSE)
  names(drug_targets)[names(drug_targets)=="DRUG_NAME"] <- "DrugName"
  names(drug_targets)[names(drug_targets)=="TARGET_CLASS"] <- "TargetClass"
  names(drug_targets)[names(drug_targets)=="GENE"] <- "Gene"
  drug_targets <- drug_targets[,c(1,4,6)]  # Keep only need three variables
  drug_targets$DrugName <- firstup(drug_targets$DrugName)   # convert first letter to uppercase to match existing data
  drug_targets <- na.omit(drug_targets)# remove NA's
  
  # now unlist special entries, I edited the original file and replaced "|" with "/"
  drug_targets<-
  drug_targets %>% 
    mutate(Gene=strsplit(as.character(Gene), "/")) %>%   # symbols=Gene
    unnest(Gene)
  drug_targets$Gene <- toupper(drug_targets$Gene)  # all to uppercase
  
  # shorten some names, for ease printing etc
  drug_targets$TargetClass <- gsub('Nuclear hormone receptor', 'Nuc receptor', drug_targets$TargetClass)
  drug_targets$TargetClass <- gsub('Transcription factor', 'Transcription', drug_targets$TargetClass)
  drug_targets$TargetClass <- gsub('Membrane receptor', 'Membrane', drug_targets$TargetClass)
  
  return(drug_targets)
}

# R provides a tail and head command to view last six and first six elements, so why not the middle six?
# optional parameter n will print more or less if present.
middle <- function(mydata,n){
  len <- nrow(mydata)
  startpoint <- round(len/2)
  if(missing(n)){
    n <- 5}
  endpoint <- startpoint+n
  mydata[startpoint:endpoint,]
}

# For each protein in the list findout how many publiactions it has been mentioned in
get_pubs <- function(protein_list){
  npubs <- count_articles(protein_list)
  return(npubs)
}


# See how many research articles are written about our proteins. Uses rentrez package.
count_articles <- function (protein_list){
  for (i in 1:length(protein_list)){
    
    print(protein_list[i])
    pname <- paste(protein_list[i],'[GENE]) AND (Homo sapiens[ORGN])',sep="")
    ids<-entrez_search(db="pubmed", term=pname,retmax=40000)
    atemp <- cbind(protein_list[i],length(ids$ids))
    
    if(i!=1){
      articles <- rbind(articles,atemp)} 
    else{
      articles <- atemp}
  }
  return(articles)
}

# Uses defunct NCBI apckage, For each protein get the proteins they interact with.
get_interactions <- function(protein_list){
  
  for (i in 1:length(protein_list)){
    print(protein_list[i])
    pname <- paste(protein_list[i],'[sym]',sep="")
    ids<-GetIDs(pname)
    plist<-GetInteractions(ids)
    plist<-unique(plist[13])
    
    cvd<-rep(protein_list[i],nrow(plist))
    ptemp <- cbind(cvd,plist)
    
    if(i!=1){
      ppi <- rbind(ppi,ptemp)} 
    else{
      ppi <- ptemp}
  }
  return(ppi)
}

# Uses defunct NCBI package, Count how many interactions each protein has.
count_interactions <- function(protein_list) {
  for (i in 1:length(protein_list)){
    print(protein_list[i])
    pname <- paste(protein_list[i],'[sym]',sep="")
    ids<-GetIDs(pname)
    plist<-GetInteractions(ids[1])
    plist<-unique(plist[13])
    ptemp <- cbind(protein_list[i],nrow(plist))
    if(i!=1){
      ppi <- rbind(ppi,ptemp)} 
    else{
      ppi <- ptemp}
  }
  
  return(ppi)
}


# intial code from Kolaczyk and Csardi book, good for intro to GO annoations
example_from_kolaczyk <- function(){
  oldw <- getOption("warn")
  options(warn = -1)
  set.seed(42)
  load("ppi_CC.RData")
  upgrade_graph(ppi.CC)
  #save(ppi.CC, file = "ppi_CC.RData")
  
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
  
  options(warn = oldw)
  
}

# Calculate some statistics about the disease gene network
get_gstatistics <- function(gt) {
  gstats <- data.frame(
    modularity=modularity(gt, membership(cluster_walktrap(gt))),
    avepath=average.path.length(gt),
    nedges=ecount(gt),
    nverts=vcount(gt),
    transit=transitivity(gt),
    degree=(degree(gt)),
    diameter=diameter(gt,weights=NA),
    connect=is.connected(gt),
    closeness=closeness(gt),
    betweenness=betweenness(gt,directed=FALSE),
    density=graph.density(gt),
    hubness=hub_score(gt)$vector,
    authority=authority.score(gt)$vector)
  #power=bonpow(gt))
  return(gstats)
}

# 
# https://jcasasr.wordpress.com/2015/02/03/plotting-the-coreness-of-a-network-with-r-and-igraph/
corenessLayout <- function(g) {
  coreness <- graph.coreness(g);
  xy <- array(NA, dim=c(length(coreness), 2));
  
  shells <- sort(unique(coreness));
  for(shell in shells) {
    v <- 1 - ((shell-1) / max(shells));
    nodes_in_shell <- sum(coreness==shell);
    angles <- seq(0,360,(360/nodes_in_shell));
    angles <- angles[-length(angles)]; # remove last element
    xy[coreness==shell, 1] <- sin(angles) * v;
    xy[coreness==shell, 2] <- cos(angles) * v;
  }
  return(xy);
}


# find_hubs() when presented with graph stats object will search for hubs and return list
# it will also add genenames to hublist. Need to create igraph object first then run get
# gs_statistics() to get the required data on "degree2 for each protein.
find_hubs <- function(gstats){
  genenames <- as.character(rownames(gstats))
  hublist <- cbind(gstats,genenames)
  cutoff <- quantile(gstats$degree, probs = c(0.70, 0.75, 0.8, 0.85, 0.9, 0.99), na.rm = T) 
  hublist <- filter(hublist,degree > cutoff[2])
  hublist <- data.frame(lapply(hublist, as.character), stringsAsFactors=FALSE)
  
  return(hublist)
}

# is_hub_target() receives a list of hubs, the drug_target data and the PPI network to see if
# these hub proteins are also targets.
is_hub_target <- function(hlist,dt,ppi){
  hub_targ_list <- dt[1,] # instantiate before use
  gnames <- hlist$genenames
  totalgenes <- c(ppi[,1],ppi[,2])
  cat("\nWe have ",length(unique(totalgenes))," unique genes in PPI network")
  cat("\nWe have ",length(unique(gnames))," unique HUB genes in PPI network")
  for (i in 1:length(gnames)){
    gene <- gnames[i] # get hub genes individually and see if they in lists of targets
    glist <- filter(dt, Gene == gene)  # This bit is OK
    if(nrow(glist) > 0){
      hub_targ_list <- rbind(hub_targ_list,glist) }
  }
  # Line below removes duplicates that appear in two variables
  hub_targ_list <-hub_targ_list[!(duplicated(hub_targ_list[c("DrugName","Gene")]) | duplicated(hub_targ_list[c("DrugName","Gene")], fromLast = TRUE)), ]
  hub_targ_list <- hub_targ_list[-1,]    # 1st entry is rubbish, so remove it
  cat("\nWe have ",length(unique(hub_targ_list$Gene))," unique genes that are hubs AND targets")
  cat("\nWe have ",length(unique(dt$Gene)) - length(unique(hub_targ_list$Gene)),   " unique genes that are targets but NOT hubs")
  cat("\nWe have ",length(unique(dt$Gene)),   " unique genes that are targets in total")
  
  return(hub_targ_list)
}


# KEGG over-representation test
kegg_analysis <- function(yourgenes){
  #cat("\nyourgenes are: ",yourgenes)
  eg = bitr(yourgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db");
  eg<-eg[,2]
  kk <- enrichKEGG(gene= eg, organism= 'hsa', pvalueCutoff = 0.05)
  
  return(kk)
}


delete_isolates <- function(gt) {
  isol <- V(gt)[degree(gt)==0]
  gt <- delete.vertices(gt, isol)
  return(gt)

}




