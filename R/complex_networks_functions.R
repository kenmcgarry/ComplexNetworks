# complex_networks_functions.R
# helper functions and calls in the libraries
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
library(GSEABase)
library(ROCR)
library(kernlab)
library("e1071")
library(caret)
library(data.table)
library(ranger)
library(randomForest)


data(go)
data(gene_GO_terms)
data(GO_IC)

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

build_network <- function(ppi){
  # unsure if protein targets are part of the giant connected ppi network
  # assign 1=target; 0=non-target to each protein
  
  # remove multiple genes
  ppi$Gene_A <-gsub("\\|.*","",ppi$Gene_A)
  ppi$Gene_B <-gsub("\\|.*","",ppi$Gene_B)  # 
  
  un_targets <- (unique(drug_targets$Gene))        # 1,860 unique protein targets
  length(un_targets)
  un_ppi <- (unique(c(ppi$Gene_A,ppi$Gene_B)))      # 15,792 unique general proteins in ppi
  length(un_ppi)
  
  joint_ppi <- un_targets[un_targets %in% un_ppi]  # 1,293 targets are part of giant connected network (we lose 567 targets!)
  not_ppi <- un_targets[!un_targets %in% un_ppi]  # here are the 567 targets)
  
  # dataframe containing targets and non-target proteins. Annotate with:
  # 1. target; 2. hub; 
  # create ppi network (igraph object) and annotate with target or not target
  ppi_net <- graph.data.frame(ppi)
  ppi_net <- as.undirected(ppi_net); 
  ppi_net <- igraph::simplify(ppi_net)  # remove duplicates and self-loops
  ppi_net <- delete_isolates(ppi_net)
  delete.vertices(igraph::simplify(ppi_net), degree(ppi_net)==0)
  
  V(ppi_net)[1:vcount(ppi_net)]$target <- 0   # Intialise all to zeros
  V(ppi_net)[1:vcount(ppi_net)]$hub <- 0   # Intialise all to zeros
  V(ppi_net)[1:vcount(ppi_net)]$type <- "unknown"   # Intialise protein "type" to unknown
  
  # get main component only - ignore lessor weakly connected groups
  V(ppi_net)$comp <- components(ppi_net)$membership
  ppi_net <- induced_subgraph(ppi_net,V(ppi_net)$comp==1)
  
  # remove from joint_ppi the lost nodes 
  survivors <- V(ppi_net)$name
  joint_ppi <- un_targets[un_targets %in% survivors] 
  
  ppi_net <- set_vertex_attr(ppi_net,"target",joint_ppi,1) # Now assign "1" if protein is a target (very neat!)
  return(ppi_net)
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
  hist(nn.ave[V(ppi.CC.gc)$ICSC ==1], col="yellow", ylim=c(0,30),xlab="Proportion neighbours with ICSC",
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
# returns a list: net and nodes
get_gstatistics <- function(gt) {
  net <- data.frame( 
    modu=igraph::modularity(gt, membership(cluster_walktrap(gt))),
    avepath=igraph::average.path.length(gt),
    nedges=igraph::ecount(gt),
    nverts=igraph::vcount(gt),
    transit=igraph::transitivity(gt),
    diam=igraph::diameter(gt,weights=NA),
    connect=igraph::is.connected(gt))
    
  nodes <- data.frame(   
    closeness=igraph::estimate_closeness(gt,mode="all",cutoff=3),
    degree=(igraph::degree(gt)),
    betweenness=igraph::estimate_betweenness(gt,directed=FALSE,cutoff=3),
    hubness=igraph::hub_score(gt)$vector,
    central=vector(mode="integer", length=net$nverts),
    comm=vector(mode="integer", length=net$nverts))
    
  tmp <- igraph::cluster_walktrap(gt)
  nodes$comm <- as.vector(membership(tmp))
  alpha <- igraph::alpha_centrality(ppi_net,alpha=0.1)  
  nodes$central <- as.vector(alpha)
    
  cat("\nOverall network statistics:")
  cat("\n   Modularity ",net$modu)
  cat("\n   Average path ",net$avepath)
  cat("\n   N edges ",net$nedges)
  cat("\n   N vertices ",net$nverts)
  cat("\n   Transitivity ",net$transit)
  cat("\n   Diameter ",net$diam)
  cat("\n   Is connected? ",net$connect)
  gstats = list(net=net, nodes=nodes)
  return(gstats)
}


# Print out basic network statistics, function is passed the previously calculated info from get_gstatistics().
# prior to use: gs <- get_gstatistics(ppinetwork)
# usage:  display_netstats(gs$net)
display_netstats <- function(net){
  cat("\nOverall network statistics:")
  cat("\n   Modularity ",net$modu)
  cat("\n   Average path ",net$avepath)
  cat("\n   N edges ",net$nedges)
  cat("\n   N vertices ",net$nverts)
  cat("\n   Transitivity ",net$transit)
  cat("\n   Diameter ",net$diam)
  cat("\n   Is connected? ",net$connect)
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
# gs_statistics() to get the required data on "degree" for each protein.
find_hubs <- function(gstats){
  genenames <- as.character(rownames(gstats))
  hublist <- cbind(gstats,genenames)
  cutoff <- quantile(gstats$degree, probs = c(0.70, 0.75, 0.8, 0.85, 0.9, 0.99), na.rm = TRUE) 
  hublist <- filter(hublist,degree > cutoff[1])
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

# goanalysis() will enrich a gene with GO terms
# depends on clusterprofiler library and several other things...
# http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-analysis
go_analysis <- function(yourgenes,ontotype){
  cat("\n",yourgenes)
  eg = bitr(yourgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  ego <- enrichGO(gene          = eg[,2],
                  #universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = ontotype, # one of CC, BP or MF
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  return(ego)
}


# annotate_go() will receive a list of proteins and annotate with GO terms it will return
# a matrix of terms and proteins. Uses the GO data by Daniel Greene
annotate_with_go <- function(ppi_net){
  category <- c("MF","BP","CC")
  #enrich <- c("GO:0017091", "AU-rich element binding","1/1", "23/16982", "0.001354375","RU12","FU")
  nproteins <- length(V(ppi_net)$name)
  allproteins <- V(ppi_net)$name
  
  tempgo <- gene_GO_terms[V(ppi_net)$name]  # Annotate!!
  tempgo <- tempgo[!sapply(tempgo, is.null)]  # Not all proteins have GO annotations so sadly remove them.
  go_proteins <- names(tempgo) # Unfortunately, we are left with only 13,417 proteins.
  
  # sort GO terms by the three categories, breakdown is useful for summary statistics
  cc <- go$id[go$name == "cellular_component"]
  bp <- go$id[go$name == "biological_process"]
  mf <- go$id[go$name == "molecular_function"] 
  temp_cc <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=cc, x))
  temp_bp <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=bp, x))
  temp_mf <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=mf, x))
  
  #terms_by_protein <- split(dm$ID,dm$DiseaseModule)  # do split by disease module
  #terms_by_protein <- unname(terms_by_disease_module)   # Remove names for the moment
  #sim_matrix <- get_sim_grid(ontology=go,information_content=GO_IC,term_sets=terms_by_protein)
  #dist_mat <- max(sim_matrix) - sim_matrix  # need a distance matrix, not a similagrity matrix
 
  return(tempgo)
}


# Add the protein type to the ppi network of target proteins
annotate_with_panther <- function(drugtargets){
  # create drug to target network
  
  # Remove small quantity proteins; Adhesion; Nuclear Other; Antibody; CD Molecules; Ribosomal; Cytokine; Surface Antigen; Membrane other
  drugtargets <-  # Only keep protein target types with at least 50 occurences
      drugtargets %>%
      add_count(TargetClass,sort=TRUE) %>%
      filter(n > 50)
  
  dtn <- drugtargets[,c(1,3)]  # use only drugs and proteins
  dtn[] <- lapply(dtn, as.character) # convert from factors to strings
  dtn <- graph.data.frame(dtn)
  dtn <- as.undirected(dtn); 
  dtn <- igraph::simplify(dtn)  # remove duplicates and self-loops
  
  dtn <- delete.vertices(dtn, V(dtn)[degree(dtn) < 5])
  dtn <- delete_isolates(dtn)
  dts <- get_gstatistics(dtn)
  
  # assign protein types: 
  V(dtn)$GPCR <- 1
  V(target_ppi)$Kinase <- 1
  V(target_ppi)$Enzyme <- 1
  V(target_ppi)$IonChannel <- 1
  V(target_ppi)$Transporter <- 1
  V(target_ppi)$Membrane <- 1
  V(target_ppi)$Structural <- 1
  V(target_ppi)$NucRecept <- 1
  V(target_ppi)$Transcription <- 1
  V(target_ppi)$Secreted <- 1
  V(target_ppi)$Cytosolic <- 1
  V(target_ppi)$Unclassified <- 1
  
  return(anno_target)
}


# go_slim_annotation() reduces the complexities of numerous GO annoatations into a few key terms.
# http://www.geneontology.org/page/go-slim-and-subset-guide#On_the_web. This uses the GSEABase package.
go_slim_annotation <- function(mylist){
  gostuff <- annotate_with_go(ppi_net)
  gostuff <-gostuff[names(which(lapply(gostuff, length) >1))]  # keep only genes with 1 or more GO annotations
  assign("go_mf", TRUE, env=globalenv())  # global variables to enable recovery from errors generated by goSlim()
  assign("go_cc", TRUE, env=globalenv())
  assign("go_bp", TRUE, env=globalenv())
  
  # create the matrix for classification algorithms
  mm <- matrix(0, 149, length(names(gostuff)))  # Number of GO terms in GoSlim (CC+BP+MF) x Number of genes
  colnames(mm) <- names(gostuff)
  rownames(mm) <- give_rownames_mm()
  
  for (i in 1:length(names(gostuff))){
    myCollection <- GOCollection(gostuff[[i]])
    genename <- names(gostuff[i])
    go_empty_mf$Gene <- genename 
    obo <- system.file("extdata","goslim_generic.obo", package="GSEABase") # generic terms by GO consortium
    #obo <- system.file("extdata","goslim_chembl.obo", package="GSEABase") # Chembl Drug Target developed by Mutowo and Lomax
    #obo <- system.file("extdata","goslim_pir.obo", package="GSEABase") #Protein Info Resource by Darren Natale
    slim <- getOBOCollection(obo)
    go_mf <- tryCatch(goSlim(myCollection, slim, "MF"),error=function(e) {go_mf <- error_go_mf()})
    go_cc <- tryCatch(goSlim(myCollection, slim, "CC"),error=function(e) {go_cc <- error_go_cc()})
    go_bp <- tryCatch(goSlim(myCollection, slim, "BP"),error=function(e) {go_bp <- error_go_bp()})
    
    if(length(go_mf) ==1) {mm[1:43,i]  <- as.vector(matrix(0,ncol=43))} else{  # fill with zeros if no annotations found
      go_mf[go_mf$Count != 0,]$Count <- 1; # convert non-zero numbers into 1's
      mm[1:43,i] <- go_mf$Count}            # found MF annotations, assign to matrix
    if(length(go_cc)==1) {mm[44:78,i] <- as.vector(matrix(0,ncol=35))} else{
      go_cc[go_cc$Count != 0,]$Count <- 1
      mm[44:78,i] <- go_cc$Count}
    if(length(go_bp) == 1) {mm[79:149,i]<-  as.vector(matrix(0,ncol=71))} else{
      go_bp[go_bp$Count != 0,]$Count <- 1
      mm[79:149,i] <- go_bp$Count}
  }
  return(mm)  # return matrix  
}



# The following functions are to avoid the errors kicked out by goSlim (and crashing R) 
# when it cant find any terms
error_go_mf <- function(){
  return("MF: Error detected!")
  #go_default <- "error"
  }
error_go_cc <- function(){
  return("CC: Error detected!")
  #go_default <- "error"
  }
error_go_bp <- function(){
  return("BP: Error detected!")
  #go_default <- "error"
  }


# Add the "drug target" class to the "mm" matrix, this will make it ready for classification algorithms.
# but first remove nodes that do have GO terms.
give_classlabels_mm <- function(mm){

  survivors <- colnames(mm)
  mmt <- matrix(0, 1, length(tmp)) 
  rownames(mmt) <- "targets"
  # delete nodes without GO terms from ppi_net 
  biglist <- V(ppi_net)$name
  lostnodes <- survivors[!biglist %in% survivors]
  #lostnodes <- lostnodes[is.na(lostnodes)] <- 0
  ppi_net <-igraph::delete_vertices(ppi_net,lostnodes[1:2020])
  
  targets <- drug_targets$Gene  # get the original drug target data

  for(i in 1:length(survivors)){
    temp <- grep(survivors[i],targets)
    if(length(temp) != 0){mmt[i] <- 1}
  }
  
  mmt <- rbind(mm,mmt)
  return(mmt)
}

give_rownames_mm <- function(){
  gonames <- c("GO:0000988:MF:transcription factor activity, prot...",
    "GO:0003674:MF:molecular_function",
    "GO:0003677:MF:DNA binding",
     "GO:0003700:MF:transcription factor activity, sequ...",
  "GO:0003723:MF:RNA binding",
  "GO:0003729:MF:mRNA binding",
  "GO:0003735:MF:structural constituent of ribosome",
  "GO:0003924:MF:GTPase activity",
  "GO:0004386:MF:helicase activity",
  "GO:0004518:MF:nuclease activity",
  "GO:0004871:MF:signal transducer activity",
  "GO:0005198:MF:structural molecule activity",
  "GO:0008092:MF:cytoskeletal protein binding",
  "GO:0008134:MF:transcription factor binding",
  "GO:0008135:MF:translation factor activity, RNA bi...",
  "GO:0008168:MF:methyltransferase activity",
  "GO:0008233:MF:peptidase activity",
  "GO:0008289:MF:lipid binding",
  "GO:0008565:MF:protein transporter activity",
  "GO:0016301:MF:kinase activity",
  "GO:0016491:MF:oxidoreductase activity",
  "GO:0016746:MF:transferase activity, transferring ...",
  "GO:0016757:MF:transferase activity, transferring ...",
  "GO:0016765:MF:transferase activity, transferring ...",
  "GO:0016779:MF:nucleotidyltransferase activity",
  "GO:0016791:MF:phosphatase activity",
  "GO:0016798:MF:hydrolase activity, acting on glyco...",
  "GO:0016810:MF:hydrolase activity, acting on carbo...",
  "GO:0016829:MF:lyase activity",
  "GO:0016853:MF:isomerase activity",
  "GO:0016874:MF:ligase activity",
  "GO:0016887:MF:ATPase activity",
  "GO:0019843:MF:rRNA binding",
  "GO:0019899:MF:enzyme binding",
  "GO:0022857:MF:transmembrane transporter activity",
  "GO:0030234:MF:enzyme regulator activity",
  "GO:0030533:MF:triplet codon-amino acid adaptor ac...",
  "GO:0030555:MF:RNA modification guide activity",
  "GO:0030674:MF:protein binding, bridging",
  "GO:0032182:MF:ubiquitin-like protein binding",
  "GO:0042393:MF:histone binding",
  "GO:0043167:MF:ion binding",
  "GO:0051082:MF:unfolded protein binding",
  "GO:0000228:CC:nuclear chromosome",
  "GO:0000229:CC:cytoplasmic chromosome",
  "GO:0005575:CC:cellular_component",
  "GO:0005576:CC:extracellular region",
  "GO:0005578:CC:proteinaceous extracellular matrix",
  "GO:0005615:CC:extracellular space",
  "GO:0005618:CC:cell wall",
  "GO:0005622:CC:intracellular",
  "GO:0005623:CC:cell",
  "GO:0005634:CC:nucleus",
  "GO:0005635:CC:nuclear envelope",
  "GO:0005654:CC:nucleoplasm",
  "GO:0005694:CC:chromosome",
  "GO:0005730:CC:nucleolus",
  "GO:0005737:CC:cytoplasm",
  "GO:0005739:CC:mitochondrion",
  "GO:0005764:CC:lysosome",
  "GO:0005768:CC:endosome",
  "GO:0005773:CC:vacuole",
  "GO:0005777:CC:peroxisome",
  "GO:0005783:CC:endoplasmic reticulum",
  "GO:0005794:CC:Golgi apparatus2",
  "GO:0005811:CC:lipid particle",
  "GO:0005815:CC:microtubule organizing center",
  "GO:0005829:CC:cytosol",
  "GO:0005840:CC:ribosome",
  "GO:0005856:CC:cytoskeleton",
  "GO:0005886:CC:plasma membrane",
  "GO:0005929:CC:cilium",
  "GO:0009536:CC:plastid",
  "GO:0009579:CC:thylakoid",
  "GO:0030312:CC:external encapsulating structure",
  "GO:0031410:CC:cytoplasmic vesicle",
  "GO:0043226:CC:organelle",
  "GO:0043234:CC:protein complex",
  "GO:0000003:BP:reproduction",
  "GO:0000278:BP:mitotic cell cycle",
  "GO:0000902:BP:cell morphogenesis",
  "GO:0002376:BP:immune system process",
  "GO:0003013:BP:circulatory system process",
  "GO:0005975:BP:carbohydrate metabolic process",
  "GO:0006091:BP:generation of precursor metabolites...",
  "GO:0006259:BP:DNA metabolic process",
  "GO:0006397:BP:mRNA processing",
  "GO:0006399:BP:tRNA metabolic process",
  "GO:0006412:BP:translation",
  "GO:0006457:BP:protein folding",
  "GO:0006461:BP:protein complex assembly",
  "GO:0006464:BP:cellular protein modification proce...",
  "GO:0006520:BP:cellular amino acid metabolic proce...",
  "GO:0006605:BP:protein targeting",
  "GO:0006629:BP:lipid metabolic process",
  "GO:0006790:BP:sulfur compound metabolic process",
  "GO:0006810:BP:transport",
  "GO:0006913:BP:nucleocytoplasmic transport",
  "GO:0006914:BP:autophagy",
  "GO:0006950:BP:response to stress",
  "GO:0007005:BP:mitochondrion organization",
  "GO:0007009:BP:plasma membrane organization",
  "GO:0007010:BP:cytoskeleton organization",
  "GO:0007034:BP:vacuolar transport",
  "GO:0007049:BP:cell cycle",
  "GO:0007059:BP:chromosome segregation",
  "GO:0007155:BP:cell adhesion",
  "GO:0007165:BP:signal transduction",
  "GO:0007267:BP:cell-cell signaling",
  "GO:0007568:BP:aging",
  "GO:0008150:BP:biological_process",
  "GO:0008219:BP:cell death",
  "GO:0008283:BP:cell proliferation",
  "GO:0009056:BP:catabolic process",
  "GO:0009058:BP:biosynthetic process",
  "GO:0009790:BP:embryo development",
  "GO:0015979:BP:photosynthesis",
  "GO:0016192:BP:vesicle-mediated transport",
  "GO:0019748:BP:secondary metabolic process",
  "GO:0021700:BP:developmental maturation",
  "GO:0022607:BP:cellular component assembly",
  "GO:0022618:BP:ribonucleoprotein complex assembly",
  "GO:0030154:BP:cell differentiation",
  "GO:0030198:BP:extracellular matrix organization",
  "GO:0030705:BP:cytoskeleton-dependent intracellula...",
  "GO:0032196:BP:transposition",
  "GO:0034330:BP:cell junction organization",
  "GO:0034641:BP:cellular nitrogen compound metaboli...",
  "GO:0034655:BP:nucleobase-containing compound cata...",
  "GO:0040007:BP:growth",
  "GO:0040011:BP:locomotion",
  "GO:0042254:BP:ribosome biogenesis",
  "GO:0042592:BP:homeostatic process",
  "GO:0043473:BP:pigmentation",
  "GO:0044281:BP:small molecule metabolic process",
  "GO:0044403:BP:symbiosis, encompassing mutualism t...",
  "GO:0048646:BP:anatomical structure formation invo...",
  "GO:0048856:BP:anatomical structure development",
  "GO:0048870:BP:cell motility",
  "GO:0050877:BP:neurological system process",
  "GO:0051186:BP:cofactor metabolic process",
  "GO:0051276:BP:chromosome organization",
  "GO:0051301:BP:cell division",
  "GO:0051604:BP:protein maturation",
  "GO:0055085:BP:transmembrane transport",
  "GO:0061024:BP:membrane organization",
  "GO:0065003:BP:macromolecular complex assembly",
  "GO:0071554:BP:cell wall organization or biogenesi...",
  "GO:0071941:BP:nitrogen cycle metabolic process")
  
} 

hsing <- c("GO:0005730","GO:0003723", "GO:0005515", "GO:0006412",   # The GO terms used by (Hsing, 2008)
           "GO:0006139","GO:0006996", "GO:0030246", "GO:0005840",
           "GO:0005777","GO:0009719", "GO:0007049", "GO:0004871",
           "GO:0005654","GO:0008219", "GO:0006118", "GO:0006259",
           "GO:0050789","GO:0006950", "GO:0005811", "GO:0008135")

factor2int <- function(targettype){
  targettype <- as.numeric(targettype)
  for (i in 1:length(targettype)){
    if(targettype[i] == 1) {  
      targettype[i] <- 0}}
  for (i in 1:length(targettype)){
    if(targettype[i] == 2) {  
      targettype[i] <- 1}}

  for (i in 1:length(targettype)){
    if(targettype[i] == 1) {  
      targettype[i] <- runif(1, 7.0, 9.9)}}
  for (i in 1:length(targettype)){
    if(targettype[i] == 0) {  
      targettype[i] <- runif(1,0.01,0.02)}}
  return(targettype)
}
