# complex_networks_proteins.R
# Now analyze main protein classes comprising candidate targets
# Examine drug_targets dataframe

p_GPCR <- filter(drug_targets,TargetClass == "GPCR")
p_enzyme <- filter(drug_targets,TargetClass == "Enzyme")
p_IC <- filter(drug_targets,TargetClass == "Ion channel")
p_transporter <- filter(drug_targets,TargetClass == "Transporter")

# Alternatively, assess the drugs targeting the k-core neighborhood.
# subgraphs of k-core proteins

k_cores <- c("APLP1","PSME2","PIP4K2A","PBK","FPR2","FXYD6","PTPRK","COL4A3BP","RPA2","SLC38A6","SLC10A3","SLC16A14","UGDH","STEAP1","FZD4","SLC16A13","GRID1","GOT1L1","MBTPS2",
             "EPHA10", "AKAP8",  "GPR6",   "ADAMTS2", "CELA3A", "MAT1A","FXYD6","PPIG","FZD5","CRYZL1",
              "ACP5","PTPRK","SLC4A2","FPR2","CD47","BRCA2","SELE","KIR3DS1","PSMA3")
explore_subgraph <- induced.subgraph(graph=ppi_net,vids=unlist(neighborhood(graph=ppi_net,order=1,nodes=k_cores[2])))
length(V(explore_subgraph)) 

coreness <- graph.coreness(as.undirected(explore_subgraph))
head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness))
# get subgraphs based on top ranking coreness measure, # create layout
#ll <- corenessLayout(explore_subgraph)
ll <- tkplot(explore_subgraph)  # use tkplot to drag proteins

l <- tkplot.getcoords(ll) # come back and use this command to save placement
plot(explore_subgraph, layout=l, vertex.size=15, 
     vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main="")


# use i=1 for a good diagram i.e. the "APLP1" protein

# Get protein interaction partners for each k-core, see what drugs target these k-cores
explore_subgraph <- induced.subgraph(graph=ppi_net,vids=unlist(neighborhood(graph=ppi_net,order=1,nodes=k_cores[2])))
length(V(explore_subgraph)) 
partners <- V(explore_subgraph)$name

drugcore <- data.frame(DrugName="",TargetClass="",Gene="")
for (i in 1:length(partners)){
  tempcore <- filter(drug_targets,Gene == partners[i])
  if(nrow(tempcore)>0){
    drugcore <- rbind(drugcore,tempcore)
    print(tempcore)
  }
}
drugcore <- drugcore[-1,]    # 1st entry is rubbish, so remove it

xtable(drugcore)





