# build a igraph object and get its statistics
# Do not plot the entire graph as its way too large - plot a small subset only
gppi <- graph.data.frame(ppi)
gppi <- as.undirected(gppi); 
gppi <- simplify(gppi)
gs    <- get_gstatistics(gppi)
head(gs)

# Need an unlist type function to unravel "GABRA4|GABRB2|GABRD" multiple entries in drug_target

hublist <- find_hubs(gs)
hubtargetlist <- is_hub_target(hublist,drug_targets,ppi)

hlist <- unique(hubtargetlist$Gene)

explore_subgraph <- induced.subgraph(graph=gppi,vids=unlist(neighborhood(graph=gppi,order=1,nodes=hlist[1])))
length(V(explore_subgraph))
plot(explore_subgraph)

coreness = graph.coreness(as.undirected(explore_subgraph))
head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness))

# get subgraphs based on top ranking coreness measure
# create layout
ll <- corenessLayout(explore_subgraph)
plot(explore_subgraph, layout=ll, vertex.size=15, vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main='Coreness')


