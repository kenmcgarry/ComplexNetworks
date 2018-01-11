# build a igraph object and get its statistics
# Do not plot the entire graph as its way too large - plot a small subset only
ghint <- graph.data.frame(ppi_hint)
ghint <- as.undirected(ghint); 
gs    <- get_gstatistics(ghint)
head(gs)

hublist <- find_hubs(gs)
hubtargetlist <- is_hub_target(hublist,drug_targets,ppi_hint)

hlist <- unique(hubtargetlist$Gene)

explore <- induced.subgraph(graph=ghint,vids=unlist(neighborhood(graph=ghint,order=1,nodes=hlist[1])))
length(V(explore))
plot(explore)

# Need to combine "drug_targets" ppi with "ppi_hint"
coreness = graph.coreness(as.undirected(explore))
head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness));

# get subgraphs based on top ranking coreness measure
# create layout
ll <- corenessLayout(explore);
# plot
plot(explore, layout=ll, vertex.size=15, vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main='Coreness');
