# build a igraph object and get its statistics
# Do not plot the entire graph as its way too large - plot a small subset only
gppi <- graph.data.frame(ppi)
gppi <- as.undirected(gppi); 
gppi <- igraph::simplify(gppi)
gs    <- get_gstatistics(gppi)
head(gs)

# Need an unlist type function to unravel "GABRA4|GABRB2|GABRD" multiple entries in drug_target

hublist <- find_hubs(gs)
hubtargetlist <- is_hub_target(hublist,drug_targets,ppi)

hlist <- unique(hubtargetlist$Gene)

explore_subgraph <- induced.subgraph(graph=gppi,vids=unlist(neighborhood(graph=gppi,order=1,nodes=hlist[10])))
length(V(explore_subgraph))
plot(explore_subgraph)

coreness = graph.coreness(as.undirected(explore_subgraph))
head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness))

# get subgraphs based on top ranking coreness measure
# create layout
ll <- corenessLayout(explore_subgraph)
plot(explore_subgraph, layout=ll, vertex.size=15, vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main='Coreness')

# create drug to target network
dtn <- drug_targets[,c(1,3)]  # use only drugs and proteins
dtn <- graph.data.frame(dtn)
dtn <- as.undirected(dtn); 
dtn <- igraph::simplify(dtn)  # remove duplicates and self-loops

dtn <- delete.vertices(dtn, V(dtn)[degree(dtn) < 11])
dtn <- delete_isolates(dtn)
dts <- get_gstatistics(dtn)
layout <- layout_nicely(dtn)

#V(dtn)[V(graph)$type == 1]$shape <- "square"
#V(dtn)[V(graph)$type == 0]$shape <- "circle"

#plot.igraph(dtn,layout=layout.fruchterman.reingold)
wc <- cluster_fast_greedy(dtn)#cluster_louvain(dtn)
modularity(wc)
mods <- membership(wc)
#plot(wc, dtn)

# Use linkcomm to obtain useful communities
dtn <- drug_targets[,c(1,3)]
com1 <- getLinkCommunities(dtn, hcmethod = "single",plot = FALSE)
mods <- getCommunityConnectedness(com1, conn = "modularity")
#plot(com1, type = "commsumm", summary = "modularity")
#plot(com1, type = "dend")

nf <- graph.feature(com1,type="nodes", indices = getNodesIn(com1,clusterids=100,type="indices"),features=30,default = 5)
plot(com1, type = "graph", vsize = nf, vshape = "circle", vlabel = TRUE)
plot(com1, type = "graph", shownodesin = 50, node.pies = TRUE)

ef <- graph.feature(com1, type = "edges", indices = getEdgesIn(com1, clusterids = 14),features = 5, default = 1)
plot(com1, type="graph", ewidth = ef)

# preprocess dtn for OCG
dtn$DrugName <- substr(dtn$DrugName, 0, 20)
dtn$Gene <- substr(dtn$Gene, 0, 20)
dtn <- data.frame(lapply(dtn, function(x) {gsub(" ", "_", x)}))
oc <- getOCG.clusters(dtn)
plot(oc, type = "graph", shownodesin = 7, scale.vertices = 0.1)


