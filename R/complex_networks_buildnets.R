# build a igraph object and get its statistics
# Do not plot the entire graph as its way too large - plot a small subset only
gppi <- graph.data.frame(ppi)
gppi <- as.undirected(gppi) 
gppi <- igraph::simplify(gppi)
gs   <- get_gstatistics(gppi)
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
#plot(explore_subgraph, layout=ll, vertex.size=15, vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main='Coreness')

# create drug to target network
dtn <- drug_targets[,c(1,3)]  # use only drugs and proteins
dtn <- graph.data.frame(dtn)
dtn <- as.undirected(dtn); 
dtn <- igraph::simplify(dtn)  # remove duplicates and self-loops

dtn <- delete.vertices(dtn, V(dtn)[degree(dtn) < 11])
dtn <- delete_isolates(dtn)
dts <- get_gstatistics(dtn)
layout <- layout_nicely(dtn)

# Start building kernal predictors
library(kernlab)

clu <- clusters(ppi_net)
ppi_net.gc <- induced.subgraph(dtn,clu$membership == which.max(clu$csize))
L <- as.matrix(graph.laplacian(ppi_net.gc))
egl <- eigen(L) 
nv <- vcount(ppi_net.gc)

# create and plot the three eigenvectors for the 1st, 2nd and 3rd largest weights
for (j in 1:3){
  evals <- egl$values[1:(nv-j)]
  fevals <- c((evals)^(-j),0)
  #plot(fevals, col="magenta",xlim=c(1,nv),xlab=c("Index i"),ylab=expression(f(gamma[i])))
  evec <- egl$vectors[, (nv-j)]
  vcolors <- character(nv)
  vcolors[evec >= 0] <- "red"
  vcolors[evec <  0] <- "blue"
  vsize <- 30 * sqrt(abs(evec))
  l <- layout.fruchterman.reingold(dtn.gc)
  plot(dtn.gc,layout=l,vertex.color=vcolors,vertex.size=vsize,vertex.label=NA)
}

# set j<-1 and rebuild data, else kernal operions wont work.
j<-1;evals <- egl$values[1:(nv-j)];fevals <- c((evals)^(-j),0);evec <- egl$vectors[, (nv-j)]
Kernal1 <- egl$vectors %*% diag(fevals) %*% t(egl$vectors)
Kernal1 <- as.kernelMatrix(Kernal1)




#V(dtn)[V(graph)$type == 1]$shape <- "square"
#V(dtn)[V(graph)$type == 0]$shape <- "circle"
#plot.igraph(dtn,layout=layout.fruchterman.reingold)
#wc <- cluster_fast_greedy(dtn)#cluster_louvain(dtn)
#modularity(wc)
#mods <- membership(wc)
#plot(wc, dtn)


