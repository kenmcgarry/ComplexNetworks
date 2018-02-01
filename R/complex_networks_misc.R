# complex_networks_misc.R
# part completed functions and code that has yet to be assigned.
# Not called up as yet, as will cause errors.

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
gs <- get_gstatistics(ppi_net)
gs$nodes$central <- abs(gs$nodes$central)



######################### PLOT NETWORKS ######################################
# link salience - might solve hub idenification problem
# https://www.nature.com/articles/ncomms1847
# https://github.com/csgillespie/poweRlaw
# Calculate power law for our protein networks - any diff between hubs and non-hubs?
# also graph coreness https://jcasasr.wordpress.com/2015/02/03/plotting-the-coreness-of-a-network-with-r-and-igraph/
# get subgraphs based on coreness measure

subgraph <- make_ego_graph(ppi_net, order=1, c("RPSA", "PNMA1", "MZT2B","SDC2","TUBGCP4"))

sg <- 4;

coreness = graph.coreness(as.undirected(subgraph[[sg]]))
head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness));
ll <- corenessLayout(subgraph[[sg]]);
plot(subgraph[[sg]], layout=ll, vertex.size=15, vertex.color=colbar[coreness], 
     vertex.label.color= "black", vertex.frame.color=colbar[coreness], main='');


# Try different layouts for best looking graph
net <- subgraph[[sg]]
# create network and plot it
layouts <- grep("^layout_",ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(4,4), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <-do.call(layout,list(net))
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }


