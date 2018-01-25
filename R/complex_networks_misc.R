# complex_networks_misc.R
# part completed functions and code that has yet to be assigned.
# Not called up as yet, as will cause errors.

# unsure protein targets are part of the giant connected ppi network
# assign 1=target; 0=non-target to each protein

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

V(ppi_net)[1:vcount(ppi_net)]$target <- 0   # Intialise all to zeros
V(ppi_net)[1:vcount(ppi_net)]$hub <- 0   # Intialise all to zeros
V(ppi_net)[1:vcount(ppi_net)]$type <- "unknown"   # Intialise protein "type" to unknown

ppi_net <- set_vertex_attr(ppi_net,"target",joint_ppi,1) # Now assign "1" if protein is a target

#digenes <- file.path('C://R-files//proteins', 'all_gene_disease_associations.tsv.gz') %>% read.delim(na.strings='')

# read in protein classification based on pharos protein families
protein_class <- read.csv("c:\\R-files\\proteins\\pharos_v4.6.2.csv", header=TRUE,stringsAsFactors = FALSE,na.strings=c("", "NA"))
protein_class <- protein_class[,c(3,8)]
names(protein_class)[names(protein_class)=="DTO.Family"] <- "TargetClass"
names(protein_class)[names(protein_class)=="HGNC.Sym"] <- "Gene"
protein_class <- protein_class[!(duplicated(protein_class[c("TargetClass","Gene")]) | duplicated(protein_class[c("TargetClass","Gene")], fromLast = TRUE)), ]
protein_class <- na.omit(protein_class)


#ppi_net <- delete.vertices(ppi_net, V(ppi_net)[degree(ppi_net) < 11])   # probably wrong thing to do!
######################### PLOT NETWORKS ######################################
# link salience - might solve hub idenification problem
# https://www.nature.com/articles/ncomms1847
# https://github.com/csgillespie/poweRlaw
# Calculate power law for our protein networks - any diff between hubs and non-hubs?
# also graph coreness https://jcasasr.wordpress.com/2015/02/03/plotting-the-coreness-of-a-network-with-r-and-igraph/
# get subgraphs based on coreness measure
subgraph <- make_ego_graph(ghint, order=1, c("RPSA", "RPS2", "RPL7"))
plot(subgraph[[3]])

coreness = graph.coreness(as.undirected(subgraph[[1]]))
head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness));

# get subgraphs based on top ranking coreness measure
subgraph <- make_ego_graph(ghint, order=1, c("RPSA", "RPS2", "RPL7","RPL10","TP53"))
plot(subgraph[[5]])

# create layout
ll <- corenessLayout(subgraph[[5]]);
# plot
plot(subgraph[[5]], layout=ll, vertex.size=15, vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main='Coreness');

# Try different layouts for best looking graph
net <- subgraph[[1]]
# create network and plot it
layouts <- grep("^layout_",ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <-do.call(layout,list(net))
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }
######################### PLOT NETWORKS ######################################


