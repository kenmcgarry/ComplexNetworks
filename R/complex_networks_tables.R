# complex_networks_tables.R
# plot the tables used in the paper.

# Table for nodes sorted on betweenness etc measures.
gs <- get_gstatistics(ppi_net)
display_netstats(gs$net)


print(xtable(head(gs$nodes[order(gs$nodes$betweenness,decreasing=TRUE),],15),
       display=c("s","g","f","f","g","g","f"), math.style.exponents = TRUE,digits=c(0,0,0,2,2,2,1)))

print(xtable(head(gs$nodes[order(gs$nodes$betweenness,decreasing=FALSE),],15),
             display=c("s","g","f","f","g","g","f"), math.style.exponents = TRUE,digits=c(0,0,0,2,2,2,1)))

print(xtable(head(gs$nodes[order(gs$nodes$hubness,decreasing=TRUE),],15),
             display=c("s","g","f","f","g","g","f"), math.style.exponents = TRUE,digits=c(0,0,0,2,2,2,1)))

