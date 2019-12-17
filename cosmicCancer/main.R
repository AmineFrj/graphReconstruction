
#----------------------------------- PART 1 -----------------------------------
library(miic)
library(bnlearn)
library(igraph)
library(qgraph)

data("cosmicCancer")

hc(cosmicCancer)
## We can replace the missing value with the mode or the mean for example

## Replace the samples which has at least one NA
cCancer <- cosmicCancer[complete.cases(cosmicCancer), ]
## Convert the ploidy to factor
cCancer$Ploidy <- as.factor(cCancer$Ploidy)
## remove variables which have one level
cCancer <- cCancer[, sapply(cCancer, nlevels) > 1]
## Since we removed some samples, 
  #some samples can have constant number of classes, so we have to remove them 
same <- sapply(cCancer, function(.col){
  all(.col[1L] == .col)
})

cCancer <- cCancer[!same]

## Perform the hc on the dataset
hc = bnlearn::hc(cCancer, score = "bic")

## Get the adj matrix 
adj.hc = amat(hc)

## Create the igraph::graph 
graph <- graph_from_adjacency_matrix(adj.hc,
                                     diag = FALSE,
                                     weight = TRUE,
                                     mode = "directed")
## To extract the color 
nodes <- as.data.frame(as.matrix(V(graph)))
nodes$col <- 3

for (i in seq(1,dim(nodes)[1]-1)) {
  row <- rownames(nodes)[i]
  if (tolower(row) == row)
    nodes[i,'col'] <- 1
  else 
    nodes[i,'col'] <- 2
}

V(graph)$color <- nodes$col

## Identify the isolated nodes to ignore them in the plot
Isolated = which(igraph::degree(graph)==0)
g2 = delete.vertices(graph, Isolated)

e <- get.edgelist(g2, names = F)

l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g2),
                                       area=8*(vcount(g2)^2),repulse.rad=(vcount(g2)^3.1))

plot(g2, layout=l,vertex.size=4,vertex.label=NA, edge.arrow.size=.3)
mtext("HC on Cosmic Cancer", side=1)

legend(x = -1.7, y = 1.3, c("Mutated genes", "Over/Under expressed genes","Ploidy"), 
       pch = 21, col = "dodgerblue4", pt.bg = c("dodgerblue1","darkorange","forestgreen"),
       pt.cex = 1.8, cex = .8, bty = "n")

## Identify variables related 'Ploidy'
for (i in seq(1,dim(adj.hc)[1])) {
  if(adj.hc[i,dim(adj.hc)[2]] == 1)
    print(rownames(adj.hc)[i])
}

## Idendify the mutated genes that are significantly related to gene expression
adj.hc <- adj.hc[,-dim(adj.hc)[1]]
adj.hc <- adj.hc[-dim(adj.hc)[1],]

mutated <- c()
for (i in seq(1,dim(adj.hc)[1])) {
  row <- rownames(adj.hc)[i]
  if (tolower(row) == row)
    mutated <- c(mutated,i)
}

adj.hc <- adj.hc[ ,-mutated]
adj.hc <- adj.hc[mutated,]
## We mutated genes in rows and gene expression in columns
  # So we just have to sum the rows
gene_mut <- as.data.frame(rowSums(adj.hc))
gene_mut$names <- rownames(adj.hc)
gene_mut <- gene_mut[order(gene_mut$`rowSums(adj.hc)`,decreasing = T),]
head(gene_mut) ## ==> TP53

centrality <- cbind(as.data.frame(centr_betw(g2)$res),as.matrix(V(g2)))
colnames(centrality) <- c("centrality","i")
centrality <- centrality[order(centrality$centrality,decreasing = T),]
head(centrality)
#------------------------------------------------------------------------------#
#----------------------------------- PART 2 -----------------------------------#
#------------------------------------------------------------------------------#

library(pcalg)

cCancer <- cosmicCancer[complete.cases(cosmicCancer), ]
cCancer$Ploidy <- as.factor(cCancer$Ploidy)
cCancer <- cCancer[, sapply(cCancer, nlevels) > 1]

same <- sapply(cCancer, function(.col){
  all(.col[1L] == .col)
})

cCancer <- cCancer[!same]

### Call the pc algorithm
dm <- as.matrix(data.matrix(cCancer))-1
View(head(dm))
## labels aka node names
V <- colnames(dm) 
levels <- sapply(cCancer[,sapply(cCancer, is.factor)], nlevels)
View(levels)
suffStat = list(dm = dm, nlev = levels, adaptDF = FALSE)
## !! L'execution prend un peu de temps ... !! ##
pc.D <- pcalg::pc( suffStat, 
                   indepTest = disCItest, ## indep.test: disc
                   alpha=0.01, labels = V, verbose = TRUE)


## Convert the result to BN object
bn <- bnlearn::as.bn(pc.D, check.cycles = F)
#class(bn)
adjacency <- amat(bn)
# Build the Graph
graph <- graph_from_adjacency_matrix(adjacency,
                                     diag = FALSE,
                                     weight = TRUE,
                                     mode = "directed")

## To extract the color 
nodes <- as.data.frame(as.matrix(V(graph)))
nodes$col <- 3

for (i in seq(1,dim(nodes)[1]-1)) {
  row <- rownames(nodes)[i]
  if (tolower(row) == row)
    nodes[i,'col'] <- 1
  else 
    nodes[i,'col'] <- 2
}


V(graph)$color <- nodes$col

## Identify the isolated nodes to ignore them in the plot
Isolated = which(igraph::degree(graph)==0)
g2 = delete.vertices(graph, Isolated)

e <- get.edgelist(g2, names = F)

l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g2),
                                       area=8*(vcount(g2)^2),repulse.rad=(vcount(g2)^3.1))

plot(g2, layout=l,vertex.size=4,vertex.label=NA, edge.arrow.size=.2)
mtext("PC on Cosmic Cancer", side=1)

legend(x = -1.8, y = 1.1, c("Mutated genes", "Over/Under expressed genes","Ploidy"), 
       pch = 21, col = "dodgerblue4", pt.bg = c("dodgerblue1","darkorange","forestgreen"),
       pt.cex = 1.8, cex = .8, bty = "n")


############################## significance level :: 0.05
######################################################################
pc.D_0.5 <- pcalg::pc( suffStat, 
                   indepTest = disCItest, ## indep.test: disc
                   alpha=0.035, labels = V, verbose = TRUE)

## Convert the result to BN object
bn <- bnlearn::as.bn(pc.D_0.5, check.cycles = F)
#class(bn)
adjacency <- amat(bn)
# Build the Graph
graph <- graph_from_adjacency_matrix(adjacency,
                                     diag = FALSE,
                                     weight = TRUE,
                                     mode = "directed")

## To extract the color 
nodes <- as.data.frame(as.matrix(V(graph)))
nodes$col <- 3

for (i in seq(1,dim(nodes)[1]-1)) {
  row <- rownames(nodes)[i]
  if (tolower(row) == row)
    nodes[i,'col'] <- 1
  else 
    nodes[i,'col'] <- 2
}


V(graph)$color <- nodes$col

## Identify the isolated nodes to ignore them in the plot
Isolated = which(igraph::degree(graph)==0)
g2_0.5 = delete.vertices(graph, Isolated)

e <- get.edgelist(g2_0.5, names = F)

l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g2_0.5),
                                       area=8*(vcount(g2_0.5)^2),repulse.rad=(vcount(g2_0.5)^3.1))

par(mfrow=c(1,2))
plot(g2, layout=l,vertex.size=4,vertex.label=NA, edge.arrow.size=.2)
mtext("PC wich alpha = 0.01", side=1)
plot(g2_0.5, layout=l,vertex.size=4,vertex.label=NA, edge.arrow.size=.2)
mtext("PC wich alpha = 0.05", side=1)

legend(x = -1.8, y = 1.1, c("Mutated genes", "Over/Under expressed genes","Ploidy"), 
       pch = 21, col = "dodgerblue4", pt.bg = c("dodgerblue1","darkorange","forestgreen"),
       pt.cex = 1.8, cex = .8, bty = "n")
###################################################################################
###################################################################################


## Identify variables related 'Ploidy'
for (i in seq(1,dim(adjacency)[1])) {
  if(adjacency[i,dim(adjacency)[2]] == 1)
    print(rownames(adjacency)[i])
}

## Idendify the mutated genes that are significantly related to gene expression
adjacency <- adjacency[,-dim(adjacency)[1]]
adjacency <- adjacency[-dim(adjacency)[1],]

mutated <- c()
for (i in seq(1,dim(adjacency)[1])) {
  row <- rownames(adjacency)[i]
  if (tolower(row) == row)
    mutated <- c(mutated,i)
}

adjacency <- adjacency[ ,-mutated]
adjacency <- adjacency[mutated,]
## We mutated genes in rows and gene expression in columns
# So we just have to sum the rows
gene_mut <- as.data.frame(rowSums(adjacency))
gene_mut$names <- rownames(adjacency)
gene_mut <- gene_mut[order(gene_mut$`rowSums(adjacency)`,decreasing = T),]
View(head(gene_mut)) ## ==> TP53


centrality <- cbind(as.data.frame(centr_betw(g2)$res),as.matrix(V(g2)))
colnames(centrality) <- c("centrality","i")
centrality <- centrality[order(centrality$centrality,decreasing = T),]
head(centrality)

#------------------------------------------------------------------------------€
#----------------------------------- PART 3 -----------------------------------€
#------------------------------------------------------------------------------€
?miic

# EXAMPLE CANCER
data(cosmicCancer)
data(cosmicCancer_stateOrder)

##  confidenceShuffle : specifying whether edge-specific confidence ratios should be evaluated and used to filter the least robust edges
##  confidenceThreshold : Threshold above which corresponding edges are discarded from the network skeleton.


# execute MIIC (reconstruct graph) 
miic.res = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
                confidenceShuffle = T, confidenceThreshold = 0.001)
# plot graph
miic.plot(miic.res, igraphLayout=igraph::layout_on_grid)

# try with different shuffle and threshold confidence parmeters
miic.res.2 = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
                confidenceShuffle = F, confidenceThreshold = 0.001)
miic.res.3 = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
                confidenceShuffle = T, confidenceThreshold = 0.01)
miic.res.4 = miic(inputData = cosmicCancer, categoryOrder = cosmicCancer_stateOrder, latent = TRUE,
                confidenceShuffle = T, confidenceThreshold = 0.0001)

#Adj matrix
adj.miic = as(miic.res, "matrix")
adj.miic.2 = as(miic.res.2, "matrix")
adj.miic.3 = as(miic.res.3, "matrix")
adj.miic.4 = as(miic.res.4, "matrix")

graph <- graph_from_adjacency_matrix(miic.res$adjMatrix, diag = FALSE, weight = TRUE, mode = "undirected")
graph2 <- graph_from_adjacency_matrix(miic.res.2$adjMatrix, diag = FALSE, weight = TRUE, mode = "undirected")
graph3 <- graph_from_adjacency_matrix(miic.res.3$adjMatrix, diag = FALSE, weight = TRUE, mode = "undirected")
graph4 <- graph_from_adjacency_matrix(miic.res.4$adjMatrix, diag = FALSE, weight = TRUE, mode = "undirected")

## To extract the color 
nodes <- as.data.frame(as.matrix(V(graph4)))
nodes$col <- 3

for (i in seq(1,dim(nodes)[1]-1)) {
  row <- rownames(nodes)[i]
  if (tolower(row) == row)
    nodes[i,'col'] <- 1
  else 
    nodes[i,'col'] <- 2
}


V(graph)$color <- nodes$col
V(graph2)$color <- nodes$col
V(graph3)$color <- nodes$col
V(graph4)$color <- nodes$col

## Identify the isolated nodes to ignore them in the plot
show_graph <- function(graph){
  Isolated = which(degree(graph)==0)
  g2 = delete.vertices(graph, Isolated)
  
  e <- get.edgelist(g2, names = F)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g2),
                                         area=8*(vcount(g2)^2),repulse.rad=(vcount(g2)^3.1))
  
  plot(g2, layout=l, vertex.size=4, vertex.label=NA, edge.arrow.size=.2)
}

par(mfrow=c(2,2))
show_graph(graph); mtext("confidenceShuffle: T, confidenceThreshold: 0.001", side=1)
legend(x = -1.8, y = 1.1, c("Mutated genes", "Over/Under expressed genes","Ploidy"), 
       pch = 21, col = "dodgerblue4", pt.bg = c("dodgerblue1","darkorange","forestgreen"),
       pt.cex = 1.5, cex = .6, bty = "n")
show_graph(graph2); mtext("confidenceShuffle: F, confidenceThreshold: 0.001", side=1)
show_graph(graph3); mtext("confidenceShuffle: T, confidenceThreshold: 0.01", side=1)
show_graph(graph4); mtext("confidenceShuffle: T, confidenceThreshold: 0.0001", side=1)

## Export the graph to cyto format ("cyto.graphml")
miic.write.network.cytoscape(g = miic.res,file = "/Users/amine/Desktop/Graphs/cyto")

## For the rest we will use the first graph: ConfidenceShuffle:T, ConfidenceThreshold: 0.001
adj.miic <- as.data.frame(adj.miic[1,1])
for (i in seq(1,dim(adj.miic)[1])) {
  if(adj.miic[i,dim(adj.miic)[2]] != 0)
    print(rownames(adj.miic)[i])
}

## Idendify the mutated genes that are significantly related to gene expression
adj.miic <- adj.miic[,-dim(adj.miic)[1]]
adj.miic <- adj.miic[-dim(adj.miic)[1],]

mutated <- c()
for (i in seq(1,dim(adj.miic)[1])) {
  row <- rownames(adj.miic)[i]
  if (tolower(row) == row)
    mutated <- c(mutated,i)
}

adj.miic <- adj.miic[ ,-mutated]
adj.miic <- adj.miic[mutated,]
## We mutated genes in rows and gene expression in columns
# So we just have to sum the rows
gene_mut <- as.data.frame(rowSums(adj.miic))
gene_mut$names <- rownames(adj.miic)
gene_mut <- gene_mut[order(gene_mut[,1],decreasing = T),]
head(gene_mut) ## ==> TP53

centrality <- cbind(as.data.frame(centr_betw(g2)$res),as.matrix(V(g2)))
colnames(centrality) <- c("centrality","i")
centrality <- centrality[order(centrality$centrality,decreasing = T),]
head(centrality)
















