library(igraph)

## hyperbolic_graph_generator too opinionated? Maybe try very high temperatures.
## try the dot product model graphs from igraph
P <- matrix(runif(9), ncol = 3)
P <- apply(P, 2, function(x) abs(x)/sqrt(sum(x^2)))

edge_probs <- function(X) {
    ## where x is a matrix of positions
    ncols <- ncol(X)
    ep <- expand.grid(1:ncols, 1:ncols)
    ep <- ep[!duplicated(ep), ]
    ep <- ep[ep[, 1] != ep[, 2], ]
    ep[, 3] <- apply(ep, 1, function(x) X[, x[1]] %*% X[, x[2]])
    ep
}

ep <- expand.grid(1:3, 1:3)
ep <- cols[!duplicated(cols), ]
ep <- cols[!(cols[,1] == cols[,2]), ]
ep$p <- apply(ep, 1, function(x) P[, x[1]] %*% P[, x[2]]) 
ep

nnodes <- 100
ndim <- 8
P <- sample_dirichlet(nnodes, rep(1.5, ndim))
##P
##hist(P)
##edge_probs(P)
g <- sample_dot_product(P)
plot(g)


## drichilet is too uniform
##X <- sample_dirichlet(n = 100, alpha = runif(n = 10, min = .5, max = 2))
##X <- sample_sphere_surface(dim = 10, n = 500, radius = .1)
X <- sample_sphere_volume(dim = 2, n = 100, radius = .5)
##pos <- t(X)
##plot(pos)
g <- sample_dot_product(X)
wtc <- cluster_walktrap(g); plot(wtc, g)
lc <- cluster_louvain(g); plot(lc, g)
##plot(g, vertex.label = "", vertex.size = 3)
kcore <- coreness(g)
cores <- sort(unique(kcore))
pal <- RColorBrewer::brewer.pal(length(cores), "YlOrRd")
V(g)$color <- sapply(kcore, function(x) pal[which(x == cores)])
plot(g, vertex.label = "", vertex.size = 3)



edge_density(g)
transitivity(g)

degs <- degree(g)
summary(degs)
hist(degs)




## Try islands
g <- sample_islands(islands.n = 5, islands.size = 15, islands.pin = .3, n.inter = 2)
plot(g, vertex.size = 1)
wtc <- cluster_walktrap(g); plot(wtc, g)
lc <- cluster_louvain(g); plot(lc, g)
imc <- cluster_infomap(g); plot(imc, g)
ebc <- cluster_edge_betweenness(g); plot(ebc, g)
