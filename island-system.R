## Code by Neil G. MacLaren, Dec 30, 2021

## Simulation inspired by Drake & Griffen 2010
## Attempting a modular network
## Try:
### variance across x_i, w/n time step
### skewness, same
### power spectrum --  not sure how yet
### Chen L solution: cor(X) w/n cluster, sd(X) w/n cluster, cor(X_c1, X_c2) between clusters

library(igraph)
library(doublewells)

### Network parameters
n_islands <- 5
size_islands <- 100
edge_prob <- .2
bridges <- 1

### Double-well system parameters
r1 <- 1
r2 <- 2
r3 <- 5

## make the network
g <- sample_islands(islands.n = n_islands, islands.size = size_islands, islands.pin = edge_prob, n.inter = bridges)
V(g)$membership <- rep(1:5, each = 100)
##plot(g, vertex.size = 0)
##C <- cluster_edge_betweenness(g) # Does this run as O(n^3)?
##C <- cluster_walktrap(g)
##V(g)$membership <- C$membership
##plot(C, g)
A <- as_adjacency_matrix(g, sparse = FALSE)

### Simulation parameters
## time
dt <- 0.01
T <- 4000#2.5e5#20000

## initial state
initialx <- rep(r1, vcount(g))

## noise
noise_level <- .1
noise <- function(n, noise_level) rnorm(n = n, mean = 0, sd = noise_level)

## forcing
## nsteps <- 20
## sl <- floor(T/nsteps)#2500#1000
## Ds <- rep(.2, length.out = sl)#c(rep(.1, 5), seq(.1, .35, length.out = nsteps-5))


### Simulation

## set up simulation
## niter <- 1 # probably want 50 or 100 for production runs
##p <- 1
##D <- Ds[p]
D <- .025
results <- matrix(0, ncol = length(initialx), nrow = T)
x <- initialx

## run simulation
for(t in 1:T) {
    ## if(t %% sl == 0 & t != T) {
    ##     p <- p + 1
    ##     D <- Ds[p]
    ## }
    
    results[t, ] <- x
    x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
}

### Results
plot(1:T, rowMeans(results))

pal <- rainbow(max(V(g)$membership))
colors <- sapply(V(g)$membership, function(x) pal[x])
matplot(1:T, results, type = "l", lwd = .5, lty = 1, col = colors)

dev.new(); plot(1:T, apply(results, 1, sd), type = "l", lty = 1, lwd = 1)
dev.new(); plot(1:T, apply(results, 1, moments::skewness), type = "l", lty = 1, lwd = 1)


times <- seq(1, T, by = 100)
samp <- results[times, ]
endpoints <- 1:n_islands * size_islands
startpoints <- (0:(n_islands - 1) * size_islands) + 1
clusters <- mapply(seq, startpoints, endpoints, SIMPLIFY = FALSE)

sds <- vector("list", length(clusters))
for(i in 1:length(clusters)) sds[[i]] <- apply(samp[, clusters[[i]]], 1, sd)

## for the correlations, I need the within time step correlation of each cluster's states with each other cluster's states

cor_wn <- vector("list", length(clusters))
for(i in 1:length(clusters)) cor_wn[[i]] <- apply(samp[, clusters[[i]]], 1, cor)
