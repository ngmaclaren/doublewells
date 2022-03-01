## Code by Neil MacLaren 2/23/2022

library(igraph)
library(doublewells)
MoranI <- ape::Moran.I

##set.seed(123)

save_plots <- FALSE

### Graph Selection
choice <- "pref_attach"
if(choice %in% c("dolphins", "netscience-lcc", "lfr-n100k6-decayrate")) {
    filename <- paste0("./data/", choice, ".mat")
    edgelist <- as.matrix(read.table(file = filename, header = FALSE, sep = " ", skip = 1)[, 1:2])
    g <- graph_from_edgelist(edgelist, directed = FALSE)
} else {
    data(list = choice)
    g <- get(choice)
}

N <- vcount(g)
A <- as_adj(g, type = "both", sparse = FALSE)

### Node Systems
K <- 10 # carrying capacity, currently uniform for all nodes
h <- 1.79 # harvesting rate, currently uniform for all nodes; this is like u (stress)
D <- 0.19 # dispersal rate, currently uniform for all nodes; this is like D (connection strength)
s <- 0.02/sqrt(10)# intensity (sd) for noise process, currently Gaussian

### Simulation settings
initialx <- rep(1, N) + noise(N, s)

τ <- 100
dt <- 0.001
T <- τ/dt
s <- s*sqrt(dt) 

## Main Loop
x <- initialx
X <- matrix(0, nrow = T, ncol = N)
for(t in 1:T) {
    X[t, ] <- x
    x <- metapopulation(x, K, h, A, D, dt, noise(N, s))
}

plotsamples <- seq(1, nrow(X), by = τ) # separates plotting samples from integration
plotX <- X[plotsamples, ]
matplot(plotsamples, plotX, type = "p", col = "black", pch = 1, cex = .5,
        xlab = "t", ylab = expression(x[i]), xaxt = "n")
xticklabels <- pretty(range(plotsamples))*dt
axis(1, at = pretty(range(plotsamples)), labels = xticklabels)

