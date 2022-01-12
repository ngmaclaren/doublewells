## I want to know how maximum entropy (Erdos-Renyi) graphs vary in terms of features that might be relevant to the other simulations we're doing.
library(igraph)

## I will match the basic graph type to what I'm using in the 100 Nodes simulations
nnodes <- 100
cprob <- 0.06
nislands <- 5
nbridges <- 1
base_prob <- cprob*nislands
pin <- base_prob * (100 - nbridges)/100
outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)

nsims <- 1000
maxeigs <- numeric(nsims)
gaps <- numeric(nsims)
ccs <- numeric(nsims)
for(sim in 1:nsims) {
    ## g <- sample_k_regular(nnodes, 6)
    g <- sample_gnp(nnodes, cprob, directed = FALSE) # for production runs, make many of these
    ## g <- sample_dot_product(sample_sphere_surface(dim = 3, n = nnodes, radius = .279))
    ## g <- sample_islands(islands.n = nislands, islands.size = nnodes/nislands, islands.pin = pin, n.inter = nbridges)
    ## g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
    ## g <- sample_smallworld(dim = 1, size = 100, nei = 3, p = .1)
    A <- as_adj(g, type = "both", sparse = FALSE) # need the adj matrix for spectral properties
    ##L <- laplacian_matrix(g, sparse = FALSE) # dangerous because L means something in R (of course, so does T and I use that)
    Λ <- eigen(A, only.values = TRUE)$values
    maxeigs[sim] <- Λ[1]
    gaps[sim] <- Λ[1] - Λ[2]
    ccs[sim] <- transitivity(g, type = "localaverageundirected")
}

dev.new()
par(mfrow = c(3, 1))
nbreaks <- 30
hist(maxeigs, breaks = nbreaks)
hist(gaps, breaks = nbreaks)
hist(ccs, breaks = nbreaks)

