## Code by Neil G. MacLaren, 3 Jan 2021

### Try
## Fixed D, relatively low (not enough to make the system move)
## Fixed u or force on u

### Networks
## 100 nodes
## E-R: [parameters]
## Dot-product: [parameters]
## Islands: [parameters]
## B-A: [parameters]
## W-S: [parameters]
## LFR: [parameters], check documentation in NetworkX

library(igraph)
library(doublewells)

graph_choices <- c(
    "regular", "max-entropy", "sphere-surface", "islands", "pref-attach", "small-world")
graph_choice <- graph_choices[1]
nnodes <- 100

cprob <- .06 # target density approx .04; this needs to be balanced with D and probably u
if(graph_choice == "max-entropy") {
    g <- sample_gnp(nnodes, cprob)
} else if(graph_choice == "sphere-surface") {
    g <- sample_dot_product(sample_sphere_surface(dim = 3, n = nnodes, radius = .279)) # .35 is ~ .1
} else if(graph_choice == "islands") {
    nislands <- 5
    nbridges <- 1
    base_prob <- cprob*nislands
    pin <- base_prob * (100 - nbridges)/100
    g <- sample_islands(islands.n = nislands, islands.size = nnodes/nislands,
                        islands.pin = pin, n.inter = nbridges)
} else if(graph_choice == "pref-attach") {
    ## this may not add to one. Parameters are balanced by hand to achieve tgt density of .04
    outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
    g <- sample_pa(nnodes, power = 1.5,
                   ##m = 2,
                   out.dist = outdist,
                   directed = FALSE)
} else if(graph_choice == "small-world") {
    ## not fine enough control over density
    g <- sample_smallworld(dim = 1, size = 100, nei = 3, p = .1)
} else if(graph_choice == "regular") {
    g <- sample_k_regular(nnodes, 6)
}

edge_density(g)
plot(g, vertex.label = "", vertex.size = 2.5)

A <- as_adj(g, type = "both", sparse = FALSE)

r1 <- 1 # lower equil
r2 <- 2 # separatrix
r3 <- 5 # upper equil
s <- .05 # noise parameter
D <- .09 # connection strength
u <- .1  # stress parameter, add it to one node
dt <- 0.01
T <- 2000

initialx <- rep(1, nnodes)
stress <- rep(0, nnodes)
stress[sample(1:nnodes, 1)] <- u
noise <- function(s) rnorm(1, 0, s)

results <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx
for(t in 1:T) {
    results[t, ] <- x
    x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(x), stress)
}

dev.new()
matplot(1:T, results, type = "l", lty = 1, lwd = .2, col = "gray50")
