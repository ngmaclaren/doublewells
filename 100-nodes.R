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
graph_choice <- graph_choices[4]

nnodes <- 100
r1 <- 1 # lower equil
r2 <- 2 # separatrix
r3 <- 5 # upper equil
s <- .01 # noise parameter
D <- .06 # connection strength
u <- .01  # stress parameter, add it to one node
dt <- 0.01
T <- 2000


cprob <- .06 # target density approx .04; this needs to be balanced with D and probably u
if(graph_choice == "max-entropy") {
    g <- sample_gnp(nnodes, cprob)
    main <- "Maximum Entropy"
} else if(graph_choice == "sphere-surface") {
    g <- sample_dot_product(sample_sphere_surface(dim = 3, n = nnodes, radius = .279)) # .35 is ~ .1
    main <- "Sphere Surface/Dot-Product"
} else if(graph_choice == "islands") {
    nislands <- 5
    nbridges <- 1
    base_prob <- cprob*nislands
    pin <- base_prob * (100 - nbridges)/100
    g <- sample_islands(islands.n = nislands, islands.size = nnodes/nislands,
                        islands.pin = pin, n.inter = nbridges)
    main <- "`Islands'"
} else if(graph_choice == "pref-attach") {
    ## this may not add to one. Parameters are balanced by hand to achieve tgt density of .04
    outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
    g <- sample_pa(nnodes, power = 1.5,
                   ##m = 2,
                   out.dist = outdist,
                   directed = FALSE)
    main <- "Preferential Attachment"
} else if(graph_choice == "small-world") {
    ## not fine enough control over density
    g <- sample_smallworld(dim = 1, size = 100, nei = 3, p = .1)
    main <- "Small World"
} else if(graph_choice == "regular") {
    g <- sample_k_regular(nnodes, 6)
    main <- "Random Regular"
}

A <- as_adj(g, type = "both", sparse = FALSE)

initialx <- rep(1, nnodes)
stress <- rep(0, nnodes)
stressnode <- sample(1:nnodes, 1)
stress[stressnode] <- u
noise <- function(s) rnorm(nnodes, 0, s)

results <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx
for(t in 1:T) {
    results[t, ] <- x
    x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(s), stress)
}

dev.new(width = 20, height = 10)
par(mfrow = c(1, 2))
V(g)$nodestate <- x/max(x)
colorfun <- colorRamp(c("blue", "orange"), space = "Lab")
nodecolor <- colorfun(V(g)$nodestate)
V(g)$color <- apply(nodecolor, 1, function(x) rgb(t(x), maxColorValue = 255))
##nodecolor <- ifelse(stress == 0, "orange", "dodgerblue")
nodesize <- ifelse(stress == 0, 3, 9)
linewidth <- ifelse(stress == 0, .2, 2)
plot(g, vertex.label = "",  main = main,
     vertex.size = nodesize
     ##vertex.color = nodecolor
     )
matplot(1:T, results, type = "l",
        lty = 1,
        lwd = linewidth,
        col = V(g)$color,
        xlab = "t", ylab = expression(x[i]), main = "Node States")
