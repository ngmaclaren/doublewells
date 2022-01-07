## Code by Neil G. MacLaren, 6 Jan 2021

## use windowed_acmethod for each node
## need to set up the forcing on u; check single-system.R for a template

library(igraph)
library(doublewells)

select_stressnode <- function(g, add_stress_to = NULL) {
    "Choose a random node to which to apply stress. Can constrain choice to either 'high' or 'low' degree nodes, in which case nodes are chosen from the top or bottom quintiles of the degree distribution, respectively. Alternatively, a random node with the 'highest' or 'lowest' (non-zero) degree can be chosen."
    if(is.null(add_stress_to)) {
        nnodes <- vcount(g)
        selectnode <- sample(1:nnodes, 1)
        return(selectnode)
    } else {
        k <- degree(g)
        breaks <- quantile(k, probs = c(.2, .8))
        if(add_stress_to == "high") {
            poss <- V(g)[which(k >= breaks[2])]
        } else if(add_stress_to == "low") {
            poss <- V(g)[which(k <= breaks[1] & k > 0)]
        } else if(add_stress_to == "highest") {
            poss <- V(g)[which(k == max(k))]
            if(length(poss) == 1) return(poss)
        } else if(add_stress_to == "lowest") {
            poss <- V(g)[which(k == min(k[which(k > 0)]))]
        }
        selectnode <- sample(poss, 1)
        return(selectnode)
    }
}

graph_choices <- c(# comments are approximate critical values of D with u = 0
    "regular", # 0.298
    "max-entropy", # 0.195
    "sphere-surface", # 0.18
    "islands", # 0.194
    "pref-attach", # 0.067
    "small-world" # 0.24
)

graph_choice <- graph_choices[6]
calc_earlywarnings <- FALSE # TRUE
add_stress_to <- "lowest" # one of "highest", "high", "low", "lowest", or NULL
linear_increase <- FALSE # TRUE, FALSE, or NULL (for no increase)

nnodes <- 100
r1 <- 1 # lower equil
r2 <- 3 # separatrix
r3 <- 5 # upper equil
s <- .005 # noise parameter
D <- NULL # connection strength; can set here or let the code below set the value to just below the approximate critical threshold for each graph type.
maxU <- 2
dt <- 0.01
T <- 5000
if(is.null(linear_increase)) { # for completeness, when forcing on u is not desired
    U <- seq(0, 0, length.out = T)
} else if(linear_increase) {
    usteps <- T
    U <- seq(0, maxU, length.out = usteps)
} else {
    usteps <- 2000
    U <- c(seq(0, maxU, length.out = usteps), rep(maxU, T - usteps))
}
   
noise <- function(n, s) rnorm(n, 0, s) # the white noise function


## target density approx .06
## this needs to be balanced with D and u, as well as the small world and random regular networks, which are limited in possible values for density
cprob <- .06
if(graph_choice == "max-entropy") {
    if(is.null(D)) D <- 0.18
    g <- sample_gnp(nnodes, cprob)
    main <- "Maximum Entropy"
} else if(graph_choice == "sphere-surface") {
    if(is.null(D)) D <- 0.17
    g <- sample_dot_product(sample_sphere_surface(dim = 3, n = nnodes, radius = .279)) # .35 is ~ .1
    main <- "Sphere Surface/Dot-Product"
} else if(graph_choice == "islands") {
    if(is.null(D)) D <- 0.18
    nislands <- 5
    nbridges <- 1
    base_prob <- cprob*nislands
    pin <- base_prob * (100 - nbridges)/100
    g <- sample_islands(islands.n = nislands, islands.size = nnodes/nislands,
                        islands.pin = pin, n.inter = nbridges)
    main <- "`Islands'"
} else if(graph_choice == "pref-attach") {
    if(is.null(D)) D <- 0.06
    ## this may not add to one. Parameters are balanced by hand to achieve tgt density of .04
    outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
    g <- sample_pa(nnodes, power = 1.5,
                   ##m = 2,
                   out.dist = outdist,
                   directed = FALSE)
    main <- "Preferential Attachment"
} else if(graph_choice == "small-world") {
    if(is.null(D)) D <- 0.23
    ## not fine enough control over density
    g <- sample_smallworld(dim = 1, size = 100, nei = 3, p = .1)
    main <- "Small World"
} else if(graph_choice == "regular") {
    if(is.null(D)) D <- 0.29
    g <- sample_k_regular(nnodes, 6)
    main <- "Random Regular"
}

A <- as_adj(g, type = "both", sparse = FALSE)

initialx <- rep(1, nnodes)
results <- matrix(0, nrow = T, ncol = nnodes)

stress <- rep(0, nnodes)
stressnode <- select_stressnode(g, add_stress_to = add_stress_to) # "lowdegree")

x <- initialx

for(t in 1:T) {
    results[t, ] <- x
    
    stress[stressnode] <- U[t] # u
    x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(nnodes, s), stress)
}

## Sim Results
dev.new(width = 20, height = 10)
par(mfrow = c(1, 2))
V(g)$nodestate <- x/max(x)
colorfun <- colorRamp(c("blue", "orange"), space = "Lab")
nodecolor <- colorfun(V(g)$nodestate)
V(g)$color <- apply(nodecolor, 1, function(x) rgb(t(x), maxColorValue = 255))
nodesize <- ifelse(stress == 0, 4, 8)
linewidth <- ifelse(stress == 0, .2, 2)
plot(
    g, vertex.label = "",  main = main,
    vertex.size = nodesize
)
matplot(
    1:T, results, type = "l",
    lty = 1,
    lwd = linewidth,
    col = V(g)$color,
    xlab = "t", ylab = expression(x[i]), main = "Node States"
)

## Early Warning
if(calc_earlywarnings) {
    acresults <- apply(results, 2, function(x) windowed_acmethod(x, wl = T/4))

    dev.new()
    matplot(
        1:nrow(acresults), acresults, type = "l",
        lty = 1,
        lwd = linewidth,
        col = V(g)$color,
        xlab = "Window", ylab = expression("AR Coefficient"),
        main = "Node-Level Early Warning Signals"
    )
}

print(paste0("Degree of stress node is ", degree(g, V(g)[stressnode])))
