## Using this more complete simulation set up to explore the different features of the simulation more systematically
## Some features that matter, besides D and u, are (1) network structure or class, (2) different instances of that class, and (3) which node within the network is chosen.
## Let's first try fixing the network and applying u to several different nodes.

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

noise <- function(n, s, f = rnorm) {
    "A function to generate `n` random noise values from a distribution with standard deviation `s`. Currently set up only for Gaussian noise, but could be expanded to support more distributions."
    f(n, 0, s)
}

graph_choices <- c(# comments are approximate critical values of D with u = 0
    "regular", # 0.298
    "max-entropy", # 0.185
    "sphere-surface", # 0.18
    "islands", # 0.194
    "pref-attach", # 0.067
    "small-world" # 0.24
)

graph_choice <- graph_choices[2]
calc_earlywarnings <- FALSE # TRUE
add_stress_to <- NULL # one of "highest", "high", "low", "lowest", or NULL
linear_increase <- TRUE # TRUE, FALSE, or NULL (for no increase)

nnodes <- 100
r1 <- 1 # lower equil
r2 <- 3 # separatrix
r3 <- 5 # upper equil
s <- 0#.005 # noise parameter
D <- NULL # connection strength; can set here or let the code below set the value to just below the approximate critical threshold for each graph type.
maxU <- 2.1 # stress/bias; same as for D
dt <- 0.01
T <- 5000

nnets <- 4
netresults <- vector("list", nnets)
allstressnodes <- vector("list", nnets)
colors <- c("steelblue", "darkorange", "orchid", "firebrick")
for(net in 1:nnets) {
    ## target density approx .06
    ## this needs to be balanced with D and u, as well as the small world and random regular networks, which are limited in possible values for density
    cprob <- .06
    if(graph_choice == "regular") {
        if(is.null(D)) D <- 0.29
        if(is.null(maxU)) maxU <- 0.7
        g <- sample_k_regular(nnodes, 6)
        main <- "Random Regular"
    } else if(graph_choice == "max-entropy") {
        if(is.null(D)) D <- 0.18
        if(is.null(maxU)) maxU <- 2
        g <- sample_gnp(nnodes, cprob)
        main <- "Maximum Entropy"
    } else if(graph_choice == "sphere-surface") {
        if(is.null(D)) D <- 0.17
        if(is.null(maxU)) maxU <- 2.5
        g <- sample_dot_product(sample_sphere_surface(dim = 3, n = nnodes, radius = .279)) # .35 is ~ .1
        main <- "Sphere Surface/Dot-Product"
    } else if(graph_choice == "islands") {
        if(is.null(D)) D <- 0.18
        if(is.null(maxU)) maxU <- 2.5
        nislands <- 5
        nbridges <- 1
        base_prob <- cprob*nislands
        pin <- base_prob * (100 - nbridges)/100
        g <- sample_islands(islands.n = nislands, islands.size = nnodes/nislands,
                            islands.pin = pin, n.inter = nbridges)
        main <- "`Islands'"
    } else if(graph_choice == "pref-attach") {
        if(is.null(D)) D <- 0.06
        if(is.null(maxU)) maxU <- 2.5
        ## this may not add to one. Parameters are balanced by hand to achieve tgt density of .04
        outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
        g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
        main <- "Preferential Attachment"
    } else if(graph_choice == "small-world") {
        if(is.null(D)) D <- 0.23
        if(is.null(maxU)) maxU <- 2.5
        ## not fine enough control over density
        g <- sample_smallworld(dim = 1, size = 100, nei = 3, p = .1)
        main <- "Small World"
    }

    ## For setting the increase in u
    if(is.null(linear_increase)) {# for completeness, when forcing on u is not desired
        U <- seq(0, 0, length.out = T)
    } else if(linear_increase) {
        usteps <- T
        U <- seq(0, maxU, length.out = usteps)
    } else {# stable period at the beginning, a period of increase to the max, then another stable period
        usteps <- 2000
        U <- c(rep(0, (T - usteps)/2), seq(0, maxU, length.out = usteps), rep(maxU, (T - usteps)/2))
    }

    A <- as_adj(g, type = "both", sparse = FALSE) # as of here, the graph is fixed, as is the progression in u. 

    ## let's adjust here
    ntrials <- 6
    resultlist <- vector("list", length = ntrials)
    stressnodes <- numeric(ntrials)

    for(trial in 1:ntrials) {
        ## This part is unchanged
        stress <- rep(0, nnodes)
        stressnode <- select_stressnode(g, add_stress_to = add_stress_to)
        stressnodes[trial] <- degree(g, V(g)[stressnode])

        initialx <- rep(1, nnodes)
        results <- matrix(0, nrow = T, ncol = nnodes)
        x <- initialx
        for(t in 1:T) {
            results[t, ] <- x
            
            stress[stressnode] <- U[t] # u
            x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(nnodes, s), stress)
        }

        ## store the results
        resultlist[[trial]] <- results
    }

    netresults[[net]] <- resultlist
    allstressnodes[[net]] <- stressnodes
}

## Now here I want a paneled figure, the panels of which will be the matplots like below
dev.new(height = 15, width = 20)
par(mfrow = c(nnets, ntrials))
for(net in 1:nnets) {
    for(trial in 1:ntrials) {
        matplot(
            1:T, netresults[[net]][[trial]], type = "l",
            lty = 1, lwd = .2, col = colors[net],
            xlim = c(0, T), ylim = c(1, 6.5),
            xlab = "t", ylab = expression(x[i]),
            main = paste0("Network ", net, ", Trial ", trial,
                          "; k_u = ", allstressnodes[[net]][[trial]])
        )
    }
}

## ## Sim Results
## dev.new(width = 20, height = 10)
## par(mfrow = c(1, 2))
## V(g)$nodestate <- x/max(x)
## colorfun <- colorRamp(c("blue", "orange"), space = "Lab")
## nodecolor <- colorfun(V(g)$nodestate)
## V(g)$color <- apply(nodecolor, 1, function(x) rgb(t(x), maxColorValue = 255))
## nodesize <- ifelse(stress == 0, 4, 8)
## linewidth <- ifelse(stress == 0, .2, 2)
## plot(
##     g, vertex.label = "",  main = main,
##     vertex.size = nodesize
## )
## matplot(
##     1:T, results, type = "l",
##     lty = 1,
##     lwd = linewidth,
##     col = V(g)$color,
##     xlab = "t", ylab = expression(x[i]), main = "Node States"
## )

## ## print(paste0("Degree of stress node is ", degree(g, V(g)[stressnode])))
