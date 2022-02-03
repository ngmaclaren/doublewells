## Code by Neil MacLaren 2/2/2022

## Use the perturbation experiment as a base
## For a node that has transitioned,
### what was <x_i> over some stable period of t?
### what was k(i)?
### if y_j = 1 if x_h and 0 otherwise, does y_j = 1 for any neighbor of i?
## no early warnings per se in this experiment

library(igraph)
library(doublewells)
library(ggplot2)

whichstate <- function(x, cutoff = 7) ifelse(x > cutoff, 1, 0)

set.seed(1)

which_param <- "u" # "D"
production_runs <- TRUE # FALSE

## Graph
generate_new <- FALSE
if(generate_new) {
    nnodes <- 100
    ##cprob <- 0.06
    ##g <- sample_gnp(nnodes, cprob)
    outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
    g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
    write_graph(g, "./data/testgraph-BA.gml", format = "gml")
} else {
    g <- read_graph("./data/testgraph-BA.gml", format = "gml")
    nnodes <- vcount(g)
}
A <- as_adj(g, type = "both", sparse = FALSE)

## Node Systems
r <- c(1, 4, 7) # double well parameters
s <- 0.005 # sd of noise process
D <- 0.21 # connection strength; needs to be a lot lower for PA graphs
p <- if(s > 0) 3*s else 0.015 # perturbation strength
u <- rep(0, nnodes) # stress vector

dt <- 0.01
T <- 2000

## Perturbation
p_start <- 500
p_dur <- 10
p_stop <- p_start + p_dur

if(!production_runs) {
    ## Set/Re-Set x state
    initialx <- rep(1, nnodes) + noise(nnodes, s)
    X <- matrix(0, nrow = T, ncol = nnodes)
    x <- initialx
    
    ## Stress
    stressnode <- select_stressnode(g, "highest")
    if(which_param == "u") u[stressnode] <- 1

    ## Simulation
    for(t in 1:T) {
        X[t, ] <- x

        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
        if(t >= p_start & t < p_stop) x <- x + p
    }

    ## Analysis
    baseline <- X[250:450, ]
    endstate <- X[T, ]

    meanx <- colMeans(baseline)
    degx <- degree(g)

    V(g)$y <- whichstate(endstate)
    nbstate <- sapply(V(g), function(i) ifelse(1 %in% neighbors(g, i)$y, 1, 0))

    results <- data.frame(
        D = D, u = u[stressnode],
        meanx = meanx, degx = degx, endx = endstate, nbstate = nbstate
    )

    matplot(
        1:T, X, type = "l", lty = 1, lwd = .5, col = "steelblue",
        xlab = "t", ylab = expression(x[i]), main = paste("D =", D, "; u =", u[stressnode])
    )
    
} else {
    if(which_param == "D") {
        ##Ds <- seq(0.615, 0.635, length.out = 10)
        Ds <- seq(.2, .4, length.out = 10)
        Us <- rep(0, length(Ds))
    } else if(which_param == "u") {
        ##Us <- seq(1.84, 2.02, length.out = 10)
        ##Us <- seq(0.05, 2, length.out = 10)
        Us <- seq(.25, 1, length.out = 10)#c(.01, .025, .05, .1, .25, .5, 1, 2.5, 5, 10)
        Ds <- rep(D, length(Us))
    }

    stressnode <- select_stressnode(g, "highest")
    nsims <- 10
    results <- vector("list", length(Us)*nsims)
    sim_n <- 1

    for(i in 1:length(Us)) {
        for(j in 1:nsims) {
            D <- Ds[i]
            u[stressnode] <- Us[i]

            ## Set/Re-Set x state
            initialx <- rep(1, nnodes) + noise(nnodes, s)
            X <- matrix(0, nrow = T, ncol = nnodes)
            x <- initialx

            ## Simulation
            for(t in 1:T) {
                X[t, ] <- x

                x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
                if(t >= p_start & t < p_stop) x <- x + p
            }

            ## Analysis
            baseline <- X[250:450, ]
            endstate <- X[T, ]

            meanx <- colMeans(baseline)
            degx <- degree(g)

            V(g)$y <- whichstate(endstate)
            nbstate <- sapply(V(g), function(i) ifelse(1 %in% neighbors(g, i)$y, 1, 0))

            result <- data.frame(
                D = D, u = u[stressnode], nsim = j, node = 1:nnodes,
                meanx = meanx, degx = degx, endx = endstate, nbstate = nbstate
            )

            results[[sim_n]] <- result
            sim_n <- sim_n + 1

            ## if(j == nsims) {
            ##     dev.new()
            ##     matplot(
            ##         1:T, X, type = "l", lty = 1, lwd = .5, col = "steelblue",
            ##         xlab = "t", ylab = expression(x[i]),
            ##         main = paste("D =", D, "; u =", u[stressnode])
            ##     )
            ## }
        }
    }

    results <- do.call(rbind, results)

    results$y <- whichstate(results$endx)

    ##m1 <- glm(y ~ D + degx + nbstate, family = binomial, data = subset(results, D < 0.64))
    dev.new(width = 15, height = 10)
    ggplot(results, aes(degx, endx, color = as.factor(nbstate))) +
        geom_point(size = 3, pch = 1) +
        scale_color_brewer(palette = "Set1",
                           labels = c("No", "Yes"),
                           name = "Any Neighbors\nin Upper State?") +
        facet_wrap("u") +
        xlab("Node Degree") + ylab("State of X at End of Simulation") +
        theme_bw() +
        theme(text = element_text(size = 14))
}


