## Core code by Neil G. MacLaren, 6 Jan 2021
## Updated throughout Jan 2022

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

Iadj <- function(X, A, t = NULL, nsteps = NULL, times = NULL) {
    ## A is a 2D array of connection weights (all w_ij ∈ {0, 1} for now)
    ## t is the chosen time
    ## X is a 2D array where the rows are t_t and the columns are x_i
    stopifnot((length(t) == 1 & length(nsteps) == 1) | length(times) > 1)
    if(length(times) == 0) times <- seq(t - nsteps + 1, t)

    X <- X[times, ]
    N <- nrow(A)
    W <- sum(A)

    muI <- colMeans(X)
    deltaX <- apply(X, 2, function(x) x[nsteps] - mean(x))
    
    numerator <- 0
    for(i in 1:nrow(A)) {
        for(j in 1:ncol(A)) {
            if(j >= i) next
            y <- A[i, j]*deltaX[i]*deltaX[j]
            numerator <- numerator + y
        }
    }
    denomenator <- sum((X[nrow(X), ] - muI)^2)

    (N/W)*(numerator/denomenator)
}

graph_choices <- c(
    "regular",
    "max-entropy",
    "sphere-surface",
    "islands",
    "pref-attach",
    "small-world"
)

graph_choice <- graph_choices[2]
calc_earlywarnings <- FALSE # TRUE
add_stress_to <- "highest" # one of "highest", "high", "low", "lowest", or NULL
linear_increase <- TRUE # TRUE, FALSE, or NULL (for no increase)
increase_D <- FALSE
increase_u <- TRUE

cutoff_value <- 2 # ideally this would be the "point of no return" value
plot_transition <- TRUE
smooth_ts <- FALSE

nnodes <- 100
r1 <- 1 # lower equil
r2 <- 4 # separatrix
r3 <- 7 # upper equil
s <- 0.01 # noise parameter
maxD <- NULL # connection strength; can set here or let the code below set the value to just below the approximate critical threshold for each graph type.
maxU <- NULL # stress/bias; same as for maxD
dt <- 0.01
T <- 5000
 
## target density approx .06
## this needs to be balanced with D and u, as well as the small world and random regular networks, which are limited in possible values for density
cprob <- .06
if(graph_choice == "regular") {
    if(is.null(maxD)) maxD <- 0.9#0.84
    if(is.null(maxU)) maxU <- 1.4
    g <- sample_k_regular(nnodes, 6)
    main <- "Random Regular"
} else if(graph_choice == "max-entropy") {
    if(is.null(maxD)) maxD <- 0.6#0.66#0.54
    if(is.null(maxU)) maxU <- 4#2.7
    g <- sample_gnp(nnodes, cprob)
    main <- "Maximum Entropy"
} else if(graph_choice == "sphere-surface") {
    if(is.null(maxD)) maxD <- 0.63#0.56
    if(is.null(maxU)) maxU <- 2.5
    g <- sample_dot_product(sample_sphere_surface(dim = 3, n = nnodes, radius = .279)) # .35 is ~ .1
    main <- "Sphere Surface/Dot-Product"
} else if(graph_choice == "islands") {
    if(is.null(maxD)) maxD <- 0.68#0.58
    if(is.null(maxU)) maxU <- 4.5
    nislands <- 5
    nbridges <- 1
    base_prob <- cprob*nislands
    pin <- base_prob * (100 - nbridges)/100
    g <- sample_islands(islands.n = nislands, islands.size = nnodes/nislands,
                        islands.pin = pin, n.inter = nbridges)
    main <- "'Islands'"
} else if(graph_choice == "pref-attach") {
    if(is.null(maxD)) maxD <- 0.35#0.18
    if(is.null(maxU)) maxU <- 10.5
    ## this may not add to one. Parameters are balanced by hand to achieve tgt density of .04
    outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
    g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
    main <- "Preferential Attachment"
} else if(graph_choice == "small-world") {
    if(is.null(maxD)) maxD <- 0.8#0.73
    if(is.null(maxU)) maxU <- 2.7
    ## not fine enough control over density
    g <- sample_smallworld(dim = 1, size = 100, nei = 3, p = .1)
    main <- "Small World"
}

A <- as_adj(g, type = "both", sparse = FALSE)

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

if(increase_D) {
    Ds <- seq(0, maxD, length.out = T)
} else {
    D <- maxD*.85 # .8 for ME, .85 for sph-s and islands, 
}

stress <- rep(0, nnodes)
if(increase_u) {
    stressnode <- select_stressnode(g, add_stress_to = add_stress_to) # "lowdegree")
}

initialx <- rep(1, nnodes)
results <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx
for(t in 1:T) {
    results[t, ] <- x

    if(increase_D) {
        D <- Ds[t]
    }

    if(increase_u) stress[stressnode] <- U[t] # u
    x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(nnodes, s), stress)
}

if(plot_transition) {
    first_node <- which(apply(results, 1, function(x) sum(x > cutoff_value)) > 0)[1]
    t_time <- which(apply(results, 1, function(x) sum(x > cutoff_value)/length(x)) > .3)[1]

    p_transitioned <- apply(results, 1, function(x) sum(x > cutoff_value)/length(x))
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

## print(paste0("Degree of stress node is ", degree(g, V(g)[stressnode])))

## ##ape::Moran.I(results[nrow(results), ], A, scaled = TRUE)$observed
## mi <- apply(results[2:T, ], 1, function(x) ape::Moran.I(x, A, scaled = TRUE)$observed)
## dev.new()
## plot(2:T, mi, type = "l", main = "I", xlab = "t", ylab = "I")

## ## This is I' applied across the network.
## nsteps <- 5
## ts <- seq(nsteps + 1, T)
## miadj <- sapply(ts, function(t) Iadj(results, A, t = t, nsteps = nsteps))
## dev.new()
## plot(ts, miadj, type = "l", main = "I'", xlab = "t", ylab = "I'")

## dev.new()
## lwds <- ifelse(1:100 %in% nbs, 2, .5)
## ltys <- ifelse(1:100 %in% nbs, 2, 1)
## matplot(1:T, results, type = "l", lty = ltys, lwd = lwds, col = "black")

## I want I' applied only to the neighbors of stressnode

## stressnode_ew <- earlywarnings::generic_ews(results[results[, stressnode] < 2.5, stressnode],
##                                             winsize = 50, detrending = "no")

##nbs <- neighborhood(g, order = 2, nodes = stressnode)[[1]]

if(increase_u) {
    dev.new(width = 18, height = 12); par(mfrow = c(2, 3))
} else {
    dev.new(width = 14, height = 14); par(mfrow = c(2, 2))
}

## 1
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

## 2
par(xpd = FALSE)
matplot(
    1:T, results, type = "l",
    lty = 1,
    lwd = linewidth,
    col = V(g)$color,
    xlab = "t", ylab = expression(x[i]), main = "Node States"
)
abline(h = cutoff_value, col = "black", lwd = 2, lty = 3)
abline(v = first_node, col = "red1", lwd = 2, lty = 3)
abline(v = t_time, col = "red3", lwd = 2, lty = 3)

## 3
ts <- 2:T
mi <- apply(results[ts, ], 1, function(x) ape::Moran.I(x, A, scaled = TRUE)$observed)
plot(ts, mi, type = "l", main = "Network I (Original)", xlab = "t", ylab = "I",
     col = adjustcolor("darkgray", alpha.f = .8), lwd = 1)#, xlim = c(1, T), ylim = c(-.5, .5))
if(smooth_ts) {
    smooth <- loess(mi ~ ts, span = .1)
    smoothed <- predict(smooth)
    lines(ts, smoothed, col = "darkorchid", lwd = 2, lty = 2)
}
abline(v = first_node, col = "red1", lwd = 2, lty = 3)
abline(v = t_time, col = "red3", lwd = 2, lty = 3)
lines(1:T, p_transitioned, lwd = 2, lty = 1, col = "goldenrod")
## par(new = TRUE)
## plot(1:T, p_transitioned, type = "l", lwd = 2, lty = 1, col = "goldenrod",
##      axes = FALSE, xlab = "", ylab = "")
## axis(side = 4, at = pretty(range(p_transitioned)), col.ticks = "goldenrod", col.axis = "goldenrod")

## 4
nsteps <- 5
##offset <- 50
ts <- seq(nsteps + 1, T)#seq(offset + 1, T)
miadj <- sapply(ts, function(t) Iadj(results, A, t = t, nsteps = nsteps))
##miadj <- sapply(ts, function(t) Iadj(results, A, times = seq(t - offset, t, by = 10)))##t = t, nsteps = nsteps))
##par(mfrow = c(1, 3))
plot(ts, miadj, type = "l", main = "Network I'", xlab = "t", ylab = "I'",
     col = adjustcolor("darkgray", alpha.f = .8), lwd = 1)#, xlim = c(1, T), ylim = c(-.5, .5))
if(smooth_ts) {
    smooth <- loess(miadj ~ ts, span = .1)
    smoothed <- predict(smooth)
    lines(ts, smoothed, col = "dodgerblue", lwd = 2, lty = 2)
}
abline(v = first_node, col = "red1", lwd = 2, lty = 3)
abline(v = t_time, col = "red3", lwd = 2, lty = 3)
lines(1:T, p_transitioned, lwd = 2, lty = 1, col = "goldenrod")
## par(new = TRUE)
## plot(1:T, p_transitioned, type = "l", lwd = 2, lty = 1, col = "goldenrod",
##      axes = FALSE, xlab = "", ylab = "")
## axis(side = 4, at = pretty(range(p_transitioned)), col.ticks = "goldenrod", col.axis = "goldenrod")

if(increase_u) {
    ## 5
    nbs <- neighbors(g, stressnode)
    miadj <- sapply(ts, function(t) Iadj(results[, nbs], A[nbs, nbs], t = t, nsteps = nsteps))
    plot(ts, miadj, type = "l", main = "Neighbors I'", xlab = "t", ylab = "I'",
         col = adjustcolor("darkgray", alpha.f = .8), lwd = 1)#, xlim = c(1, T), ylim = c(-.5, .5))
    if(smooth_ts) {
        smooth <- loess(miadj ~ ts, span = .1)
        smoothed <- predict(smooth)
        lines(ts, smoothed, col = "dodgerblue", lwd = 2, lty = 2)
    }
    abline(v = first_node, col = "red1", lwd = 2, lty = 3)
    abline(v = t_time, col = "red3", lwd = 2, lty = 3)
    lines(1:T, p_transitioned, lwd = 2, lty = 1, col = "goldenrod")
    ## par(new = TRUE)
    ## plot(1:T, p_transitioned, type = "l", lwd = 2, lty = 1, col = "goldenrod",
    ##      axes = FALSE, xlab = "", ylab = "")
    ## axis(side = 4, at = pretty(range(p_transitioned)), col.ticks = "goldenrod", col.axis = "goldenrod")
    
    ## 6
    nbhd <- neighborhood(g, order = 1, nodes = stressnode)[[1]]
    miadj <- sapply(ts, function(t) Iadj(results[, nbhd], A[nbhd, nbhd], t = t, nsteps = nsteps))
    plot(ts, miadj, type = "l", main = "Neighborhood I'", xlab = "t", ylab = "I'",
         col = adjustcolor("darkgray", alpha.f = .8), lwd = 1)#, xlim = c(1, T), ylim = c(-.5, .5))
    if(smooth_ts) {
        smooth <- loess(miadj ~ ts, span = .1)
        smoothed <- predict(smooth)
        lines(ts, smoothed, col = "dodgerblue", lwd = 2, lty = 2)
    }
    abline(v = first_node, col = "red1", lwd = 2, lty = 3)
    abline(v = t_time, col = "red3", lwd = 2, lty = 3)
    lines(1:T, p_transitioned, lwd = 2, lty = 1, col = "goldenrod")
    ## par(new = TRUE)
    ## plot(1:T, p_transitioned, type = "l", lwd = 2, lty = 1, col = "goldenrod",
    ##      axes = FALSE, xlab = "", ylab = "")
    ## axis(side = 4, at = pretty(range(p_transitioned)), col.ticks = "goldenrod", col.axis = "goldenrod")
    
}
## for(v in nbhd) {
##     ew <- earlywarnings::generic_ews(results[results[, v] < 2.5, v],
##                                       winsize = 50, detrending = "no")
## }


## anb <- sample(nbs, 1)
## aneighbor_ew <-  earlywarnings::generic_ews(results[results[, anb] < 2.5, anb],
##                                             winsize = 50, detrending = "no")
