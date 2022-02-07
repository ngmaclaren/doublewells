## Code by Neil MacLaren 2/7/2022

## 1. covariance/eigenvalue method with decreasing numbers of nodes
##    Try also removing nodes from the covariance matrix as they transition to the upper state
##    Here, just collect a matrix like this, X[<window>, <nodes>], then calculate its covariance matrix, and from there the dominant eigenvalue thereof
sampled_eigenmethod <- function(X, A, samples, nodes) {
    "From an output matrix, X, with connections between x_i given in A, take x_i at `samples` points in time and return the dominant eigenvalue of the resulting covariance matrix."
    X <- X[samples, nodes]
    A <- A[nodes, nodes]
    C <- cov(X)
    eig <- eigen(C, symmetric = TRUE, only.values = TRUE)[[1]][1]
    return(eig)
}
## 2. Moran's I with node "deletion" as well
##    Leave the "real" network in place, but remove nodes (adjusting N and W as well) from what goes into the Moran's I calculation
##    I think might be able to just adjust the ape::Moran.I() call by removing nodes inside the call, both to X[stepT, <here>] and to A[<here>, <here>].
##    ape::Moran.I(X[stepT, nodes], A[nodes, nodes])$observed
sampled_MoranI <- function(X, A, nodes, t) ape::Moran.I(X[t, nodes], A[nodes, nodes])$observed

## 3. sd(x_i) with sentinel nodes
##    e.g., start by looking at the highest degree nodes, then next highest degree, and so on
##    In scale-free network may need to include whether or not i is a neighbor of a transitioned node
##    This one will take some more work, maybe. First, just find the degree distribution, then watch only the highest degree nodes with the regular sd(x_i)[<window>] method.
sampled_sdmethod <- function(X, t, wl) apply(X[(wl+1):t, ], 2, sd)

choose_sentinels <- function(g, n, v = V(g)) { ## <-- This function is not doing what I want
    ## g is the graph, n is how many sentinels
    deg <- sort(degree(g, v = v), decreasing = TRUE)
    maxdegs <- head(deg, n)
    nodes <- sample(V(g)[which(degree(g) %in% maxdegs)], n)
    V(g)[nodes]
}

library(igraph)
library(doublewells)

whichstate <- function(x, cutoff = 7) ifelse(x >= cutoff, 1, 0)
movingup <- function(x, cutoff = 2.5) ifelse(x >= cutoff, 1, 0)

set.seed(123)

which_net <- "pa" # "me"
which_param <- "D" # "u"
which_ew <- "sd" # "ac"
save_plot <- FALSE # TRUE
##if(which_net == "me") cutoff <- .25*nnodes else cutoff <- .75*nnodes

stopmsg <- "Will not converge with param settings: maybe 'u' needs to be applied to more than one node."
if(which_net == "pa" & which_param == "u") stop(stopmsg)

## Graph
generate_new <- FALSE 
if(generate_new) {
    nnodes <- 100
    if(which_net == "me") {
        cprob <- 0.06
        g <- sample_gnp(nnodes, cprob)
        write_graph(g, "./data/testgraph.gml", format = "gml")
    } else if(which_net == "pa") {
        outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
        g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
        write_graph(g, "./data/testgraph-BA.gml", format = "gml")
    }
} else {
    if(which_net == "me") {
        g <- read_graph("./data/testgraph.gml", format = "gml")
    } else if(which_net == "pa") {
        g <- read_graph("./data/testgraph-BA.gml", format = "gml")
    }
    nnodes <- vcount(g)
}
A <- as_adj(g, type = "both", sparse = FALSE)

## Node Systems
r <- c(1, 4, 7) # double well parameters
s <- 0.005 # sd of noise process
### the below may need to be fixed to support varying u
if(which_net == "me") D <- 0.54 else if(which_net == "pa") D <- 0.21 # connection strength
p <- if(s > 0) 3*s else 0.015 # perturbation strength
u <- rep(0, nnodes) # stress vector

dt <- 0.01
##T <- 2000

initialx <- rep(1, nnodes) + noise(nnodes, s)
##X <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx
    
## Stress
stressnode <- select_stressnode(g, "highest")
if(which_param == "u") {
    if(which_net == "me") {
        u[stressnode] <- 1.8 # with D = 0.54
    } else if(which_param == "pa") {
        u[stressnode] <- 0.3
    }
}

## Simulation
cutoff <- .25*nnodes
atrisk <- 1:nnodes
stepT <- 500
wl <- 250#.5*stepT
n_sentinels <- 5
sentinels <- choose_sentinels(g, n = n_sentinels, v = atrisk)
##ews <- list()
Ds <- numeric()
Us <- numeric()
eigs <- numeric()
fulleigs <- numeric()
moranIs <- numeric()
fullIs <- numeric()
ssds <- list()
##Is <- numeric()
i <- 1
n_atrisk <- numeric()
##for(i in 1:1) {
while(length(atrisk) > cutoff) {
    X <- matrix(0, nrow = stepT, ncol = nnodes)

    for(t in 1:stepT) {
        X[t, ] <- x
        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
    }

    ## Analysis
    atrisk <- 1:nnodes
    upwardbound <- movingup(x)
    atrisk <- atrisk[which(upwardbound == 0)]
    n_atrisk[i] <- length(atrisk)
    if(n_atrisk[i] <= cutoff) break

    ## Debugging
    ##if(i > 2) if(n_atrisk[i] < n_atrisk[i - 1]) break

    ## Calculate the early warning indicators here
    ## The nodes that feed into the early warning indicators are contained in `atrisk`
    samples <- (wl+1):stepT
    ## Covariance matrix --> dominant eigenvalue
    eig <- sampled_eigenmethod(X, A, samples = samples, nodes = atrisk)
    fulleig <- sampled_eigenmethod(X, A, samples = samples, nodes = V(g))
    ## Moran's I
    moranI <- sampled_MoranI(X, A, nodes = atrisk, t = stepT)
    fullI <- ape::Moran.I(x, A)$observed
    ## Individual-level
    ## only those in `atrisk` with the highest degree
    if(any(!(sentinels %in% atrisk))) {
        sentinels <- sentinels[which(sentinels %in% atrisk)]
        n_new <- n_sentinels - length(sentinels)
        available <- atrisk[which(!(atrisk %in% sentinels))]
        ##deg <- degree(g, v = atrisk)
        sentinels <- c(sentinels, choose_sentinels(g, n = n_new, v = available))
    }
    ssd <- sampled_sdmethod(X[, sentinels], t = stepT, wl = wl)

    ## Debugging
    ## print(D)
    ## print(x)
    ## print(length(atrisk))
    ## print(sentinels)

    ## if(which_ew == "sd") {
    ##     ew <- apply(X[(wl+1):stepT, ], 2, sd)
    ## } else if(which_ew == "ac") {
    ##     window <- (wl+1):stepT
    ##     ew <- apply(X, 2, function(x) cor(x[window], x[window-10]))
    ## }
    ## ew[-atrisk] <- NA
    ## I <- ape::Moran.I(X[stepT, ], A)$observed # calculations are made here, so need to write the functions above

    ## ews[[i]] <- ew
    ## Is[i] <- I
    eigs[i] <- eig
    fulleigs[i] <- fulleig
    moranIs[i] <- moranI
    fullIs[i] <- fullI
    ssds[[i]] <- ssd
    Ds[i] <- D
    Us[i] <- u[stressnode]

    if(which_param == "D") {
        D <- D + 0.001
    } else if(which_param == "u") {
        if(which_net == "me") {
            u[stressnode] <- u[stressnode] + 0.01
        } else if(which_net == "pa") {
            u[stressnode] <- u[stressnode] + 0.1
        }
    }
    
    i <- i + 1
}

ssds <- do.call(rbind, ssds)
##ews <- do.call(rbind, ews)
k <- degree(g)

## Initial test plots
## dev.new(width = 10, height = 5)
## plot(Ds, n_atrisk[-length(n_atrisk)]) ## showing the stepping down of "at risk" nodes

## eigs_withna <- eigs
## eigs_withna[which(eigs_withna > .01)] <- NA
## plot(Ds, eigs_withna) 
## showing that the covariance matrix --> eigenvalue method works, but that the magnitude of the relevant signal is quite small compared to the "transition is happening" signal
## try log?
dev.new(width = 10, height = 5)
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(Ds, log(eigs), type = "o", pch = 1, col = "steelblue",
     xlab = "D", ylab = expression(ln(lambda[1])), sub = "'At Risk' Nodes Only") ## log scale shows it beautifully
points(Ds, log(fulleigs), type = "o", pch = 3, col = "navy")
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], pch = 2, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("At Risk Nodes", "All Nodes", "# At Risk Nodes"), bty = "n",
       pch = c(1, 3, 2), col = c("steelblue", "navy", "darkorchid"), inset = c(-0.2, 0))

dev.new(width = 10, height = 5)
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(Ds, moranIs, type = "o", col = "red2",
     xlab = "D", ylab = "Moran's I", sub = "'At Risk' Nodes Only") ## showing that whatever I did above with Moran's I (which I should check) doesn't work very well after the initial step, but does up until that point
points(Ds, fullIs, type = "o", pch = 3, col = "red4")
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], pch = 2, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("At Risk Nodes", "All Nodes", "# At Risk Nodes"), bty = "n",
       pch = c(1, 3, 2), col = c("red2", "red4", "darkorchid"), inset = c(-0.2, 0))

dev.new(width = 10, height = 5)
par(mar = c(5, 4, 4, 8), xpd = TRUE)
matplot(Ds, log(ssds), type = "o", pch = 1, col = "olivedrab",
        xlab = "D", ylab = expression(ln(sd(x[i]))), sub = "Sentinel Nodes")
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], pch = 2, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("Sentinel Nodes", "# At Risk Nodes"), bty = "n",
       pch = c(1, 2), col = c("olivedrab", "darkorchid"), inset = c(-0.2, 0))


## ## First Plot
## filenamestem <- "./img/n-earlywarnings"
## filename <- paste0(filenamestem, "-", which_net, "-", which_ew, "-", which_param, ".pdf")

## if(save_plot) {
##     pdf(filename, height = 7, width = 21)
## } else dev.new(width = 21, height = 7)

## colorfun <- colorRamp(c("navy", "dodgerblue1"), space = "Lab")
## nodecolor <- colorfun(k/max(k))
## V(g)$color <- apply(nodecolor, 1, function(x) rgb(t(x), maxColorValue = 255))
## par(mai = c(1, .8, .8, .8) + .02)
## if(which_param == "D") {
##     matplot(Ds, ews, type = "p", pch = 1, lty = 1, col = V(g)$color,
##             xlab = "D", ylab = bquote(.(which_ew)(x[i])))
## } else if(which_param == "u") {
##     matplot(Us, ews, type = "p", pch = 1, lty = 1, col = V(g)$color,
##             xlab = "u", ylab = bquote(.(which_ew)(x[i])))
## }
## legend("left", pch = c(1, 1, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
##        col = c("dodgerblue1", "navy", "darkorchid"), bty = "n",
##        legend = c("Highest degree", "Lowest degree", "# Nodes 'At Risk'"))
## par(new = TRUE)
## if(which_param == "D") {
##     plot(Ds, n_atrisk[-length(n_atrisk)], type = "l", lwd = 2, lty = 1, col = "darkorchid",
##          axes = FALSE, xlab = "", ylab = "")
## } else if(which_param == "u") {
##     plot(Us, n_atrisk[-length(n_atrisk)], type = "l", lwd = 2, lty = 1, col = "darkorchid",
##          axes = FALSE, xlab = "", ylab = "")
## }
## axis(side = 4, at = pretty(range(n_atrisk[-length(n_atrisk)])))
## mtext("# Nodes", side = 4, line = 3)
## if(save_plot) dev.off()

## ## Second plot
## filenamestem <- "./img/n-earlywarnings"
## filename <- paste0(filenamestem, "-", which_net, "-I-", which_param, ".pdf")
## ## if(which_net == "me") {
## ##     filename <- paste0(filenamestem, which_ew, "-"
## ## } else if(which_net == "pa") {
## ##     filename <- "./img/n-earlywarnings-I-PA.pdf"
## ## }

## if(save_plot) {
##     pdf(filename, height = 7, width = 21)
## } else dev.new(width = 21, height = 7)

## par(mai = c(1, .8, .8, .8) + .02)
## if(which_param == "D") {
##     plot(Ds, Is, type = "p", lwd = 2, col = "red3",
##          xlab = "D", ylab = "Moran's I")
## } else if(which_param == "u") {
##     plot(Us, Is, type = "p", lwd = 2, col = "red3",
##          xlab = "u", ylab = "Moran's I")
## }
## legend("left", pch = c(1, NA), lty = c(NA, 1), lwd = c(2, 2),
##        col = c("red3", "darkorchid"), bty = "n",
##        legend = c("Moran's I", "# Nodes 'At Risk'"))
## par(new = TRUE)
## if(which_param == "D") {
##     plot(Ds, n_atrisk[-length(n_atrisk)], type = "l", lwd = 2, lty = 1, col = "darkorchid",
##          axes = FALSE, xlab = "", ylab = "")
## } else if(which_param == "u") {
##     plot(Us, n_atrisk[-length(n_atrisk)], type = "l", lwd = 2, lty = 1, col = "darkorchid",
##          axes = FALSE, xlab = "", ylab = "")
## }
## axis(side = 4, at = pretty(range(n_atrisk[-length(n_atrisk)])))
## mtext("# Nodes", side = 4, line = 3)
## if(save_plot) dev.off()
## ##print(n_atrisk)
