## Code by Neil MacLaren 2/7/2022

library(igraph)
library(doublewells)

set.seed(123)

which_net <- "pa" # "me"
which_param <- "D" # "u"
save_plots <- TRUE # FALSE

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

## Simulation Parameters
cutoff <- .25*nnodes
atrisk <- V(g)
stepT <- 500
wl <- 250#.5*stepT
n_sentinels <- 5
sentinels <- choose_sentinels(g, n = n_sentinels)
Ds <- numeric()
Us <- numeric()
eigs <- numeric()
fulleigs <- numeric()
moranIs <- numeric()
fullIs <- numeric()
ssds <- list()
n_atrisk <- numeric()

## Main Loop
i <- 1
##for(i in 1:20) {
while(length(atrisk) > cutoff) {
    X <- matrix(0, nrow = stepT, ncol = nnodes)

    for(t in 1:stepT) {
        X[t, ] <- x
        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
    }

    ## Analysis
    atrisk <- V(g)[which(lowerstate(x) == 1)] # discretizes x_i on purpose: we may not know what the numeric state of a node is, yet still be able to tell whether it is in the lower (healthy) or upper (diseased) state.
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
    ## only those in `atrisk` with the highest degree, but we prefer the neighbors of transitioned nodes
    if(any(!(sentinels %in% atrisk))) {
        transitioned <- sentinels[which(!(sentinels %in% atrisk))]
        sentinels <- sentinels[which(sentinels %in% atrisk)]
        for(v in transitioned) {
            nbs <- neighbors(g, v)
            avail <- nbs[which(lowerstate(x[nbs]) == 1)]
            if(length(avail) > 0) {
                ## select from among those with highest degree
                sentinels <- c(sentinels, choose_sentinels(g, n = 1, v = avail))
            } else next
        }
        n_needed <- n_sentinels - length(sentinels)
        if(n_needed > 0) {
            ## select randomly from among the remaining lowerstate nodes
            avail <- atrisk[which(!(atrisk %in% sentinels))]
            sentinels <- c(sentinels, choose_sentinels(g, n = n_needed, v = avail))
        }
        ## sentinels <- sentinels[which(sentinels %in% atrisk)]
        ## n_new <- n_sentinels - length(sentinels)
        ## available <- atrisk[which(!(atrisk %in% sentinels))]
        ## sentinels <- c(sentinels, choose_sentinels(g, n = n_new, v = available))
    }
    ssd <- sampled_sdmethod(X[, sentinels], t = stepT, wl = wl)

    ## Debugging
    ## print(D)
    ## print(x)
    ## print(length(atrisk))
    ## print(sentinels)

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

## Collating
ssds <- do.call(rbind, ssds)

## Plotting
if(save_plots) {
    pdf("./img/ew-covmat-eig.pdf", width = 10, height = 5)
} else {
    dev.new(width = 10, height = 5)
}
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(Ds, log(eigs), type = "o", pch = 1, col = "steelblue",
     xlab = "D", ylab = expression(ln(lambda[1])))
points(Ds, log(fulleigs), type = "p", pch = 3, col = "navy")
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], pch = 2, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("At Risk Nodes", "All Nodes", "# At Risk Nodes"), bty = "n",
       pch = c(1, 3, 2), col = c("steelblue", "navy", "darkorchid"), inset = c(-0.2, 0))
if(save_plots) dev.off()

if(save_plots) {
    pdf("./img/ew-moranI.pdf", width = 10, height = 5)
} else {
    dev.new(width = 10, height = 5)
}
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(Ds, moranIs, type = "o", col = "red2",
     xlab = "D", ylab = "Moran's I")
points(Ds, fullIs, type = "o", pch = 3, col = "red4")
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], pch = 2, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("At Risk Nodes", "All Nodes", "# At Risk Nodes"), bty = "n",
       pch = c(1, 3, 2), col = c("red2", "red4", "darkorchid"), inset = c(-0.2, 0))
if(save_plots) dev.off()

if(save_plots) {
    pdf("./img/ew-sentinelnodes.pdf", width = 10, height = 5)
} else {
    dev.new(width = 10, height = 5)
}
par(mar = c(5, 4, 4, 8), xpd = TRUE)
matplot(Ds, log(ssds), type = "o", pch = 1, col = "olivedrab",
        xlab = "D", ylab = expression(ln(sd(x[i]))))
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], pch = 2, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
legend("bottomright", legend = c("Sentinel Nodes", "# At Risk Nodes"), bty = "n",
       pch = c(1, 2), col = c("olivedrab", "darkorchid"), inset = c(-0.2, 0))
if(save_plots) dev.off()
