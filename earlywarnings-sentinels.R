## Code by Neil MacLaren 2/10/2022

library(igraph)
library(doublewells)

set.seed(123)

save_plots <- FALSE

generate_new <- FALSE
if(generate_new) {
    nnodes <- 100
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
### the below may need to be fixed to support varying u
D <- 0.21 # connection strength
p <- if(s > 0) 3*s else 0.015 # perturbation strength
u <- rep(0, nnodes) # stress vector

dt <- 0.01

initialx <- rep(1, nnodes) + noise(nnodes, s)
##X <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx

## Simulation Parameters
stepsize <- 1e-4
cutoff <- .25*nnodes
atrisk <- V(g)
stepT <- 500
wl <- 250
samples <- (stepT - wl + 1):stepT
n_sentinels <- 5
Ds <- numeric()
eigs <- numeric() # only the eigs of the sentinel nodes
n_atrisk <- numeric()

## Main Loop
i <- 1
while(length(atrisk) > cutoff) {
    X <- matrix(0, nrow = stepT, ncol = nnodes)

    for(t in 1:stepT) {
        X[t, ] <- x
        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
    }

    ## Exit condition
    atrisk <- V(g)[which(lowerstate(x) == 1)]
    if(length(atrisk) <= cutoff) break

    ## Else, store the number of "at risk" nodes and continue
    n_atrisk[i] <- length(atrisk)

    ## Calculate EWI
    ## sentinels <- sentinel_ranking(g, x, n = n_sentinels)
    sentinels <- sentinel_ranking_ts(g, X, t = stepT, wl = wl, n = n_sentinels)
    eig <- sampled_eigenmethod(X, samples = samples, nodes = sentinels)

    ## Store
    eigs[i] <- eig
    Ds[i] <- D

    ## Iterate
    D <- D + stepsize
    i <- i + 1
}

## Analysis
## This will be a bit more complicated than before.
## Let's plot first to make sure everything's working.
if(save_plots) {
    pdf("./img/ew-newsentinelnodes.pdf", width = 15, height = 5)
} else {
    dev.new(width = 15, height = 5)
}
par(mar = c(5, 4, 4, 13), xpd = TRUE)
plot(Ds, log(eigs), type = "o", pch = 1, col = "sandybrown",
        xlab = "D", ylab = expression(ln(lambda[1])))
legend("bottomright", bty = "n", inset = c(-0.21, 0),
       legend = c("Sentinel Node EWI", "# Nodes in Lower State"),
       col = c("sandybrown", "darkorchid"), pch = c(1, 2))
par(new = TRUE)
plot(Ds, n_atrisk, pch = 2, col = "darkorchid",#[-length(n_atrisk)]
     axes = FALSE, xlab = "", ylab = "")
axis(4, at = pretty(range(n_atrisk)))
if(save_plots) dev.off()


df <- data.frame(Ds, eigs, n_atrisk)

##library(dplyr)
ag <- df %>% dplyr::group_by(n_atrisk) %>%
    dplyr::summarize(
               D = mean(Ds),
               n = dplyr::n_distinct(eigs),
               tau = cor(Ds, eigs, method = "kendall"))
ag %>% dplyr::arrange(-n) %>% print(n = 15)
plot(tau ~ n, data = subset(ag, n > median(n)))#, xlim = rev(range(ag$n_atrisk)), ylim = range(ag$tau, na.rm = TRUE))
plot(tau ~ D, data = subset(ag, n > median(n)))

## for(i in 100:26) {
##     with(subset(df, n_atrisk == i),
##          print(cor(Ds, eigs, method = "kendall")))
## }


## aggregate(cbind(Ds, eigs) ~ n_atrisk, data = df, FUN = cor, method = "kendall")

## with(subset(df, n_atrisk == 100),
##      cor(Ds, eigs, method = "kendall"))
