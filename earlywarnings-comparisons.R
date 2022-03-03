## Code by Neil MacLaren 2/15/2022

library(igraph)
library(doublewells)
MoranI <- ape::Moran.I
    
set.seed(123)

save_plots <- FALSE

### Graph Selection
## pref_attach, max_entropy, sphere_surface, me_islands, random_regular, small_world
## mine, jpr, pira, sn_auth;;; mine went to D = 1.044?!
## computational cost of "petster" ~1800 nodes in petster GCC is too much
choice <- "dolphins"
## if(choice %in% c("dolphins", "netscience-lcc", "lfr-n100k6-decayrate")) {
##     filename <- paste0("./data/", choice, ".mat")
##     edgelist <- as.matrix(read.table(file = filename, header = FALSE, sep = " ", skip = 1)[, 1:2])
##     g <- graph_from_edgelist(edgelist, directed = FALSE)
## } else {
##     data(list = choice)
##     g <- get(choice)
## }
data(list = choice)
g <- get(choice)

nnodes <- vcount(g)
A <- as_adj(g, type = "both", sparse = FALSE)

### Node Systems
r <- c(1, 4, 7) # double well parameters
s <- 0.005 # sd of noise process; 0.005 for production
## if(choice %in% c("pref_attach", "jpr", "lfr-n100k6-decayrate")) {# connection strength
##     D <- 0.24 # production value is 0.20
## } else if(choice %in% c("small_world", "random_regular")) {
##     D <- 0.7
## } else D <- 0.5
D <- 0.2
p <- if(s > 0) 3*s else 0.015 # perturbation strength
u <- rep(0, nnodes) # stress vector

τ <- 75 # duration in time units
dt <- 0.01 # step size for integration
T <- τ/dt

initialx <- rep(1, nnodes) + noise(nnodes, s) # initial values of x_i

nsamples <- 250 # number of samples for early warning indicators; 250 for production
sample_spacing <- .1 # in τ scale
sample_spacing <- sample_spacing/dt # in Δt scale
samples <- seq(T, T - (nsamples*sample_spacing), by = -sample_spacing) # vector of indices to sample

lag <- .1 # lag for autocorrelation, τ scale; .1 for production; .2 is better, and .3 and .5
lag <- lag/dt # lag for autocorrelation, Δt scale

n_sentinels <- 5 # set the number of sentinel nodes

### Simulation Parameters
x <- initialx # an updating vector of x_i

stepsize <- 1e-3 # for incrementing D
cutoff <- .1*nnodes # to stop simulation; production value is .1
in_lowerstate <- V(g) # initial vector of nodes in the lower state (all of them)
##stepT <- 5000
##wl <- 250
##samples <- (stepT - wl + 1):stepT
n_lowerstate <- numeric() # storage vector for the number of nodes in the lower state
Ds <- numeric() # storage vector for the state of the bifurcation parameter

## Storage vectors for all early warning indicators
maxeig <- list(all = numeric(), lower = numeric(), sentinel = numeric())
maxsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
avgsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
moranI <- list(all = numeric(), lower = numeric(), sentinel = numeric())
maxac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
avgac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
sentinel_history <- list()

### Main Loop
i <- 1
while(length(in_lowerstate) > cutoff) {
    X <- matrix(0, nrow = T, ncol = nnodes) # storage matrix for the integration step states of x

    ## Integration
    for(t in 1:T) {
        X[t, ] <- x
        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
    }

    ## Exit condition
    in_lowerstate <- V(g)[which(lowerstate(x) == 1)]
    if(length(in_lowerstate) <= cutoff) break
    if(length(in_lowerstate) < n_sentinels) break

    ## Else, store the number of "at risk" nodes and continue
    n_lowerstate[i] <- length(in_lowerstate)

    ## Determine sentinels
    sentinels <- sentinel_ranking_ts(g, X, samples, n = n_sentinels)

    X.all <- X[samples, ]
    X.lower <- X[samples, in_lowerstate]
    X.sentinels <- X[samples, sentinels]

    ## Calculate and store
    maxeig$all[i] <- sampled_eigenmethod(X, samples = samples, nodes = V(g))
    maxeig$lower[i] <- sampled_eigenmethod(X, samples = samples, nodes = in_lowerstate)
    maxeig$sentinel[i] <- sampled_eigenmethod(X, samples = samples, nodes = sentinels)

    maxsd$all[i] <- max(apply(X.all, 2, sd))
    maxsd$lower[i] <- max(apply(X.lower, 2, sd))
    maxsd$sentinel[i] <- max(apply(X.sentinels, 2, sd))

    avgsd$all[i] <- mean(apply(X.all, 2, sd))
    avgsd$lower[i] <- mean(apply(X.lower, 2, sd))
    avgsd$sentinel[i] <- mean(apply(X.sentinels, 2, sd))

    moranI$all[i] <- MoranI(x, A)$observed
    if(is_connected(induced_subgraph(g, V(g)[in_lowerstate]))) {
        moranI$lower[i] <- MoranI(x[in_lowerstate], A[in_lowerstate, in_lowerstate])$observed
    } else moranI$lower[i] <- NA
    if(ecount(induced_subgraph(g, sentinels)) > 0) {
        moranI$sentinel[i] <- MoranI(x[sentinels], A[sentinels, sentinels])$observed
    } else moranI$sentinel[i] <- NA

    ac <- list()
    ac$all <- sampled_acmethod(X, samples, lag = lag)
    ac$lower <- sampled_acmethod(X[, in_lowerstate], samples, lag = lag)
    ac$sentinel <- sampled_acmethod(X[, sentinels], samples, lag = lag)

    for(j in 1:length(ac)) maxac[[j]][i] <- max(ac[[j]])
    for(j in 1:length(ac)) avgac[[j]][i] <- mean(ac[[j]])

    Ds[i] <- D
    sentinel_history[[i]] <- sentinels

    ## Iterate
    D <- D + stepsize
    i <- i + 1
}

## Convert lists to data frames
maxeig <- do.call(cbind, maxeig)
colnames(maxeig) <- paste("maxeig", colnames(maxeig), sep = "_")
maxsd <- do.call(cbind, maxsd)
colnames(maxsd) <- paste("maxsd", colnames(maxsd), sep = "_")
avgsd <- do.call(cbind, avgsd)
colnames(avgsd) <- paste("avgsd", colnames(avgsd), sep = "_")
moranI <- do.call(cbind, moranI)
colnames(moranI) <- paste("moranI", colnames(moranI), sep = "_")
maxac <- do.call(cbind, maxac)
colnames(maxac) <- paste("maxac", colnames(maxac), sep = "_")
avgac <- do.call(cbind, avgac)
colnames(avgac) <- paste("avgac", colnames(avgac), sep = "_")
sentinel_history <- do.call(rbind, sentinel_history)

## Analysis
## Goal is to compare Kendall's τ across each early warning indicator stored above.
## Rank correlation is within steady states (i.e., same number of nodes in the lower state) between D (always increasing) and the early warning indicator.
df <- data.frame(
    n_lowerstate = n_lowerstate, # [-length(n_lowerstate)]
    Ds = Ds
)
df <- cbind(df, maxeig, maxsd, avgsd, moranI, maxac, avgac)
columns <- colnames(df)[3:length(colnames(df))]

dfsplit <- split(df, factor(df$n_lowerstate))
kendalls <- lapply(dfsplit, function(x) {
    cor(x[, c("Ds", columns)], method = "kendall", use = "pairwise.complete.obs")[1, -1]
})
kendalls <- as.data.frame(do.call(rbind, kendalls))
kendalls$n_lowerstate <- as.integer(rownames(kendalls))
kendalls$n_steps <- sapply(dfsplit, nrow)

## discard minor transitions
kendalls <- kendalls[kendalls$n_steps > 15, ] # 15 for production 
                                        #floor(quantile(kendalls$n_steps, probs = .75)), ]
                                        #mean(kendalls$n_steps)), ]

## Plot results
kendalls <- kendalls[order(kendalls$n_steps), ]
as.matrix(round(colMeans(kendalls[, columns], na.rm = TRUE), 3))
vcount(g)
edge_density(g)

ht <- 7
wd <- 14
if(save_plots) {
    pdf(paste0("./img/", choice, "-rough-comparison.pdf"), width = wd, height = ht)
} else dev.new(width = wd, height = ht)
##figcols <- columns[-grep("moran", columns)]
##palette(rainbow(length(figcols)))
ewis <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
par(
    mar = c(5, 4, 4, 4),#8),
    xpd = TRUE,
    mfrow = c(2, 3)#length(ewis))
)
for(i in 1:length(ewis)) {
    figcols <- columns[grep(ewis[i], columns)]
    palette(c("red2", "dodgerblue", "darkorchid", "darkorange", "gray"))
    matplot(kendalls$n_steps, kendalls[, figcols], type = "p",
            pch = 1, col = 1:length(figcols), cex = 1.5,
            lwd = 1.5, lty = 1,
            xlab = "# Steps in Steady State", ylab = expression(tau),
            xlim = rev(range(kendalls$n_steps)), ylim = c(0, 1))
    legend("bottomleft", bty = "n",# inset = c(-0.29, 0),
           legend = figcols, col = 1:length(figcols), pch = 1, pt.cex = 1.5, pt.lwd = 1.5)
}
if(save_plots) dev.off()

ht <- 5
wd <- 10
if(save_plots) {
    pdf(paste0("./img/", choice, "-overview.pdf"), height = ht, width = wd)
} else dev.new(width = wd, height = ht)
par(mfrow = c(1, 2))
plot(g, vertex.label = "", vertex.size = 3)
hist(degree(g), main = "", xlab = "Degree", ylab = "Frequency", breaks = 20)
if(save_plots) dev.off()

wd <- 15
ht <- 5
if(save_plots) {
    pdf(paste0("./img/", choice, "-earlywarnings.pdf"), width = wd, height = ht)
} else dev.new(width = wd, height = 5)
par(mar = c(5, 4, 4, 13), xpd = TRUE)
## some stuff to plot the early warning indicator(s) over time
ewi <- "avgac"
columns <- colnames(df)[grep(ewi, colnames(df))]
colors <- c("sienna", "navy", "olivedrab")
matplot(df$Ds, df[, columns], type = "p", pch = c(1, 3, 4), lwd = 1, col = colors,
        xlab = "D", ylab = "Early Warning Indicator Value")
legend("bottomright", bty = "n", inset = c(-0.21, 0),
       legend = c(columns, "# Nodes"),
       pch = c(1, 3, 4, 2), col = c(colors, "darkorchid"))
par(new = TRUE)
plot(df$Ds, df$n_lowerstate, pch = 2, lwd = 2, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
axis(4, at = pretty(range(df$n_lowerstate)))
mtext("# Nodes in Lower State", side = 4, cex = 1, font = 1, line = 3)
if(save_plots) dev.off()
