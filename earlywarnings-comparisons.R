## Code by Neil MacLaren 2/15/2022

library(igraph)
library(doublewells)
MoranI <- ape::Moran.I

##set.seed(123)

save_plots <- FALSE

## Graph Selection
generate_new <- FALSE
use_empirical <- FALSE
which_network <- "mine" # "jpr" "surfers"
if(generate_new) {
    nnodes <- 100
    outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
    g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
    write_graph(g, "./data/testgraph-BA.gml", format = "gml")
} else if(!use_empirical) {
    g <- read_graph("./data/testgraph-BA.gml", format = "gml")
    nnodes <- vcount(g)
} else if(which_network == "mine") {
    data("mine", package = "networkdata")
    g <- mine
} else if(which_network == "jpr") {
    data("jpr", package = "networkdata")
    g <- jpr
} else if(which_network == "surfers") {
    data("surfersb", package = "networkdata")
    g <- surfersb
}

    
A <- as_adj(g, type = "both", sparse = FALSE)

## Node Systems
r <- c(1, 4, 7) # double well parameters
s <- 0.005 # sd of noise process
D <- 0.21 # connection strength
p <- if(s > 0) 3*s else 0.015 # perturbation strength
u <- rep(0, nnodes) # stress vector

dt <- 0.01

initialx <- rep(1, nnodes) + noise(nnodes, s)
x <- initialx

## Simulation Parameters
stepsize <- 1e-3 
cutoff <- .25*nnodes
in_lowerstate <- V(g)
stepT <- 5000
wl <- 250
samples <- (stepT - wl + 1):stepT
n_lowerstate <- numeric()
n_sentinels <- 5
Ds <- numeric()

## Storage vectors for all early warning indicators
maxeig <- list(all = numeric(), lower = numeric(), sentinel = numeric())
maxsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
avgsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
moranI <- list(all = numeric(), lower = numeric(), sentinel = numeric())
maxac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
avgac <- list(all = numeric(), lower = numeric(), sentinel = numeric())

## Main Loop
i <- 1
while(length(in_lowerstate) > cutoff) {
    X <- matrix(0, nrow = stepT, ncol = nnodes)

    for(t in 1:stepT) {
        X[t, ] <- x
        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
    }

    ## Exit condition
    in_lowerstate <- V(g)[which(lowerstate(x) == 1)]
    if(length(in_lowerstate) <= cutoff) break

    ## Else, store the number of "at risk" nodes and continue
    n_lowerstate[i] <- length(in_lowerstate)

    ## Determine sentinels
    sentinels <- sentinel_ranking_ts(g, X, t = stepT, wl = wl, n = n_sentinels)

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
    moranI$lower[i] <- MoranI(x[in_lowerstate], A[in_lowerstate, in_lowerstate])$observed
    if(ecount(induced_subgraph(g, sentinels)) > 0) {
        moranI$sentinel[i] <- MoranI(x[sentinels], A[sentinels, sentinels])$observed
    } else moranI$sentinel[i] <- NA

    ac <- list()
    ac$all <- sampled_acmethod(X, stepT, wl, lag = 1)
    ac$lower <- sampled_acmethod(X[, in_lowerstate], stepT, wl, lag = 1)
    ac$sentinel <- sampled_acmethod(X[, sentinels], stepT, wl, lag = 1)

    for(j in 1:length(ac)) maxac[[j]][i] <- max(ac[[j]])
    for(j in 1:length(ac)) avgac[[j]][i] <- mean(ac[[j]])

    Ds[i] <- D

    ## Iterate
    D <- D + stepsize
    i <- i + 1
}

## May need to adjust column names here
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
kendalls <- kendalls[kendalls$n_steps > 15, ]
                                        #floor(quantile(kendalls$n_steps, probs = .75)), ]
                                        #mean(kendalls$n_steps)), ]

## Plot results
kendalls <- kendalls[order(kendalls$n_steps), ]

if(save_plots) {
    pdf("./img/rough-comparison.pdf", width = 8, height = 7)
} else dev.new(width = 20, height = 4)
##figcols <- columns[-grep("moran", columns)]
##palette(rainbow(length(figcols)))
ewis <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
par(
    mar = c(5, 4, 4, 4),#8),
    xpd = TRUE,
    mfrow = c(1, length(ewis))
)
for(i in 1:length(ewis)) {
    figcols <- columns[grep(ewis[i], columns)]
    palette(c("red2", "dodgerblue", "darkorchid", "darkorange", "gray"))
    matplot(kendalls$n_steps, kendalls[, figcols], type = "p",
            pch = 1, col = 1:length(figcols), cex = 1.5,
            lwd = 1.5, lty = 1,
            xlab = "# Steps in Steady State", ylab = expression(tau),
            xlim = rev(range(kendalls$n_steps)), ylim = c(-.3, 1))
    legend("topright", bty = "n",# inset = c(-0.29, 0),
           legend = figcols, col = 1:length(figcols), pch = 1, pt.cex = 1.5, pt.lwd = 1.5)
}
if(save_plots) dev.off()

as.matrix(round(colMeans(kendalls[, columns], na.rm = TRUE), 3))

### Old code below here ### 

## agg <- data.frame(n_lowerstate = unique(df$n_lowerstate))






## for(col in columns) {
##     ag <- aggregate(
##         df[, col], by = list(n_lowerstate = df[, "n_lowerstate"]),
##         function(x) {
##             cor(df[, "Ds"], x, method = "kendall", use = "pairwise.complete.obs"))
##     colnames(ag)[2] <- col
##     agg <- merge(agg, ag, by = "n_lowerstate")
## }

    

## apply(df[, columns], 2, function(x) cor(Ds

## for(i in 1:length(columns)) {
    

## ## Finish this after confirm column names above
## agg <- df %>% group_by(n_lowerstate) %>%
##     summarize(
##         D = mean(Ds),
##         n = n_distinct(Ds)
##     )
        
