## Code by Neil MacLaren 2/25/2022

library(igraph)
library(doublewells)
    
##set.seed(123)

save_plots <- FALSE # TRUE

## Make an organized file that runs the analysis on a specified set of networks.
## Looking for networks with multistage transitions.
## The networks should be:
## - Scale free with several exponents (2 & 3, maybe 5)
## - LFR
## - [ ] dolphins
## - [ ] jpr
## - pira

choices <- c("dolphins", "pira")
data(list = choices)

choice <- "dolphins" # "pira"

g <- get(choice)

nnodes <- vcount(g)
A <- as_adj(g, type = "both", sparse = FALSE)

### Node Systems
r <- c(1, 4, 7) # double well parameters
s <- 0.01 # sd of noise process; 0.005 for production
D <- 0.5#0.01
p <- if(s > 0) 3*s else 0.015 # perturbation strength
u <- rep(0, nnodes) # stress vector

τ <- 75 # duration in time units
dt <- 0.01 # step size for integration
T <- τ/dt
s <- s*sqrt(dt) # reassign s to scale for dt

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

stepsize <- .01# 5e-3 # for incrementing D
cutoff <- .1*nnodes # to stop simulation; production value is .1
in_lowerstate <- V(g) # initial vector of nodes in the lower state (all of them)
n_lowerstate <- numeric() # storage vector for the number of nodes in the lower state
Ds <- numeric() # storage vector for the state of the bifurcation parameter

## Storage vectors for all early warning indicators
maxeig <- list(all = numeric(), lower = numeric(), sentinel = numeric())
maxsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
avgsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
maxac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
avgac <- list(all = numeric(), lower = numeric(), sentinel = numeric())

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

    ac <- list()
    ac$all <- sampled_acmethod(X, samples, lag = lag)
    ac$lower <- sampled_acmethod(X[, in_lowerstate], samples, lag = lag)
    ac$sentinel <- sampled_acmethod(X[, sentinels], samples, lag = lag)

    for(j in 1:length(ac)) maxac[[j]][i] <- max(ac[[j]])
    for(j in 1:length(ac)) avgac[[j]][i] <- mean(ac[[j]])

    Ds[i] <- D

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
maxac <- do.call(cbind, maxac)
colnames(maxac) <- paste("maxac", colnames(maxac), sep = "_")
avgac <- do.call(cbind, avgac)
colnames(avgac) <- paste("avgac", colnames(avgac), sep = "_")

## Analysis
df <- data.frame(
    n_lowerstate = n_lowerstate, # [-length(n_lowerstate)]
    Ds = Ds
)
df <- cbind(df, maxeig, maxsd, avgsd, maxac, avgac)
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

## Plot results
kendalls <- kendalls[order(kendalls$n_steps), ]
as.matrix(round(colMeans(kendalls[, columns], na.rm = TRUE), 3))
vcount(g)
edge_density(g)

wd <- 7
ht <- 7
ewi <- "avgac"
columns <- colnames(df)[grep(ewi, colnames(df))]
colors <- c("sienna", "navy", "olivedrab")
if(save_plots) {
    pdf(paste0("./img/", choice, "-earlywarnings.pdf"), width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mar = c(5, 4, 4, 4)+.1, xpd = TRUE) #13
plot(df$Ds, df$n_lowerstate,
     ##type = "l", lwd = 2, lty = 1,
     type = "p", pch = 1,
     col = "darkorchid",
     xlim = range(df$Ds), ylim = c(0, vcount(g)),
     xlab = "D", ylab = "# Nodes in Lower State")
## legend("bottomleft", bty = "n", legend = "# Nodes in Lower State",
##        ##lty = 1, lwd = 2,
##        pch = 1,
##        col = "darkorchid") #pt.cex = 1.5, pch = 19
legend("bottomleft", ncol = 4, bty = "y",
       pch = 1:4, col = c("darkorchid", colors),
       legend = c("State", "All", "Lower", "Sentinel")
                  ##expression(X[all]), expression(X[lower]), expression(X[sentinel]))
       )
par(new = TRUE)
## if(stepsize == 1e-2) {
##     plotsamples <- 1:nrow(df)
## } else if(stepsize == 5e-3) {
##     plotsamples <- seq(1, nrow(df), by = 5)
## } else plotsamples <- seq(1, nrow(df), by = 10)
plotsamples <- 1:nrow(df)
matplot(df$Ds[plotsamples], df[plotsamples, columns], type = "p",
        pch = 2:4, lwd = 1, lty = 2, col = colors,
        xlim = range(df$Ds), ylim = c(0, 1), axes = FALSE,
        xlab = "", ylab = "")
## legend("topright", bty = "n",
##        legend = c(expression(X[all]), expression(X[lower]), expression(X[sentinel])),
##        pch = 2:4, col = colors)
axis(4, at = pretty(c(0, 1)))
mtext(#expression(paste("Average ", "Autocorrelation ", (X[i]))),
    "Average Autocorrelation",
    side = 4, cex = 1, font = 1, line = 3)
if(save_plots) dev.off()
