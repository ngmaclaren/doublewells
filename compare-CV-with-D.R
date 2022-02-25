## Code by Neil MacLaren, 2/24/2022

library(igraph)
library(doublewells)

save_plots <- FALSE # TRUE
use_noise <- FALSE # TRUE
run_simulation <- FALSE # TRUE

## choices <- c(
##     "max_entropy", "me_islands", "pref_attach", "LFR",
##     "scale_free_2", "scale_free_2.1", "scale_free_2.2", "scale_free_2.4", "scale_free_2.6",
##     "scale_free_2.8", "scale_free_3", "scale_free_5", "scale_free_10", "scale_free_100",
##     "dolphins", "jpr", "pira", "sn_auth", "NNS"
## )
choices <- empiricals # list of empirical networks from doublewells package
data(list = choices)

if(run_simulation) {
    results <- list()
    for(i in 1:length(choices)) {
        choice <- choices[i]
        print(choice)
        g <- get(choice)
        A <- as_adj(g, type = "both", sparse = FALSE) # this eliminates edge weights -- only returns adjacency
        
        nnodes <- vcount(g)
        nedges <- ecount(g)
        k <- degree(g)
        CV <- sd(k)/mean(k)

        ## Node Systems
        r <- c(1, 4, 7) # double well parameters
        if(use_noise) s <- 0.005 else s <- 0 # sd of noise process
        D <- 0.15 # this is the initial value of D
        p <- if(s > 0) 3*s else 0.015 # perturbation strength
        u <- rep(0, nnodes) # stress vector

        τ <- 75 # duration in time units
        dt <- 0.01 # step size for integration
        T <- τ/dt

        initialx <- rep(1, nnodes) + noise(nnodes, s) # initial values of x_i

        ## Simulation Parameters
        x <- initialx # an updating vector of x_i

        stepsize <- 1e-2 # for incrementing D
        cutoff <- .1*nnodes # to stop simulation; production value is .1
        in_lowerstate <- V(g) # initial vector of nodes in the lower state (all of them)

        ## Storage Vectors
        n_lowerstate <- numeric() # storage vector for the number of nodes in the lower state
        Ds <- numeric() # storage vector for the state of the bifurcation parameter

        ## Integrate and Store
        j <- 1
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

            ## Store
            n_lowerstate[j] <- length(in_lowerstate)
            Ds[j] <- D

            ## Iterate
            D <- D + stepsize
            j <- j + 1
        }

        ## Make Results Data Frame
        df <- data.frame(graph = choice, CV = CV, D = Ds, n_lowerstate = n_lowerstate)
        results[[i]] <- df
    }
    results <- do.call(rbind, results)
    save(results, file = "./data/results.rda")
} else {
    load("./data/results.rda")
    stepsize <- results$D[2] - results$D[1]
}
    

parts <- split(results, factor(results$graph))

find_range <- function(df, stepsize = 1e-3) {
    nnodes <- df$n_lowerstate[1] # the initial D is below that at which any of these networks starts to show transitions
    firststate <- floor(nnodes - (nnodes*.1))
    first <- which(df$n_lowerstate < firststate)[1]
    if(is.na(first)) return(0)
    ## otherwise...
    last <- nrow(df)
    firstD <- df$D[first]
    lastD <- df$D[last] + stepsize
    lastD - firstD
}

dat <- lapply(parts, function(x) {
    c(rangeD = find_range(x), CV = unique(x$CV), nstages = length(unique(x$n_lowerstate)))
})

df <- do.call(rbind, dat)
##empirical_rows <- c("dolphins", "jpr", "NNS", "pira", "sn_auth")

ht <- 7
wd <- 14
if(save_plots) {
    pdf("./img/compare-CV-with-D.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mfrow = c(1, 2))
plot(nstages ~ CV, data = df, pch = 19, col = "firebrick",
     xlab = "CV of Degree Distribution", ylab = "Number of Stages in Double-Wells Transition")
plot(rangeD ~ CV, data = df, pch = 19, col = "steelblue",
     xlab = "CV of Degree Distribution", ylab = "Magnitude of Range of D Over Which Transitions Occur")
if(save_plots) dev.off()
