## Code by Neil MacLaren, 2/24/2022
## Parameters and bifurcation parameter for Prosenjit Kundu, 3/4/2022

library(igraph)
library(doublewells)

save_plots <- TRUE # FALSE
use_noise <- FALSE # TRUE
run_simulation <- FALSE # TRUE

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
        r <- c(1, 2, 5) # double well parameters
        if(use_noise) s <- 0.005 else s <- 0 # sd of noise process
        D <- 0.001 # this is the initial value of D
        p <- if(s > 0) 3*s else 0.015 # perturbation strength
        u <- .3 # this is the initial value
        U <- rep(u, nnodes) # stress vector

        τ <- 75 # duration in time units
        dt <- 0.01 # step size for integration
        T <- τ/dt

        initialx <- rep(.1, nnodes) + noise(nnodes, s) # initial values of x_i

        ## Simulation Parameters
        x <- initialx # an updating vector of x_i

        stepsize <- 1e-2 # for incrementing u
        cutoff <- .1*nnodes # to stop simulation; production value is .1
        in_lowerstate <- V(g) # initial vector of nodes in the lower state (all of them)

        ## Storage Vectors
        n_lowerstate <- numeric() # storage vector for the number of nodes in the lower state
        us <- numeric() # storage vector for the state of the bifurcation parameter

        ## Integrate and Store
        j <- 1
        while(length(in_lowerstate) > cutoff) {
            X <- matrix(0, nrow = T, ncol = nnodes) # storage matrix for the integration step states of x
            
            ## Integration
            for(t in 1:T) {
                X[t, ] <- x
                x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), U)
            }
            
            ## Exit condition
            in_lowerstate <- V(g)[which(lowerstate(x) == 1)]
            if(length(in_lowerstate) <= cutoff) break

            ## Store
            n_lowerstate[j] <- length(in_lowerstate)
            us[j] <- u

            ## Iterate
            u <- u + stepsize
            U <- rep(u, nnodes) # stress vector
            j <- j + 1
        }

        ## Make Results Data Frame
        df <- data.frame(graph = choice, vcount = vcount(g),
                         CV = CV, u = us, n_lowerstate = n_lowerstate)
        results[[i]] <- df
    }
    results <- do.call(rbind, results)
    save(results, file = "./data/cv-results.rda")
} else {
    load("./data/cv-results.rda")
    stepsize <- results$u[2] - results$u[1]
}
    

parts <- split(results, factor(results$graph))

find_range <- function(df, stepsize = 1e-3) {
    nnodes <- unique(df$vcount)
    firststate <- floor(nnodes - (nnodes*.1))
    first <- which(df$n_lowerstate < firststate)[1]
    if(is.na(first)) return(0)
    ## otherwise...
    last <- nrow(df)
    first_u <- df$u[first]
    last_u <- df$u[last] + stepsize
    last_u - first_u
}

dat <- lapply(parts, function(x) {
    c(range_u = find_range(x), CV = unique(x$CV), nstages = length(unique(x$n_lowerstate)))
})

df <- do.call(rbind, dat)
##empirical_rows <- c("dolphins", "jpr", "NNS", "pira", "sn_auth")

ht <- 7
wd <- 14
if(save_plots) {
    pdf("./img/compare-CV-with-u.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mfrow = c(1, 2))
plot(nstages ~ CV, data = df, pch = 19, col = "firebrick",
     xlab = "CV of Degree Distribution", ylab = "Number of Stages in Double-Wells Transition")
plot(range_u ~ CV, data = df, pch = 19, col = "steelblue",
     xlab = "CV of Degree Distribution", ylab = "Magnitude of Range of u Over Which Transitions Occur")
if(save_plots) dev.off()
