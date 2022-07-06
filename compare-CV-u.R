## Code by Neil MacLaren, 2/24/2022
## for preprint at https://arxiv.org/abs/2205.11592v1

library(igraph)
library(doublewells)

save_plots <- FALSE # TRUE
use_noise <- FALSE # TRUE
run_simulation <- FALSE # TRUE
x_is_kdiff <- TRUE # FALSE

outfile <- "./data/kdiff-c-results.rda"

choices <- empiricals # list of empirical networks from doublewells package
data(list = choices)

stepsize <- 1e-2 # for incrementing u

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
                                        # Change Monday, April 18, 2022
        k_diff <- quantile(k, probs = 0.95) - quantile(k, probs = 0.05)

        ## Node Systems
        r <- c(1, 2, 5) # double well parameters
        if(use_noise) s <- 0.005 else s <- 0 # sd of noise process
        D <- 0.001 # this is the initial value of D
        p <- if(s > 0) 3*s else 0.015 # perturbation strength
        u <- .01 # this is the initial value
        U <- rep(u, nnodes) # stress vector

        τ <- 30#75 # duration in time units
        dt <- 0.001#0.01 # step size for integration
        T <- τ/dt

        initialx <- rep(.1, nnodes)# initial values of x_i

        ## Simulation Parameters
        x <- initialx # an updating vector of x_i

        cutoff <- .1*nnodes # to stop simulation; production value is .1
        state_cutoff <- 1.5
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
            in_lowerstate <- V(g)[which(lowerstate(x, cutoff = state_cutoff) == 1)]
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
                         CV = CV, kdiff = k_diff,
                         u = us, n_lowerstate = n_lowerstate)
        results[[i]] <- df
    }
    results <- do.call(rbind, results)
    save(results, file = outfile)
} else {
    load(outfile)
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
    c(range_u = find_range(x, stepsize),
      CV = unique(x$CV), kdiff = unique(x$kdiff),
      nstages = length(unique(x$n_lowerstate)))
})

df <- as.data.frame(do.call(rbind, dat))
df$N <- sapply(rownames(df), function(x) vcount(get(x)))
df$kdiff_raw <- sapply(rownames(df), function(x) {
    g <- get(x)
    k <- degree(g)
    max(k) - min(k)
})

if(x_is_kdiff) {
    ht <- 6; wd <- 6.5
    if(save_plots) {
        pdf("./img/compare-kdiff-with-c.pdf", height = ht, width = wd)
    } else dev.new(height = ht, width = wd)
    par(mar = c(5, 6, 2, 2) + 0.1)
    plot(range_u ~ kdiff_raw, data = df, pch = 1, cex = 2, lwd = 1.75, col = "black",
         axes = FALSE,
         xlim = c(0, 100), ylim = range(df$range_u),
         xlab = "",
         ylab = ""
         )
    box()
    axis(1, cex.axis = 1.75, las = 1)
    axis(2, cex.axis = 1.75, las = 1)
    title(ylab = "Magnitude of range of u", line = 4.5, cex.lab = 1.75)
    title(xlab = expression(italic(k)[max] - italic(k)[min]), cex.lab = 1.75, line = 3.5)
    mtext("(a)", line = -1.3, adj = .01, font = 1, cex = 1.5)
    if(save_plots) dev.off()
} else{
    ht <- 7
    wd <- 7
    if(save_plots) {
        pdf("./img/compare-CV-with-c.pdf", height = ht, width = wd)
    } else dev.new(height = ht, width = wd)
    plot(range_u ~ CV, data = df, pch = 1, lwd = 2, col = "black",
         cex = (df$N/100)*2,
         xlab = "CV of Degree Distribution",
         ylab = "Magnitude of Range of u Over Which Transitions Occur")
    if(save_plots) dev.off()
}
