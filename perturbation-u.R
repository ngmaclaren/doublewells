## Code by Neil MacLaren, 27 Jan 2022

## In this version, set D to a relatively high level and increase u.

library(igraph)
library(doublewells)

calc_ew <- function(X, A, samples) {
    if(which_ew == "I'") {
        ## select five time points that happen before the perturbation.
        I. <- Iadj(X, A, times = samples)
        ew <- I.
    } else if(which_ew == "eig") {
        ## This is the dominant eigenvalue of the covariance matrix at the same samples
        C <- cov(X[samples, ])
        eig <- eigen(C, symmetric = TRUE, only.values = TRUE)[[1]][1]
        ew <- eig
    }
    return(ew)
}

ews <- c("I'", "eig")
which_ew <- ews[2]

production_runs <- TRUE

## Graph
generate_new <- FALSE
if(generate_new) {
    nnodes <- 100
    cprob <- 0.06
    g <- sample_gnp(nnodes, cprob)
} else {
    g <- read_graph("./data/testgraph.gml", format = "gml")
    nnodes <- vcount(g)
}
A <- as_adj(g, type = "both", sparse = FALSE)
hdnodes <- V(g)[which(degree(g) >= 10)]
mdnode <- V(g)[which(degree(g) == max(degree(g)))]

## Node Systems
r <- c(1, 4, 7) # double well parameters
s <- 0.005 # sd of noise process
D <- 0.54 # connection strength; this is the only parameter varied in the production runs
p <- 3*s # perturbation strength
u <- rep(0, nnodes) # stress vector

dt <- 0.01
T <- 1000

## Stress
stressnode <- select_stressnode(g, "highest")
u[stressnode] <- 1.95

## Perturbation
p_start <- 500
p_dur <- 10
p_stop <- p_start + p_dur

if(!production_runs) {
    ## Set/Re-Set X state
    initialx <- rep(1, nnodes) + noise(nnodes, s)
    X <- matrix(0, nrow = T, ncol = nnodes)
    x <- initialx

    ## Simulation
    for(t in 1:T) {
        X[t, ] <- x

        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
        if(t >= p_start & t < (p_start + p_dur)) x <- x + p
    }

### Analysis
    ## Recovery
    baseline <- X[350:450, stressnode]
    bl_mean <- mean(baseline)
    bl_sd <- sd(baseline)
    time_to_recover <- which(X[p_stop:T, stressnode] <= (bl_mean + 1.96*bl_sd))[1]
    recovered <- time_to_recover + p_stop
    if(is.na(time_to_recover)) time_to_recover <- "Failed to Recover"

    ## Early Warning
    samples <- seq(100, p_start, by = 100)
    ew <- calc_ew(X, A, samples)

    colors <- ifelse(u == 0, "steelblue", "firebrick")
    lwds <- ifelse(u == 0, .5, 1)
    matplot(1:T, X, type = "l", lty = 1, lwd = lwds, col = colors,
            xlab = "t", ylab = expression(x[i]),
            main = paste0("D = ", D, "; Time to Recover: ", time_to_recover))
    mtext(paste0(which_ew, " = ", round(ew, 4)), side = 3, line = 0)
    abline(h = bl_mean + 1.96*bl_sd, lty = 1, lwd = 1, col = "slategray")
    abline(v = recovered, lty = 1, lwd = 1, col = "slategray")
} else {
    ##Ds <- c(.2, .3, .4, .5, .525, .55, .575, .59, .6, .61)
    Us <- c(.1, .5, 1, 1.25, 1.5, 1.6, 1.7, 1.8, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95)
    stressnode <- select_stressnode(g, "highest")
    nsims <- 100
    results <- vector("list", length(Us)*nsims)
    ## D, nsim, I., rtime
    sim_n <- 1
    for(i in 1:length(Us)) {
        for(j in 1:nsims) {
            ## Stress
            u[stressnode] <- Us[i]
            
            ## Set/Re-Set X state
            initialx <- rep(1, nnodes) + noise(nnodes, s)
            X <- matrix(0, nrow = T, ncol = nnodes)
            x <- initialx

            ## Simulation
            for(t in 1:T) {
                X[t, ] <- x

                x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
                if(t >= p_start & t < (p_start + p_dur)) x <- x + p
            }

            baseline <- X[350:450, stressnode]
            bl_mean <- mean(baseline)
            bl_sd <- sd(baseline)
            time_to_recover <- which(X[p_stop:T, stressnode] <= (bl_mean + 1.96*bl_sd))[1]

            samples <- seq(100, p_start, by = 100)
            ew <- calc_ew(X, A, samples)
            ##I. <- Iadj(X, A, times = samples)

            result <- data.frame(u = Us[i], nsim = j, ew = ew, rtime = time_to_recover)
            results[[sim_n]] <- result
            sim_n <- sim_n + 1
        }
    }

    results <- do.call(rbind, results)

    dev.new(width = 14, height = 7)
    par(mfrow = c(1, 2))
    plot(rtime ~ u, data = results, type = "p", pch = 19, cex = .75, col = "slategray",
         xlab = "u", ylab = "Recovery Time")
    plot(rtime ~ ew, data = results, type = "p", pch = 19, cex = .75, col = "slategray",
         xlab = which_ew, ylab = "Recovery Time")
}
