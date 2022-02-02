## Code by Neil MacLaren, 27 Jan 2022

## In this version, use an option to specify whether varying D or u
## Calculate all EW parameters for each simulation                                     <-- I am here
## Clarify any comments etc. to refer to runs, batches, and simulations as the nested typology

library(igraph)
library(doublewells)

calc_ew <- function(X, A, samples, which = "eig") {
    if(which == "I'") {
        ## select five time points that happen before the perturbation.
        I. <- Iadj(X, A, times = samples)
        ew <- I.
    } else if(which == "eig") {
        ## This is the dominant eigenvalue of the covariance matrix at the same samples
        C <- cov(X[samples, ])
        eig <- eigen(C, symmetric = TRUE, only.values = TRUE)[[1]][1]
        ew <- eig
    } else if(which == "I") {
        I <- ape::Moran.I(X[500, ], A)$observed # need to select just the value
        ew <- I
    }
    return(ew)
}

set.seed(1)

ews <- c("I", "I'", "eig")
which_ew <- ews[3]

which_param <- "D" # "u"
    
production_runs <- TRUE # FALSE

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
if(which_param == "u") u[stressnode] <- 1.93

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
    ew <- calc_ew(X, A, samples, which = which_ew)

    colors <- ifelse(u == 0, "steelblue", "firebrick")
    lwds <- ifelse(u == 0, .5, 1)
    ##main <- expression(paste(u[i], "=", 
    main <- paste("Recovery Time =", time_to_recover)
    ##main <- bquote(D == .(D) u[i] == .(u[stressnode]) Recovery Time == .(time_to_recover))
    matplot(
        1:T, X, type = "l", lty = 1, lwd = lwds, col = colors,
        xlab = "t", ylab = expression(x[i]),
        ## main = paste0("D = ", D, "; ", u[i], " = ", u[stressnode],
        ##               "; Time to Recover: ", time_to_recover)
        main = main
    )
    mtext(bquote(u[i] == .(u[stressnode])), adj = .75)
    mtext(bquote(D == .(D)), adj = 0.25)
    ##mtext(paste0(which_ew, " = ", round(ew, 4)), side = 3, line = 0)
    abline(h = bl_mean + 1.96*bl_sd, lty = 1, lwd = 1, col = "slategray")
    abline(v = recovered, lty = 1, lwd = 1, col = "slategray")
} else {
    if(which_param == "D") {
        Ds <- c(.2, .3, .4, .5, .525, .55, .575, .59, .6, .61)
        Us <- rep(0, length(Ds))
    } else if(which_param == "u") {
        Us <- c(1, 1.25, 1.5, 1.6, 1.7, 1.8, 1.9, 1.91, 1.93, 1.95)#c(.1, .5, 1, 1.25, 1.5, 1.6, 1.7, 1.8, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95)
        Ds <- rep(D, length(Us))
    }
    
    stressnode <- select_stressnode(g, "highest")
    nsims <- 100
    results <- vector("list", length(Us)*nsims)
    ## D, nsim, I., rtime
    sim_n <- 1
    for(i in 1:length(Us)) {
        for(j in 1:nsims) {
            ## Stress
            D <- Ds[i]
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
            I <- calc_ew(X, A, samples, which = "I")
            I. <- calc_ew(X, A, samples, which = "I'")
            eig <- calc_ew(X, A, samples, which = "eig")

            result <- data.frame(D = D, u = Us[i], nsim = j, rtime = time_to_recover,
                                 I = I, I. = I., eig = eig)
            results[[sim_n]] <- result
            sim_n <- sim_n + 1
        }
    }

    results <- do.call(rbind, results)

    if(which_param == "D") {
        filename <- "./img/perturbation-production-D.pdf"
    } else if(which_param == "u") {
        filename <- "./img/perturbation-production-u.pdf"
    }
    
    ##dev.new(width = 16, height = 4)
    ##pdf(filename, width = 16, height = 4)
    ##par(mfrow = c(1, 4), cex.axis = 1.25, cex.lab = 1.25)
    if(which_param == "D") {
        plot(I. ~ D, data = results, type = "p", pch = 19, cex = .75, col = "slategray",
             xlab = "D", ylab = "I")#"Recovery Time")
    } else if(which_param == "u") {
        plot(rtime ~ u, data = results, type = "p", pch = 19, cex = .75, col = "slategray",
             xlab = "u", ylab = "Recovery Time")
    }
    ## plot(rtime ~ I, data = results, type = "p", pch = 19, cex = .75, col = "slategray",
    ##      xlab = "Moran's I", ylab = "")
    ## plot(rtime ~ I., data = results, type = "p", pch = 19, cex = .75, col = "slategray",
    ##      xlab = "I'", ylab = "")
    ## plot(rtime ~ eig, data = results, type = "p", pch = 19, cex = .75, col = "slategray",
    ##      xlab = expression(lambda[1]), ylab = "")
    ## dev.off()
}
