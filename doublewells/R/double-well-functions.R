## Models

                                        # A vector of names of empirical data sets stored in this
                                        # package.
empiricals <- c(
    "karate", "pira", "netsci", "jazz","drugusers", "hall", "highschoolboys", "surfersb",
    "weaverbirds", "dolphins", "nestbox", "lizards", "bats", "elephantseals", "tortoises",
    "housefinches", "voles"
)

                                        # Basic double well model
dw <- function(x, r) {
    "This is the basic double well differential equation. For finding minima/maxima."
    -(x - r[1])*(x - r[2])*(x - r[3])
}

double_well <- function(x, r1, r2, r3, dt, noise = NULL, stress = 0) {
    "Calculates the next step in a double well simulation with no coupling and optional stress and noise. The state variable is x and r1, r2, and r3 are parameters. The `noise` argument if not NULL should be a function that generates a single random value."
    if(is.null(noise)) {
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + stress)*dt
    } else {
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + stress)*dt + noise
    }
    
    nextx <- x + deltax
    return(nextx)
}

double_well_coupled <- function(x, r1, r2, r3, D, A, dt, noise = NULL, stress = rep(0, length(x))) {
    "Calculates the next step in a double well simulation and network-based coupling. Variables are generally as for `double_well`, with `x` a vector of current states, D as the coupling strength, and A as the adjacency matrix. The `noise` argument if not NULL should be a function that generates a vector of random values of length equal to the length of x."
    if(is.null(noise)) {
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*colSums(A*x) + stress)*dt
    } else {
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*colSums(A*x) + stress)*dt + noise
    }
    
    nextx <- x + deltax
    nextx
}

noise <- function(n, s, f = rnorm) {
    "A function to generate `n` random noise values from a distribution with standard deviation `s`. Currently set up only for Gaussian noise, but could be expanded to support more distributions."
    f(n, 0, s)
}

## Early Warnings

                                        # These three functions calculate early warning indicators
                                        # based on the a matrix X, which records the state of all
                                        # chosen nodes over the course of the simulation, and a
                                        # vector of samples, which should be row indices of X.
sampled_eigenmethod <- function(X, samples, nodes, var_only = FALSE) {
    "From an output matrix, X, take X at `samples` time steps, make a covariance matrix from that data, and return the dominant eigenvalue of the covariance matrix."
                                        # If only one node, cannot compute a covariance matrix
    if(length(nodes) <= 1) {
        return(NA)
    } else {
                                        # Subset X for the sampled rows (time or integration steps)
                                        # and the chosen nodes.
        X <- X[samples, nodes]
                                        # Calculate the variance/covariance matrix, C.
        C <- cov(X)
                                        # If only interested in the variance, set the off-diagonals
                                        # to zero.
        if(var_only) {
            C[lower.tri(C) | upper.tri(C)] <- 0
        }
                                        # Calculate and return the dominant eigenvalue of C.
        eig <- eigen(C, symmetric = TRUE, only.values = TRUE)[[1]][1]
        return(eig)
    }
}

sampled_sdmethod <- function(X, samples) {
    "Calculate the standard deviation sd(x_i) for output matrix, X, over the given samples. Unlike `sampled_eigenmethod`, this function expects X to be provided with the correct (selected) columns."
                                        # This function could be called for univariate data, in which
                                        # case X should not be a matrix
    if(!is.matrix(X)) {
                                        # and the typical way of calling sd() is sufficient.
        return(sd(X[samples]))
                                        # If a matrix, wrap sd() in a call to apply() to calculate
                                        # the column (node) standard deviations.
    } else return(apply(X[samples, ], 2, sd))
}

sampled_acmethod <- function(X, samples, lag = 1) {
    "Calculate the correlation between x_{i, t0:(T-lag)} with itself at some lag, x_{i, (t0+lag):T}."
                                        # This function could be called for univariate data, in which
                                        # case X should not be a matrix
    if(!is.matrix(X)) {
                                        # and acf() can be calculated from the data vector.
        acf(X[samples], lag.max = lag, plot = FALSE)$acf[2]
    } else {
                                        # Otherwise, the same acf() call needs to be wrapped in
                                        # apply() to independently calculate the correlation within
                                        # each column (node).
        apply(X[samples, ], 2, function(x) acf(x, lag.max = lag, plot = FALSE)$acf[2])
    }
}

## Helper Functions

get_gcc <- function(g) {
    "Return the largest (weakly) connected component of a graph."
    require(igraph)
                                        # Weakly connected is sufficient for directed graphs.
    if(is_directed(g)) {
        comps <- components(g, mode = "weak")
                                        # The mode argument is ignored for undirected graphs, so omit
                                        # it.
    } else comps <- components(g)
                                        # To make the GCC, pull the ID of the component which is the
                                        # largest by number of nodes
    gcc_id <- which.max(comps$csize)
                                        # then the node IDs which belong to that component ID.
    vids <- V(g)[comps$membership == gcc_id]
                                        # Make an induced subgraph with those nodes.
    g <- induced_subgraph(g, vids)
    g
}

CVk <- function(g) {
    "A helper function to calculate the coefficient of variation (standard deviation divided by the mean) of the degree distribution of a graph, g."
    k <- degree(g)
    sd(k)/mean(k)
}

Kendall_correlations <- function(df, constrain = TRUE, cutoff = 15, check_alts = FALSE,
                                 b_param = "Ds") {
    "This function relies on the exact output data frame of the `doublewells::simulation` function: it takes that data frame as input and returns a list with (1) a matrix of mean Kendall correlations, and (2) the standard deviation of those correlations. By default (`constrain = TRUE`) this function will only return computed mean/sd of correlations based on a 'sufficient' amount of data---that is, when the correlation was computed from at least `cutoff` pairs of values for D and the early warning indicator."
                                        # Select the columns with early warning indicator values
    columns <- colnames(df)[!(colnames(df) %in% c("n_lowerstate", "Ds"))]
                                        # Split-Apply-Combine
                                        # SPLIT the data frame so that each grouping has the same
                                        # number of nodes in the lower state.
    dfsplit <- split(df, factor(df$n_lowerstate))
                                        # It makes sense to ignore data where a Kendall correlation
                                        # was calculated only for a few rows.
    if(constrain) dfsplit <- dfsplit[which(sapply(dfsplit, nrow) >= cutoff)]
                                        # APPLY 
    kendalls <- lapply(dfsplit, function(x) {
                                        # warns that sd == 0 in some cases (when there is only one
                                        # row in the split object).
        suppressWarnings(
                                        # Calculate the Kendall (rank order) correlation between the
                                        # bifurcation parameter and each early warning indicator.
            cor(x[, c(b_param, columns)], method = "kendall",
                                        # Ignore NAs, but keeping as much of the data as possible.
                                        # Take only the top row, dropping the first column (which has
                                        # the correlation of D with itself.
                use = "pairwise.complete.obs")[1, -1]
        )
    })
                                        # COMBINE
    kendalls <- as.data.frame(do.call(rbind, kendalls))
                                        # n_lowerstate was stored as a rowname in the kendalls df
    kendalls$n_lowerstate <- as.integer(rownames(kendalls))
                                        # and the number of steps (in linearly increasing D) was the
                                        # number of rows in each split object.
    kendalls$n_steps <- sapply(dfsplit, nrow)
                                        # Collect both the means and st devs of the Kendall
                                        # correlations. As of now, only the means are used down
                                        # stream as these sds are not based on independent
                                        # observations.
    results <- list(
        means = as.matrix(round(colMeans(kendalls[, columns], na.rm = TRUE), 3)),
        sds = suppressWarnings(# warns when sd == 0
            as.matrix(round(apply(kendalls[, columns], 2, sd, na.rm = TRUE), 3))
        )
    )
    return(results)
}

## Algorithm Support Functions

choose_sentinels <- function(g, X, samples, n = 5, state_check = samples[length(samples)],
                             state_cutoff = NULL, upper_state_cutoff = NULL,
                             from_upper = FALSE, low_input = FALSE) {
    "Choose `n` sentinel nodes in a graph `g` based on the AVERAGE state of the x_i in `X` over the M samples."
                                        # With a data matrix X that contains all the integration
                                        # steps, determine the macro state of each node at the
                                        # correct time (e.g., TU = 50).
    x <- X[state_check, ]
                                        # If nodes began in the upper state, mark upper state nodes
                                        # as "available", otherwise mark lower state nodes as
                                        # "available."
    if(from_upper) {
        avail <- V(g)[which(upperstate(x, cutoff = upper_state_cutoff) == 1)]
    } else {
        avail <- V(g)[which(lowerstate(x, cutoff = state_cutoff) == 1)]
    }
                                        # If there are fewer nodes available than n, return all
                                        # available nodes as sentinels.
    if(length(avail) < n) {
        return(avail)
    } else {
                                        # Otherwise, calculate ranks for all available nodes. A
                                        # node's rank is R_i = \sum x_j, where j is a neighbor of i
                                        # and x_j is the mean state of x_j over the M samples.
        ranks <- sapply(avail, function(v) sum(colMeans(as.matrix(X[samples, neighbors(g, v)]))))
                                        # Make a convenience data frame
        df <- data.frame(avail = as.numeric(avail), ranks)
                                        # Sort the df based on whether low input or high input (the
                                        # default) nodes are desired
        if(low_input) {
            df <- df[order(df$ranks, decreasing = FALSE), ]
        } else {
            df <- df[order(df$ranks, decreasing = TRUE), ]
        }
                                        # Return the vector of sentinel nodes
        return(V(g)[df$avail[1:n]])
    }
}

choose_lowrank <- function(g, X, samples, cutoff = 2.5, state_check = samples[length(samples)]) {
    "Choose the lowest ranking of the available nodes as the sentinels. This function uses the same ranking scheme as choose_sentinels. Hard-coded to split ranks at the median."
    x <- X[state_check, ]
    avail <- V(g)[which(lowerstate(x, cutoff = cutoff) == 1)]
    ranks <- sapply(avail, function(v) sum(colMeans(as.matrix(X[samples, neighbors(g, v)]))))
    df <- data.frame(avail = as.numeric(avail), ranks)
    V(g)[df$avail[which(df$ranks < median(df$ranks))]]
}

choose_random <- function(g, sentinels, n = 5) {
    "Choose `n` nodes completely at random from g."
                                        # Select random completely at random from the
                                        # available nodes. All nodes are available.
    random <- V(g)[sample(V(g), n)]
                                        # Return the node vector
    return(random)
}

                                        # R': largecorr
choose_largecorr <- function(X, samples, lower, upper, n = 5, from_upper = FALSE) {
    if(from_upper) {
        stop("`Large Corr' not implemented from upper state.")
    } else {
                                        # Lower-state nodes are available to transition.
        avail <- as.numeric(lower)
                                        # Rank them by Pearson cross-correlation, the sum of
                                        # cross-correlation coefficients with each other node.
        ranks <- sapply(lower, function(u) {
            sum(cor(X[samples, u], X[samples, -u]) * mean(X[samples, -u]))
        })
                                        # Collect ranks and return sentinels as for choose_sentinels
        df <- data.frame(avail, ranks)
        df <- df[order(df$ranks, decreasing = TRUE), ]
    }
    if(nrow(df) < n) {
        return(df$avail)
    } else {
        return(df$avail[1:n])
    }
}

                                        # R'': largesd
                                        # Same basic procedure, but ranking on sd()
choose_largesd <- function(X, samples, lower, upper, n = 5, from_upper = FALSE) {
    if(from_upper) {
        stop("`Large SD' not implemented from upper state.")
    } else {
        avail <- as.numeric(lower)
        ranks <- apply(X[samples, avail], 2, sd)
    }
    df <- data.frame(avail, ranks)
    df <- df[order(df$ranks, decreasing = TRUE), ]
    if(nrow(df) < n) {
        return(df$avail)
    } else {
        return(df$avail[1:n])
    }
}

## Counting Functions

upperstate <- function(x, cutoff = 5.7) {
    "A counting function: return 1 if the node is in the upper state and zero otherwise. The default value of `cutoff` is based on r = [1, 4, 7]."
    ifelse(x >= cutoff, 1, 0)
}

lowerstate <- function(x, cutoff = 2.3) {
    "A counting function: return 1 if the node is in the lower state and zero otherwise. The current default cutoff is set to 2.3, which is /approximately/ the value of x_i at which the node state will be 'pulled' towards the separatrix when r = (1, 4, 7), the current default."
    ifelse(x <= cutoff, 1, 0)
}

## Analysis

                                        # Calculate ρ from Aparicio et al. 2021
                                        # The below calculations must be done for each level of D, so
                                        # this function is meant to be used inside `simulation`. If
                                        # the calc_rho option is chosen, `simulation` outputs a data
                                        # frame with results for each level of D.
calc_rho <- function(X, samples) {
    maxs <- apply(X[samples, ], 2, max)
    mins <- apply(X[samples, ], 2, min)
    means <- apply(X[samples, ], 2, mean)
                                        # the value for each node is (max - min)/mean
    (maxs - mins)/means
}

simulation <- function(g
                                        # These defaults are set for consistent behavior across
                                        # simulations in this study.
                                        # Evenly spaced double wells parameters
                     ,
                       r = c(1, 4, 7)
                                        # Start the bifurcation parameter at a small value.
                     ,
                       D.init = 0.01
                                        # Stress vector, not used in the primary analysis so the
                                        # default is set to zero
                     ,
                       u = rep(0, vcount(g)) 
                                        # Add a small amount of noise
                     ,
                       s = 0.05
                                        # Take a reasonably large sample of x_i values for
                                        # calculating early warning indicators.
                     ,
                       nsamples = 250
                                        # Sample spacing is consistent throughout this study
                     ,
                       sample_spacing = .1
                                        # Lag is chosen based on the correlation length in this data
                     ,
                       lag = 1
                                        # Increase the bifurcation parameter by a small amount.
                                        # Changing this parameter makes more or less data available
                                        # to calculate Kendall correlations.
                     ,
                       stepsize = 5e-3
                                        # Choose a small number of sentinel nodes.
                     ,
                       n_sentinels = 5
                                        # Consider the simulations complete when only 10% (or fewer)
                                        # nodes remain in the lower state at equilibrium.
                     ,
                       node_cutoff = .1*vcount(g)
                                        # Calculate what constitutes being in the lower state
                                        # automatically, based on r.
                     ,
                       state_cutoff = NULL
                                        # Allow 50 TU to come to equilibrium and 25 TU for sampling
                     ,
                       TU = 75
                                        # For integration
                     ,
                       dt = 0.01
                                        # Check whether the system has come to equilibrium, and
                                        # output relevant information. For debugging simulations.
                     ,
                       check_equil = FALSE
                                        # Check for lower state 25 TU before the end of the
                                        # simulation, allowing for 250 samples while the system is
                                        # expected to be at equilibrium.
                     ,
                       state_check = TU - 25
                                        # Should simulation check random and upper state nodes?
                     ,
                       check_alts = FALSE
                                        # Does the simulation start from the upper equilibrium?
                     ,
                       from_upper = FALSE
                     ,
                       D.stop = NULL
                     , 
                       calc_rho = FALSE
                       ) {
    "This function simulates a multi-node double well system over increasing values of the connection strength, D, between nodes. It returns a data frame with the state of the system after each round of integration and the value of each early warning indicator in that iteration. The process takes some time, particularly for larger networks (max attempted has about 375 nodes)."
                                        # Calculate the cutoff value for what constitutes being in
                                        # the lower state based on the value of double wells
                                        # parameter, r. This function finds the local minimum between
                                        # r_1 and r_2; above that local minimum, x_i will be "pulled"
                                        # towards, and then past, the separatrix (r_2).
    if(is.null(state_cutoff)) state_cutoff <- optimize(dw, c(r[1], r[2]), r = r)$minimum
                                        # If checking alternatives, need the local maximum as well
    if(check_alts | from_upper) {
        upper_state_cutoff <- optimize(dw, c(r[2], r[3]), r = r, maximum = TRUE)$maximum
    }
                                        # Record the total number of nodes in g for later use.
    nnodes <- vcount(g)
                                        # And make the adjacency matrix.
    A <- as_adj(g, type = "both", sparse = FALSE)

    ## Node Systems
                                        # Set D
    D <- D.init

    ## Sim Params
                                        # Rescale time-related variables. This is done so that the
                                        # user can think in terms of time units, but the simulation
                                        # can calculate with the results of integration.
    T <- TU/dt
    state_check <- state_check/dt
                                        # Rescaling needs to be done for the noise parameter as well.
    s <- s*sqrt(dt)

    ## Early Warnings Params
    sample_spacing <- sample_spacing/dt
                                        # It is easiest to select samples once by "counting
                                        # backwards" from the end of the simulation.
    samples <- seq(from = T, by = -sample_spacing, length.out = nsamples)
                                        # All x_i are initially set to 1 (lower) or 7 (upper)
    if(from_upper) {
        initialx <- rep(max(r), nnodes)
    } else {
        initialx <- rep(min(r), nnodes)
    }
    
    ## Storage Vectors
    n_lowerstate <- numeric()
    Ds <- numeric()
    maxeig <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    maxsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    avgsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    maxac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    avgac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    sentinel_history <- list() # store these and similar lists for diagnostic purposes
    if(check_equil) equils <- list()

                                        # If checking alternatives
    if(check_alts) {
        maxeig$upper <- numeric()
        maxeig$lowrank <- numeric()
        maxeig$random <- numeric()
        if(!from_upper) {
            maxeig$largecorr <- numeric()
            maxeig$largesd <- numeric()
        }
        maxeig$lowinput_sentinel <- numeric()
        maxsd$upper <- numeric()
        maxsd$lowrank <- numeric()
        maxsd$random <- numeric()
        if(!from_upper) {
            maxsd$largecorr <- numeric()
            maxsd$largesd <- numeric()
        }
        maxsd$lowinput_sentinel <- numeric()
        avgsd$upper <- numeric()
        avgsd$lowrank <- numeric()
        avgsd$random <- numeric()
        if(!from_upper) {
            avgsd$largecorr <- numeric()
            avgsd$largesd <- numeric()
        }
        avgsd$lowinput_sentinel <- numeric()
        maxac$upper <- numeric()
        maxac$lowrank <- numeric()
        maxac$random <- numeric()
        if(!from_upper) {
            maxac$largecorr <- numeric()
            maxac$largesd <- numeric()
        }
        maxac$lowinput_sentinel <- numeric()
        avgac$upper <- numeric()
        avgac$lowrank <- numeric()
        avgac$random <- numeric()
        if(!from_upper) {
            avgac$largecorr <- numeric()
            avgac$largesd <- numeric()
        }
        avgac$lowinput_sentinel <- numeric()
        random_history <- list()
        if(!from_upper) {
            largecorr_history <- list()
            largesd_history <- list()
        }
    }

    if(calc_rho) rhos <- list()

    ## Main Simulation Loop
    i <- 1
                                        # The simulation will continue as long as more than 10% of
                                        # nodes are in the lower state when node states are checked
                                        # on the previous run. However, exit criteria are actually
                                        # calculated (and acted upon) inside the loop.
    while(TRUE) {
                                        # Always restart at the same initial values for x_i
        x <- initialx
                                        # and start with a matrix of zeros to store integration
                                        # results.
        X <- matrix(0, nrow = T, ncol = nnodes)
                                        # Main integration loop
        for(t in 1:T) {
                                        # Store the present state
            X[t, ] <- x
                                        # and calculate the next state.
            x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
        }
                                        # If checking equilibrium state, check states of all nodes at
                                        # the TU each /should/ be at equilibrium to start sampling
                                        # (TU - 25) and at the end. In very few cases, some nodes
                                        # will not have reached equilibrium by TU - 25. However, in
                                        # testing it was not found that this discrepancy influenced
                                        # results.
        if(check_equil) {
            at_equil <- X[(TU-25)/dt, ]
            at_end <- X[TU/dt, ]
            equils[[i]] <- data.frame(at_equil, at_end)
        }


        ## Exit condition
                                        # Stop the simulation here (before recording any results in
                                        # this iteration) if EITHER 10% or fewer nodes remain in the
                                        # lower state at equilibrium OR there are not enough lower
                                        # state nodes for the number of sentinels required.
        in_lowerstate <- V(g)[which(lowerstate(X[state_check, ], cutoff = state_cutoff) == 1)]
        if(from_upper | check_alts) {
            in_upperstate <- V(g)[which(upperstate(X[state_check, ],
                                                   cutoff = upper_state_cutoff) == 1)]
        }
        if(from_upper) {
            if(length(in_upperstate) <= node_cutoff) break
        } else {
            if(length(in_lowerstate) <= node_cutoff) break
        }

        if(calc_rho) rhos[[i]] <- calc_rho(X, samples)

        ## Determine sentinels
                                        # Sentinels are chosen based on the average state of their
                                        # neighbors across samples of X.
                                        # These are "High Input" nodes.
        sentinels <- choose_sentinels(
            g, X, samples, n = n_sentinels, state_check = state_check, from_upper = from_upper,
            state_cutoff = state_cutoff, upper_state_cutoff = upper_state_cutoff, low_input = FALSE
        )
                                        # Collect additional node vectors if checking alternatives
        if(check_alts) {
            lowrank <- choose_lowrank(g, X, samples, cutoff = state_cutoff, state_check = state_check)
            random <- choose_random(g, sentinels, n_sentinels)
            if(!from_upper) {
                largecorrs <- choose_largecorr(X, samples, from_upper = from_upper,
                                               lower = in_lowerstate, upper = in_upperstate)
                largesds <- choose_largesd(X, samples, from_upper = from_upper,
                                           lower = in_lowerstate, upper = in_upperstate)
            }
                                        # These are "Low Input" nodes
            lowinput_sentinels <- choose_sentinels(
                g, X, samples, n = n_sentinels, state_check = state_check, from_upper = from_upper,
                state_cutoff = state_cutoff, upper_state_cutoff = upper_state_cutoff, low_input = TRUE
            )
        }
        
        ## Calculate and Store
                                        # Storage lists for node-level calculations that need to be
                                        # aggregated
        sds <- list(); acs <- list()
                                        # All
        maxeig$all[i] <- sampled_eigenmethod(X, samples = samples, nodes = V(g))
        sds$all <- sampled_sdmethod(X, samples)
        acs$all <- sampled_acmethod(X, samples, lag = lag)
                                        # Lower State
        if(length(in_lowerstate) > 0) {
            maxeig$lower[i] <- sampled_eigenmethod(X, samples = samples, nodes = in_lowerstate)
            sds$lower <- sampled_sdmethod(X[, in_lowerstate], samples)
            acs$lower <- sampled_acmethod(X[, in_lowerstate], samples, lag = lag)
        } else {
            maxeig$lower[i] <- NA
            sds$lower <- rep(NA, ncol(X))
            acs$lower <- rep(NA, ncol(X))
        }
                                        # Sentinels
        sentinels_ <- na.omit(sentinels)
        if(length(sentinels_) > 0) {
            maxeig$sentinel[i] <- sampled_eigenmethod(X, samples = samples, nodes = sentinels_)
            sds$sentinel <- sampled_sdmethod(X[, sentinels_], samples)
            acs$sentinel <- sampled_acmethod(X[, sentinels_], samples, lag = lag)
        } else {
            maxeig$sentinel[i] <- NA
            sds$sentinel <- rep(NA, ncol(X))
            acs$sentinel <- rep(NA, ncol(X))
        }
                                        # Alternates
        if(check_alts) {
                                        # Upper State
            if(length(in_upperstate) > 0) {
                maxeig$upper[i] <- sampled_eigenmethod(X, samples = samples, nodes = in_upperstate)
                sds$upper <- sampled_sdmethod(X[, in_upperstate], samples)
                acs$upper <- sampled_acmethod(X[, in_upperstate], samples, lag = lag)
            } else {
                maxeig$upper[i] <- NA
                sds$upper <- rep(NA, ncol(X))
                acs$upper <- rep(NA, ncol(X))
            }
                                        # Low Rank
            if(length(lowrank) > 0) {
                maxeig$lowrank[i] <- sampled_eigenmethod(X, samples = samples, nodes = lowrank)
                sds$lowrank <- sampled_sdmethod(X[, lowrank], samples)
                acs$lowrank <- sampled_acmethod(X[, lowrank], samples, lag = lag)
            } else {
                maxeig$lowrank[i] <- NA
                sds$lowrank <- rep(NA, ncol(X))
                acs$lowrank <- rep(NA, ncol(X))
            }
                                        # Anti-Sentinel (Random)
            maxeig$random[i] <- sampled_eigenmethod(X, samples = samples, nodes = random)
            sds$random <- sampled_sdmethod(X[, random], samples)
            acs$random <- sampled_acmethod(X[, random], samples, lag = lag)
                                        # Large Correlation Sentinels
            if(!from_upper) {
                maxeig$largecorr[i] <- sampled_eigenmethod(X, samples = samples, nodes = largecorrs)
                sds$largecorr <- sampled_sdmethod(X[, largecorrs], samples)
                acs$largecorr <- sampled_acmethod(X[, largecorrs], samples, lag = lag)
            }
                                        # Large Standard Deviation Sentinels
            if(!from_upper) {
                maxeig$largesd[i] <- sampled_eigenmethod(X, samples = samples, nodes = largesds)
                sds$largesd <- sampled_sdmethod(X[, largesds], samples)
                acs$largesd <- sampled_acmethod(X[, largesds], samples, lag = lag)
            }
                                        # Reverse Direction Sentinels
            lowinput_sentinels_ <- na.omit(lowinput_sentinels)
            if(length(lowinput_sentinels_) > 0) {
                maxeig$lowinput_sentinel[i] <- sampled_eigenmethod(X, samples = samples,
                                                                   nodes = lowinput_sentinels_)
                sds$lowinput_sentinel <- sampled_sdmethod(X[, lowinput_sentinels_], samples)
                acs$lowinput_sentinel <- sampled_acmethod(X[, lowinput_sentinels_], samples,
                                                          lag = lag)
            } else {
                maxeig$lowinput_sentinel[i] <- NA
                sds$lowinput_sentinel <- rep(NA, ncol(X))
                acs$lowinput_sentinel <- rep(NA, ncol(X))
            }
        }
                                        # Aggregate
        for(j in 1:length(sds)) {
            maxsd[[names(sds)[j]]][i] <- max(sds[[j]])
            avgsd[[names(sds)[j]]][i] <- mean(sds[[j]])
        }
        for(j in 1:length(acs)) {
            maxac[[names(acs)[j]]][i] <- max(acs[[j]])
            avgac[[names(acs)[j]]][i] <- mean(acs[[j]])
        }
                                        # Store remaining values needed for analysis
        Ds[i] <- D
        n_lowerstate[i] <- length(in_lowerstate)
        sentinel_history[[i]] <- sentinels
        if(check_alts) {
            random_history[[i]] <- random
            if(!from_upper) largecorr_history[[i]] <- largecorrs
        }

        ## Iterate
        D <- D + stepsize
        i <- i + 1

        if(i %% 10 == 0) print(D)

        ## Stopping criteria for from_upper
        if(from_upper) {
            if(!is.null(D.stop)) {
                if(D <= D.stop) break
            } else if(D <= 0) break
        }
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

    if(calc_rho) {
        rho <- do.call(rbind, rhos)
        ##rho$D <- Ds
        return(rho)
    }
    ## Collect data for analysis, and return the results.
    df <- data.frame(
        n_lowerstate = n_lowerstate,
        Ds = Ds
    )
    df <- cbind(df, maxeig, maxsd, avgsd, maxac, avgac)

    return(df)
}
