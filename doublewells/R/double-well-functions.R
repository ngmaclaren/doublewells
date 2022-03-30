## Models

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

harvestmodel <- function(x, r, K, cp, h, dt, s, noise = FALSE) {
    "Calculates the next step in the harvest model system used by Dakos et al. (2012). Dakos et al. (2012) state that in a deterministic system the model should bifurcate at around c = 2.604. Here, c = cp because c is a function in R."
    whitenoise <- function(s) rnorm(1, 0, s)
    
    if(noise) {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt + whitenoise(s)
    } else {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt
    }

    nextx <- x + diff
    nextx
}

metapopulation <- function(x, K, h, A, D, dt = 0.01, noise = NULL) {
    "From Weinans et al 2019. x_i is the abundance at node i, K_i the carrying capacity, h_i the maximum harvesting rate (c in Weinans), d_ij is migration between i and j (symmetric). In Weinans et al, the noise process σ_xi*dW_i is a Wiener process (i.e., Gaussian noise) with mean 0 and variance σ. Weinans et al parameter settings assumed three patches with K_i = [10 13 8], c_i = [3, 2, 2.3], d = [[0, .2, .08] [.2 .08 0] [.08 .08 0] ]. x (N in Weinans), K, and c are vectors, A a symmetric adjacency matrix weighted by a connection parameter D."
    N <- length(x)
                                        # Make matrices for x_i and x_j
                                        # These matrices support the summation below

                                        # For the js, col 1 is x[1], col 2 is x[2], and so on
    x_j <- matrix(rep(x, nrow(A)), byrow = TRUE, nrow = nrow(A))
                                        # For the is, row 1 is x[1], row 2 is x[2], and so on
    x_i <- matrix(rep(x, ncol(A)), byrow = FALSE, ncol = ncol(A))

    # The diff eq model itself. The summation 
    deltax <- (x*(1 - (x/K)) - ((h*(x^2))/(1 + (x^2))) + D*colSums(A*(x_j - x_i)))*dt + noise
    nextx <- x + deltax
    nextx
}

mutualistic <- function(x, A, D, model_params = list(B = 0.1, K = 5, C = 1, E = 5, H = 0.9, I = 0.1),
                        dt = 0.01, s = 0.01, noise = rnorm(length(x), mean = 0, sd = s*sqrt(dt))) {
    "Ref is Kundu et al. 2022."
    B <- model_params$B
    K <- model_params$K
    C <- model_params$C
    E <- model_params$E
    H <- model_params$H
    I <- model_params$I

                                        # Same set up as above, but with more steps to match the model
    x_j <- matrix(rep(x, nrow(A)), byrow = TRUE, nrow = nrow(A))
    x_i <- matrix(rep(x, ncol(A)), byrow = FALSE, ncol = ncol(A))
    xi_xj <- (x_i*x_j)/(E + (H*x_i) + (I*x_j))
        
    deltax <- (B + x*((1 - (x/K))*((x/C) - 1)) + D*colSums(A*xi_xj))*dt + noise
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
                                        # vector of samples, which should be indices of X. It is
                                        # easiest to determine what those samples should be outside
                                        # of these functions (see simulation(), below, for an example.
sampled_eigenmethod <- function(X, samples, nodes, var_only = FALSE) {#A, 
    "From an output matrix, X, with connections between x_i given in A, take X at `samples` time steps, make a covariance matrix from that data, and return the dominant eigenvalue of the covariance matrix."
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

                                        # These two functions are written in a slightly different
                                        # way, in that they expect a matrix X that is already
                                        # subsetted for the right columns/nodes. In the future, these
                                        # functions may be rewritten to match sampled_eigenmethod.
sampled_sdmethod <- function(X, samples) {
    "Calculate the standard deviation sd(x_i) for output matrix, X, over the time steps in a window. Time step `t` is the last step of the window, with the window calculated `(t - wl + 1):t`."
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
    "Calculate the correlation between x_{i, t0:t} with itself at some lag, x_{i, (t0-lag):(t-lag). Time step `t` is the last step of the window, with the window calculated such that `t0` = t - wl + 1."
                                        # This function could be called for univariate data, in which
                                        # case X should not be a matrix
    if(!is.matrix(X)) {
                                        # and cor() can be calculated from the data vector.
        cor(X[samples], X[samples - lag])
    } else {
                                        # Otherwise, the same cor() call needs to be wrapped in
                                        # apply() to independently calculate the correlation within
                                        # each column (node).
        apply(X, 2, function(x) cor(x[samples], x[samples - lag]))
    }
}

## Helper Functions

select_stressnode <- function(g, add_stress_to = NULL) {
    "Choose a random node to which to apply stress. Can constrain choice to either 'high' or 'low' degree nodes, in which case nodes are chosen from the top or bottom quintiles of the degree distribution, respectively. Alternatively, a random node with the 'highest' or 'lowest' (non-zero) degree can be chosen."
    require(igraph)
                                        # This function is not used in the current version of the
                                        # analysis.

                                        # The purpose of this function is to choose a random node
                                        # from among the available options meeting some criteria.

                                        # If there are no criteria, the function chooses randomly
                                        # from among all nodes.
    if(is.null(add_stress_to)) {
        nnodes <- vcount(g)
        selectnode <- sample(1:nnodes, 1)
        return(selectnode)
    } else {
                                        # Otherwise, it sets quantile boundaries on the degree
                                        # distribution: the bottom 20% will be considered low degree,
                                        # the top 20% will be considered high degree, and the minimum
                                        # and maximum degree will be "lowest" and "highest",
                                        # respectively. In any category, the function will choose
                                        # randomly from among those nodes fitting the criteria (for
                                        # example, if two nodes have the highest degree,
                                        # select_stressnode() will choose randomly from between those
                                        # two).
        k <- degree(g)
        breaks <- quantile(k, probs = c(.2, .8))
        if(add_stress_to == "high") {
            poss <- V(g)[which(k >= breaks[2])]
        } else if(add_stress_to == "low") {
            poss <- V(g)[which(k <= breaks[1] & k > 0)]
        } else if(add_stress_to == "highest") {
            poss <- V(g)[which(k == max(k))]
            if(length(poss) == 1) return(poss)
        } else if(add_stress_to == "lowest") {
            poss <- V(g)[which(k == min(k[which(k > 0)]))]
        }
        selectnode <- sample(poss, 1)
        return(selectnode)
    }
}

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
    k <- degree(g)
    sd(k)/mean(k)
}

Kendall_correlations <- function(df, constrain = TRUE, cutoff = 15) {
    "This function relies on the exact output data frame of the simulation() function: it takes that data frame as input and returns a list with (1) a matrix of mean Kendall correlations, and (2) the standard deviation of those correlations. By default (`constrain = TRUE`) this function will only return computed mean/sd of correlations based on a 'sufficient' amount of data---that is, when the correlation was computed from at least `cutoff` pairs of values for D and the early warning indicator."
                                        # Select the columns with early warning indicator values
    columns <- colnames(df)[
        c(grep("_all", colnames(df)),
          grep("_lower", colnames(df)),
          grep("_sentinel", colnames(df)))
    ]
                                        # Split-Apply-Combine
                                        # SPLIT the data frame so that each grouping has the same
                                        # number of nodes in the lower state.
    dfsplit <- split(df, factor(df$n_lowerstate))
                                        # APPLY 
    kendalls <- lapply(dfsplit, function(x) {
                                        # warns that sd == 0 in some cases (when there is only one
                                        # row in the split object).
        suppressWarnings(
                                        # Calculate the Kendall (rank order) correlation between the
                                        # bifurcation parameter and each early warning indicator.
            cor(x[, c("Ds", columns)], method = "kendall",
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
                                        # It makes sense to ignore data where a Kendall correlation
                                        # was calculated only for a few rows.
    if(constrain) kendalls <- kendalls[kendalls$n_steps >= cutoff, ]
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

## Random Networks

generate_network <- function(choice = c("random_regular", "max_entropy", "sphere_surface",
                                        "me_islands", "pref_attach", "small_world"),
                             nnodes = 100, cprob = .06,
                             rr.k = 6,
                             sph.dim = 3, sph.radius = .279,
                             nislands = 5, nbridges = 1,
                             pa.power = 1.5,
                             pa.outdist = c(0, 10, (8+1/3), (6+2/3), 5, (3+1/3), (1+2/3), (5/6))*cprob,
                             sw.dim = 1, sw.nei = 3, sw.p = .1) {
    "Given a chosen network model type, output an instance of that network with parameters that produce a network with 100 nodes and an expected density of about 0.06."
    require(igraph)
                                        # Only one graph will be produced, based on the user's choice.
    choice <- match.arg(choice)
    
    switch(
        choice,
                                        # Random Regular
        random_regular = sample_k_regular(nnodes, rr.k),
                                        # Maximum Entropy, or Erdos-Renyi
        max_entropy = sample_gnp(nnodes, cprob),
                                        # Sphere Surface (via the sample_dot_product() function)
        sphere_surface = sample_dot_product(
            sample_sphere_surface(dim = sph.dim, n = nnodes, radius = sph.radius)),
                                        # A max entropy graph with community structure
        me_islands = sample_islands(
            islands.n = nislands, islands.size = nnodes/nislands,
            islands.pin = cprob*nislands * (100 - nbridges)/100,
            n.inter = nbridges),
                                        # A preferential attachement network
        pref_attach = sample_pa(nnodes, power = pa.power, out.dist = pa.outdist, directed = FALSE),
                                        # And a small world network
        small_world = sample_smallworld(dim = sw.dim, size = nnodes, nei = sw.nei, p = sw.p)
    )
}

## Algorithm Support Functions

choose_sentinels <- function(g, X, samples, n = 5, cutoff = 2.5,
                             state_check = samples[length(samples)]) {
    "Choose `n` sentinel nodes in a graph `g` based on the AVERAGE state of the x_i in `X` over the time step window, length `wl`, ending at `t`."
                                        # By default, this function will check for the "macro" state
                                        # (e.g., is a node in the lower state?) at the end of the
                                        # simulated sequence. However, in the simulation() function
                                        # below, a different time point is passed.
    x <- X[state_check, ]
                                        # Nodes are available to be sentinels if, at the state check
                                        # time point, those nodes are below the given cutoff in x_i.
    avail <- V(g)[which(lowerstate(x, cutoff = cutoff) == 1)]
                                        # Nodes will be ranked based on their AVERAGE STATE over all
                                        # samples of x_i. By default, this function only checks the
                                        # state of nodes that are neighbors of those nodes in the
                                        # lower state. So, a node's default score is the sum of the
                                        # average states of its neighbors.
    ranks <- sapply(avail, function(v) sum(colMeans(as.matrix(X[samples, neighbors(g, v)]))))
                                        # Then, sort nodes based on their score in `ranks`
    df <- data.frame(avail = as.numeric(avail), ranks)
    df <- df[order(df$ranks, decreasing = TRUE), ]
                                        # and choose the highest scoring `n` nodes as the sentinels
    V(g)[df$avail[1:n]]
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

dw <- function(x, r) {
    "This is the basic double well differential equation. For finding minima/maxima."
    -(x - r[1])*(x - r[2])*(x - r[3])
}

simulation <- function(g
                                        # These defaults are set for consistent behavior across
                                        # simulations in this study.

                                        # Evenly spaced double wells parameters
                     , r = c(1, 4, 7)
                                        # Start the bifurcation parameter at a small value.
                     , D.init = 0.01
                                        # Add a small amount of noise
                     , s = 0.05
                                        # Take a reasonably large sample of x_i values for
                                        # calculating early warning indicators.
                     , nsamples = 250
                                        # Sample spacing is consistent throughout this study
                     , sample_spacing = .1
                                        # Lag is chosen based on the correlation length in this data
                     , lag = 0.1
                                        # Increase the bifurcation parameter by a small amount.
                                        # Changing this parameter makes more or less data available
                                        # to calculate Kendall correlations.
                     , stepsize = 5e-3
                                        # Choose a small number of sentinel nodes.
                     , n_sentinels = 5
                                        # Consider the simulations complete when only 10% (or fewer)
                                        # nodes remain in the lower state at equilibrium.
                     , node_cutoff = .1*vcount(g)
                                        # Calculate what constitutes being in the lower state
                                        # automatically, based on r.
                     , state_cutoff = NULL
                                        # Allow 50 TU to come to equilibrium and 25 TU for sampling
                     , TU = 75
                                        # For integration
                     , dt = 0.01
                                        # Check whether the system has come to equilibrium, and
                                        # output relevant information. For debugging simulations.
                     , check_equil = FALSE
                                        # Check for lower state 25 TU before the end of the
                                        # simulation, allowing for 250 samples while the system is
                                        # expected to be at equilibrium.
                     , state_check = TU - 25
                       ) {
    "This function simulates a multi-node double well system over increasing values of the connection strength, D, between nodes. It returns a data frame with the state of the system after each round of integration and the value of each early warning indicator in that iteration. The process takes some time, particularly for larger networks (max attempted has about 375 nodes)."
                                        # Calculate the cutoff value for what constitutes being in
                                        # the lower state based on the value of double wells
                                        # parameter, r. This function finds the local minimum between
                                        # r_1 and r_2; above that local minimum, x_i will be "pulled"
                                        # towards, and then past, the separatrix (r_2).
    if(is.null(state_cutoff)) state_cutoff <- optimize(dw, c(r[1], r[2]), r = r)$minimum
                                        # Record the total number of nodes in g for later use.
    nnodes <- vcount(g)
                                        # And make the adjacency matrix.
    A <- as_adj(g, type = "both", sparse = FALSE)

    ## Node Systems
                                        # Set D
    D <- D.init
                                        # Stress vector, not used in this analysis so set to zero
    u <- rep(0, nnodes) 

    ## Sim Params
                                        # Rescale time-related variables. This is done so that the
                                        # user can think in terms of time units, but the simulation
                                        # can calculate with the results of integration.
    T <- TU/dt
    state_check <- state_check/dt
                                        # Rescaling needs to be done for the noise parameter as well.
    s <- s*sqrt(dt)

    ## Early Warnings Params
    lag <- lag/dt

    ## sample_spacing <- 0.1
    sample_spacing <- sample_spacing/dt
                                        # It is easiest to select samples once by "counting
                                        # backwards" from the end of the simulation.
    samples <- seq(from = T, by = -sample_spacing, length.out = nsamples)
                                        # All x_i are initially set to 1 with no noise
    initialx <- rep(1, nnodes)
                                        # And therefore all nodes are initially in the lower state.
    in_lowerstate <- V(g)

    ## Storage Vectors
    n_lowerstate <- numeric()
    Ds <- numeric()
    maxeig <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    maxsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    avgsd <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    maxac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    avgac <- list(all = numeric(), lower = numeric(), sentinel = numeric())
    sentinel_history <- list()
    if(check_equil) equils <- list()

    ## Main Simulation Loop
    i <- 1
                                        # The simulation will continue as long as more than 10% of
                                        # nodes are in the lower state when node states are checked
                                        # on the previous run. However, exit criteria are actually
                                        # calculated (and acted upon) inside the loop.
    while(length(in_lowerstate) > node_cutoff) {
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
        if(length(in_lowerstate) <= node_cutoff) break
        if(length(in_lowerstate) < n_sentinels) break

        ## Determine sentinels
                                        # Sentinels are chosen based on the average state of their
                                        # neighbors across samples of X.
        sentinels <- choose_sentinels(g, X, samples, n = n_sentinels,
                                      cutoff = state_cutoff, state_check = state_check)

        ## Calculate and Store
                                        # These functions are documented above.
        maxeig$all[i] <- sampled_eigenmethod(X, samples = samples, nodes = V(g))
        maxeig$lower[i] <- sampled_eigenmethod(X, samples = samples, nodes = in_lowerstate)
        maxeig$sentinel[i] <- sampled_eigenmethod(X, samples = samples, nodes = sentinels)

                                        # The standard deviation indicator is the most
                                        # straightforward and is calculated directly here, but in the
                                        # future the function above may be used instead.
        sds <- list()
        sds$all <- apply(X[samples, ], 2, sd)
        sds$lower <- apply(X[samples, in_lowerstate], 2, sd)
        sds$sentinel <- apply(X[samples, sentinels], 2, sd)

        for(j in 1:length(sds)) maxsd[[j]][i] <- max(sds[[j]])
        for(j in 1:length(sds)) avgsd[[j]][i] <- mean(sds[[j]])

        acs <- list()
        acs$all <- sampled_acmethod(X, samples, lag = lag)
        acs$lower <- sampled_acmethod(X[, in_lowerstate], samples, lag = lag)
        acs$sentinel <- sampled_acmethod(X[, sentinels], samples, lag = lag)

        for(j in 1:length(acs)) maxac[[j]][i] <- max(acs[[j]])
        for(j in 1:length(acs)) avgac[[j]][i] <- mean(acs[[j]])
                                        # Store remaining values needed for analysis
        Ds[i] <- D
        n_lowerstate[i] <- length(in_lowerstate)
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
    maxac <- do.call(cbind, maxac)
    colnames(maxac) <- paste("maxac", colnames(maxac), sep = "_")
    avgac <- do.call(cbind, avgac)
    colnames(avgac) <- paste("avgac", colnames(avgac), sep = "_")

    sentinel_history <- do.call(rbind, sentinel_history)

    ## Collect data for analysis, and return the results.
    df <- data.frame(
        n_lowerstate = n_lowerstate, # [-length(n_lowerstate)]
        Ds = Ds
    )
    df <- cbind(df, maxeig, maxsd, avgsd, maxac, avgac)

    if(check_equil) {
        return(list(
            results = df, sentinels = sentinel_history, equil = equils))
    } else {
        return(list(results = df, sentinels = sentinel_history))
    }
}
