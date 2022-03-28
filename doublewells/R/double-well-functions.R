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
    "Calculates the next step in a double well simulation and network-based coupling. Variables are as for `double_well`, with `x` a row vector of current states, D as the coupling strength, and A as the adjacency matrix. The `noise` argument if not NULL should be a function that generates a vector of random values of length equal to the length of x."
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

## metapopulation <- function(x, K, h, A, D, dt = 0.01, noise = NULL) {
##     "From Weinans et al 2019. x_i is the abundance at node i, K_i the carrying capacity, h_i the maximum harvesting rate (c in Weinans), d_ij is migration between i and j (symmetric). In Weinans et al, the noise process σ_xi*dW_i is a Wiener process (i.e., Gaussian noise) with mean 0 and variance σ. Weinans et al parameter settings assumed three patches with K_i = [10 13 8], c_i = [3, 2, 2.3], d = [[0, .2, .08] [.2 .08 0] [.08 .08 0] ]. x (N in Weinans), K, and c are vectors, A a symmetric adjacency matrix weighted by a connection parameter D."
##     N <- length(x)
##     Current <- A*matrix(rep(x, times = N), byrow = FALSE, nrow = N)
##     Change <- A*matrix(rep(x, times = N), byrow = TRUE, ncol = N)
##     Dispersal <- colSums(D*(Current - Change))
##     deltax <- (x*(1 - (x/K)) - ((h*(x^2))/(1 + (x^2))) + Dispersal)*dt + noise
##     nextx <- x + deltax
##     nextx
## }

metapopulation <- function(x, K, h, A, D, dt = 0.01, noise = NULL) {
    "From Weinans et al 2019. x_i is the abundance at node i, K_i the carrying capacity, h_i the maximum harvesting rate (c in Weinans), d_ij is migration between i and j (symmetric). In Weinans et al, the noise process σ_xi*dW_i is a Wiener process (i.e., Gaussian noise) with mean 0 and variance σ. Weinans et al parameter settings assumed three patches with K_i = [10 13 8], c_i = [3, 2, 2.3], d = [[0, .2, .08] [.2 .08 0] [.08 .08 0] ]. x (N in Weinans), K, and c are vectors, A a symmetric adjacency matrix weighted by a connection parameter D."
    N <- length(x)

    x_i <- matrix(rep(x, nrow(A)), byrow = TRUE, nrow = nrow(A))
    x_j <- matrix(rep(x, ncol(A)), byrow = FALSE, ncol = ncol(A))

    ## Current <- A*matrix(rep(x, times = N), byrow = FALSE, nrow = N)
    ## Change <- A*matrix(rep(x, times = N), byrow = TRUE, ncol = N)
    ## Dispersal <- colSums(D*(Current - Change))
    ## deltax <- (x*(1 - (x/K)) - ((h*(x^2))/(1 + (x^2))) + Dispersal)*dt + noise
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

    x_i <- matrix(rep(x, nrow(A)), byrow = TRUE, nrow = nrow(A))
    x_j <- matrix(rep(x, ncol(A)), byrow = FALSE, ncol = ncol(A))
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

windowed_sdmethod <- function(results, wl) {
    "Calculate an early warning indicator in a single variable system using the standard deviation of the detrended `x` over a sliding window."
    windows <- matrix(
        c(seq(1, length(results) - wl), seq(wl + 1, length(results))),
        byrow = FALSE, ncol = 2
    )
    
    sd_results <- apply(windows, 1, function(w) {## why is this so slow?
        window <- w[1]:w[2]
        vec <- results[window]
        samplefit <- lm(vec ~ window)
        resid <- samplefit$residuals
        sd(resid)
    })

    sd_results
}

windowed_lagmethod <- function(results, wl, lag) {
    "Calculate an early warning indicator in a single variable system using the correlation of x with itself at a specified lag over a sliding window."
    wstarts <- seq(1, length(results) - wl - lag)
    wstops <- wstarts + wl
    windows <- matrix(c(wstarts, wstops), byrow = FALSE, ncol = 2)

    lag_results <- apply(windows, 1, function(w) {
        vec <- results[w[1]:w[2]]
        nextvec <- results[(w[1] + lag):(w[2] + lag)]
        cor(vec, nextvec, method = "pearson")
    })

    lag_results
}

windowed_acmethod <- function(results, wl) {
    "Calculate an early warning indicator in a single variable system using the coefficient of a lag = 1 autocorrelation model. Modified from Dakos et al. (2012)'s earlywarnings package."
    ## no `mw` because I pass in a window length
    nwindows <- length(results) - wl + 1
    Results <- matrix(NA, nrow = wl, ncol = nwindows)
    for(i in 1:nwindows) Results[, i] <- results[i:(i + wl - 1)]

    acresults <- apply(Results, 2, function(x)
        ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, intercept = FALSE)$ar)

    acresults
}

Iadj <- function(X, A, t = NULL, nsteps = NULL, times = NULL) {
    "Calculate a variant of Moran's I that removes the expected value of x_i. A is a 2D array of connection weights (all w_ij ∈ {0, 1} for now), t is the chosen time, X is a 2D array where the rows are t_t and the columns are x_i."
    if(length(t) > 0 & length(nsteps) > 0) {
        times <- seq(t - nsteps + 1, t)
    } else if(length(times) > 1) {
        nsteps <- length(times)
    }

    stopifnot((length(t) == 1 & length(nsteps) == 1) | length(times) > 1)
    ##if(length(times) == 0) times <- seq(t - nsteps + 1, t)

    X <- X[times, ]
    N <- nrow(A)
    W <- sum(A)/2

    muI <- colMeans(X)
    deltaX <- apply(X, 2, function(x) x[nsteps] - mean(x))
    
    numerator <- 0
    for(i in 1:nrow(A)) {
        for(j in 1:ncol(A)) {
            if(j >= i) next
            y <- A[i, j]*deltaX[i]*deltaX[j]
            numerator <- numerator + y
        }
    }
    denomenator <- sum((X[nrow(X), ] - muI)^2)

    (N/W)*(numerator/denomenator)
}

sampled_eigenmethod <- function(X, samples, nodes, var_only = FALSE) {#A, 
    "From an output matrix, X, with connections between x_i given in A, take X at `samples` time steps, make a covariance matrix from that data, and return the dominant eigenvalue of the covariance matrix."
    X <- X[samples, nodes]
    ##A <- A[nodes, nodes] # why?
    C <- cov(X)
    if(var_only) {
        C[lower.tri(C) | upper.tri(C)] <- 0
    }
    eig <- eigen(C, symmetric = TRUE, only.values = TRUE)[[1]][1]
    return(eig)
}

sampled_MoranI <- function(X, A, nodes, t) {
    "An adjustment to the 'ape' package's Moran's I: for output matrix, X, with connections in adjacency matrix, A, find the value of Moran's I considering only what is effectively an induced subgraph of the original graph. The index is calculated over the states x_i of the selected nodes."
    ape::Moran.I(X[t, nodes], A[nodes, nodes])$observed
}

sampled_sdmethod <- function(X, samples) {
    "Calculate the standard deviation sd(x_i) for output matrix, X, over the time steps in a window. Time step `t` is the last step of the window, with the window calculated `(t - wl + 1):t`."
    if(!is.matrix(X)) {
        return(sd(X[samples]))
    } else return(apply(X[samples, ], 2, sd))
}

sampled_acmethod <- function(X, samples, lag = 1) {
    "Calculate the correlation between x_{i, t0:t} with itself at some lag, x_{i, (t0-lag):(t-lag). Time step `t` is the last step of the window, with the window calculated such that `t0` = t - wl + 1."
    if(!is.matrix(X)) {
        cor(x[samples], x[samples - lag])
    } else {
        apply(X, 2, function(x) cor(x[samples], x[samples - lag]))
    }
}

## Helper Functions

select_stressnode <- function(g, add_stress_to = NULL) {
    "Choose a random node to which to apply stress. Can constrain choice to either 'high' or 'low' degree nodes, in which case nodes are chosen from the top or bottom quintiles of the degree distribution, respectively. Alternatively, a random node with the 'highest' or 'lowest' (non-zero) degree can be chosen."
    require(igraph)
    
    if(is.null(add_stress_to)) {
        nnodes <- vcount(g)
        selectnode <- sample(1:nnodes, 1)
        return(selectnode)
    } else {
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
    require(igraph)
    comps <- components(g)
    gcc_id <- which.max(comps$csize)
    vids <- V(g)[comps$membership == gcc_id]
    g <- induced_subgraph(g, vids)
    g
}

CVk <- function(g) {
    k <- degree(g)
    sd(k)/mean(k)
}

checks <- function(g) {
    cat("Directed: ", is_directed(g), "\n",
        "Weighted: ", is_weighted(g), "\n",
        "Simple: ", is_simple(g), "\n",
        "Bipartite: ", is_bipartite(g), "\n",
        "N Nodes: ", vcount(g), "\n",
        "N Edges: ", ecount(g), "\n",
        "CV_k: ", CVk(g), "\n",
        "Size of GCC: ", vcount(get_gcc(g)), "\n")
}

checkplot <- function(g) plot(g, vertex.label = "", vertex.size = 3)

Kendall_correlations <- function(df, constrain = TRUE, cutoff = 15) {
    "This function relies on the exact output data frame of the simulation() function: it takes that data frame as input and returns a list with (1) a matrix of mean Kendall correlations, and (2) the standard deviation of those correlations. By default (`constrain = TRUE`) this function will only return computed mean/sd of correlations based on a 'sufficient' amount of data---that is, when the correlation was computed from at least `cutoff` pairs of values for D and the early warning indicator."
    ##columns <- colnames(df)[3:length(colnames(df))]
    columns <- colnames(df)[
        c(grep("_all", colnames(df)),
          grep("_lower", colnames(df)),
          grep("_sentinel", colnames(df)))
    ]
    dfsplit <- split(df, factor(df$n_lowerstate))

    kendalls <- lapply(dfsplit, function(x) {
        suppressWarnings(# warns that sd == 0
            cor(x[, c("Ds", columns)], method = "kendall",
                use = "pairwise.complete.obs")[1, -1]
        )
    })
    kendalls <- as.data.frame(do.call(rbind, kendalls))
    kendalls$n_lowerstate <- as.integer(rownames(kendalls))
    kendalls$n_steps <- sapply(dfsplit, nrow)

    if(constrain) kendalls <- kendalls[kendalls$n_steps >= cutoff, ]

    results <- list(
        means = as.matrix(round(colMeans(kendalls[, columns], na.rm = TRUE), 3)),
        sds = suppressWarnings(# warns when sd == 0
            as.matrix(round(apply(kendalls[, columns], 2, sd, na.rm = TRUE), 3))
        )
    )
    return(results)
}

## Random Networks
generate_network <-
    function(choice = c("random_regular", "max_entropy", "sphere_surface",
                        "me_islands", "pref_attach", "small_world"),
             nnodes = 100, cprob = .06,
             rr.k = 6,
             sph.dim = 3, sph.radius = .279,
             nislands = 5, nbridges = 1,
             pa.power = 1.5,
             pa.outdist = c(0, 10, (8+1/3), (6+2/3), 5, (3+1/3), (1+2/3), (5/6))*cprob,
             sw.dim = 1, sw.nei = 3, sw.p = .1) {
    require(igraph)
    
    choice <- match.arg(choice)
    
    switch(choice,
           random_regular = sample_k_regular(nnodes, rr.k),
           max_entropy = sample_gnp(nnodes, cprob),
           sphere_surface = sample_dot_product(
               sample_sphere_surface(dim = sph.dim, n = nnodes, radius = sph.radius)),
           me_islands = sample_islands(
               islands.n = nislands, islands.size = nnodes/nislands,
               islands.pin = cprob*nislands * (100 - nbridges)/100,
               n.inter = nbridges),
           pref_attach = sample_pa(nnodes, power = pa.power, out.dist = pa.outdist, directed = FALSE),
           small_world = sample_smallworld(dim = sw.dim, size = nnodes, nei = sw.nei, p = sw.p)
           )
}

## Algorithm Support Functions
choose_sentinels <- function(g, n, v = V(g)) {
    "Choose `n` sentinel nodes in a graph `g` based on some criteria. Current criterion is only the degree of the node, checking for the `n` highest degree nodes in the list of available nodes,`v`."
    deg <- sort(degree(g, v = v), decreasing = TRUE)
    maxdegs <- head(deg, n)
    nodes <- sample(V(g)[which(degree(g) %in% maxdegs)], n)
    V(g)[nodes]
}

sentinel_ranking <- function(g, x, n = 5, cutoff = 2.5) {
    "Based on the graph, g, and the states of nodes, x, rank the nodes and return a set of n sentinel nodes."
    ## determine which nodes are in the lower state
    avail <- V(g)[which(lowerstate(x, cutoff = cutoff) == 1)]
    ## for each node i in the lower state, calculate a ranking that is the sum of the states of all adjacent nodes. <-- This isn't what I was talking about in the meeting, but I believe it is closer to Naoki's original idea.
    ## ↓ ranks nodes only on their present state, i.e. at t = T
    ranks <- sapply(avail, function(v) sum(x[neighbors(g, v)]))
    ## sort that vector of ranks
    df <- data.frame(avail = as.numeric(avail), ranks)
    df <- df[order(df$ranks, decreasing = TRUE), ]
    ## return a vector of nodes (class igraph.vs), choosing the top n to be the sentinel nodes.
    V(g)[df$avail[1:n]]
}

sentinel_ranking_ts <- function(g, X, samples, n = 5, cutoff = 2.5, state_check = samples[length(samples)],
                                use.neighbors = TRUE) {
    "Choose `n` sentinel nodes in a graph `g` based on the AVERAGE state of the x_i in `X` over the time step window, length `wl`, ending at `t`."
    x <- X[state_check, ]#X[nrow(X), ]
    avail <- V(g)[which(lowerstate(x, cutoff = cutoff) == 1)]
    ## ↓ ranks nodes on their average state over the window
    if(use.neighbors) {
        ranks <- sapply(avail, function(v) sum(colMeans(as.matrix(X[samples, neighbors(g, v)]))))
    } else ranks <- colMeans(X[samples, avail])
    df <- data.frame(avail = as.numeric(avail), ranks)
    df <- df[order(df$ranks, decreasing = TRUE), ]
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
                     , r = c(1, 4, 7)
                     , D.init = 0.01
                     , s = 0.05
                     , nsamples = 250
                     , lag = 0.1
                     , stepsize = 5e-3
                     , n_sentinels = 5
                     , node_cutoff = .1*vcount(g)
                     , state_cutoff = NULL
                     , TU = 75
                     , dt = 0.01
                     , check_equil = FALSE
                     , state_check = TU - 25
                       ) {
    "This function simulates a multi-node double well system over increasing values of the connection strength, D, between nodes. It returns a data frame with the state of the system after each round of integration and the value of each early warning indicator in that iteration."
    if(is.null(state_cutoff)) state_cutoff <- optimize(dw, c(r[1], r[2]), r = r)$minimum
    
    nnodes <- vcount(g)
    A <- as_adj(g, type = "both", sparse = FALSE)

    ## Node Systems
    ##r <- c(1, 4, 7)
    D <- D.init
    u <- rep(0, nnodes) # stress vector, not used in this analysis so set to zero

    ## Sim Params
    T <- TU/dt
    state_check <- state_check/dt

    s <- s*sqrt(dt)

    ## Early Warnings Params
    lag <- lag/dt

    sample_spacing <- 0.1
    sample_spacing <- sample_spacing/dt
    samples <- seq(from = T, by = -sample_spacing, length.out = nsamples) # for all EWs , T - (nsamples*sample_spacing)

    initialx <- rep(1, nnodes)# + noise(nnodes, s)

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
    while(length(in_lowerstate) > node_cutoff) {
        x <- initialx
        X <- matrix(0, nrow = T, ncol = nnodes)

        for(t in 1:T) {
            X[t, ] <- x
            x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
        }

        if(check_equil) {
            at_equil <- X[(TU-25)/dt, ]
            at_end <- X[TU/dt, ]
            equils[[i]] <- data.frame(at_equil, at_end)
        }

        ## Exit condition
        in_lowerstate <- V(g)[which(lowerstate(X[state_check, ], cutoff = state_cutoff) == 1)]
        if(length(in_lowerstate) <= node_cutoff) break
        if(length(in_lowerstate) < n_sentinels) break

        ## Determine sentinels
        sentinels <- sentinel_ranking_ts(g, X, samples, n = n_sentinels,
                                         cutoff = state_cutoff, state_check = state_check)

        ## Calculate and Store
        maxeig$all[i] <- sampled_eigenmethod(X, samples = samples, nodes = V(g))
        maxeig$lower[i] <- sampled_eigenmethod(X, samples = samples, nodes = in_lowerstate)
        maxeig$sentinel[i] <- sampled_eigenmethod(X, samples = samples, nodes = sentinels)

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

    ## Analysis
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

## Old functions

## double_well_gao <- function(x, r1, r2, r3, D, A, beta_eff, dt) {# beta_eff has to be calculated from A
##     "Calculates the next step in a double well simulation assuming full determinism and Gao method coupling. Variables are as for `double_well`, with D as the coupling strength, `beta_eff` and `x_eff` summarize the degree distribution of the network in different ways."
##     x_eff <- mean(rowSums(A)*t(x))/mean(rowSums(A))
##     deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*beta_eff*x_eff)*dt
##     nextx <- x + deltax
##     nextx
## }

## load_network_model <- function(empirical = FALSE, new = FALSE, network = "scale-free") {
##     "Generate or load a saved standard version of one of several network types or empirical networks. Because matching network density across the different random network models takes experimentation, all of these random network models have 100 nodes and appropriate parameters to support a target density of 0.06."
##     if(!empirical) {
##         if(new) {
##             nnodes <- 100
##             cprob <- 0.06
##             if(network == "scale-free") {
##                 outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
##                 g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
##                 return(g)
##             } else if(network == "max-entropy") {
##                 g <- sample_gnp(nnodes, cprob)
##                 return(g)
##             } else if(network == "islands") {
##                 nislands <- 5
##                 nbridges <- 1
##                 base_prob <- cprob*nislands
##                 pin <- base_prob * (100 - nbridges)/100
##                 g <- sample_islands(islands.n = nislands, islands.size = nnodes/nislands,
##                                     islands.pin = pin, n.inter = nbridges)
##                 return(g)
##             } else if(network == "random-regular") {
##                 nnodes <- 100
##                 g <- sample_k_regular(nnodes, 6)
##                 return(g)
##             } else if(network == "sphere-surface") {
##                 g <- sample_dot_product(sample_sphere_surface(dim = 3, n = nnodes, radius = .279))
##                 return(g)
##             } else if(network == "small-world") {
##                 g <- sample_smallworld(dim = 1, size = 100, nei = 3, p = .1)
##                 return(g)
##             }
##         } else {
##             if(network == "scale-free") {
##                 g <- read_graph("../data/scale-free.gml", format = "gml")
##                 return(g)
##             } else if(network == "max-entropy") {
##                 return(read_graph("../data/max-entropy.gml", format = "gml"))
##             } else if(network == "islands") {
##                 return(read_graph("../data/me-islands.gml", format = "gml"))
##             } else if(network == "random-regular") {
##                 return(read_graph("../data/random-regular.gml", format = "gml"))
##             } else if(network == "sphere-surface") {
##                 return(read_graph("../data/sphere-surface.gml", format = "gml"))
##             } else if(network == "../data/small-world") {
##                 return(read_graph("../data/small-world.gml", format = "gml"))
##             }
##         }
##     } else {
##         if(network == "mine") {
##             return(read_graph("../data/mine.gml", format = "gml"))
##         } else if(network == "jpr") {
##             return(read_graph("../data/jpr.gml", format = "gml"))
##         } else if(network == "surfers") {
##             return(read_graph("../data/surfers.gml", format = "gml"))
##         }
##     }
## }
