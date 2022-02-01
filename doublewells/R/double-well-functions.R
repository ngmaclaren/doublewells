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
    "Calculates the next step in a double well simulation assuming full determinism and network-based coupling. Variables are as for `double_well`, with `x` a row vector of current states, D as the coupling strength, and A as the adjacency matrix. The `noise` argument if not NULL should be a function that generates a vector of random values of length equal to the length of x."
    ## x is now a row vector
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

sampled_eigmethod <- function(X, A, samples) {
    "From an output matrix, X, with connections between x_i given in A, take x_i at `samples` points in time and return the dominant eigenvalue of the resulting covariance matrix."
    C <- cov(X[samples, ])
    eig <- eigen(C, symmetric = TRUE, only.values = TRUE)[[1]][1]
    return(eig)
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


## Old functions

double_well_gao <- function(x, r1, r2, r3, D, A, beta_eff, dt) {# beta_eff has to be calculated from A
    "Calculates the next step in a double well simulation assuming full determinism and Gao method coupling. Variables are as for `double_well`, with D as the coupling strength, `beta_eff` and `x_eff` summarize the degree distribution of the network in different ways."
    x_eff <- mean(rowSums(A)*t(x))/mean(rowSums(A))
    deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*beta_eff*x_eff)*dt
    nextx <- x + deltax
    nextx
}
