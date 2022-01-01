double_well <- function(x, r1, r2, r3, dt, noise = NULL, stress = 0) {
    "Calculates the next step in a double well simulation with no coupling and optional stress and noise. The state variable is x and r1, r2, and r3 are parameters. The `noise` argument if not NULL should be a function that generates a single random value."
    if(is.null(noise)) {
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + stress)*dt
    } else {
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + noise + stress)*dt
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
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*colSums(A*x) + noise + stress)*dt
    }
    
    nextx <- x + deltax
    nextx
}

double_well_gao <- function(x, r1, r2, r3, D, A, beta_eff, dt) {# beta_eff has to be calculated from A
    "Calculates the next step in a double well simulation assuming full determinism and Gao method coupling. Variables are as for `double_well`, with D as the coupling strength, `beta_eff` and `x_eff` summarize the degree distribution of the network in different ways."
    x_eff <- mean(rowSums(A)*t(x))/mean(rowSums(A))
    deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*beta_eff*x_eff)*dt
    nextx <- x + deltax
    nextx
}

harvestmodel <- function(x, r, K, cp, h, dt, noise = NULL) {
    ## Dakos et al. 2012 PLoS ONE
    if(is.null(noise)) {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt
    } else {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt + noise
    }

    nextx <- x + diff
    nextx
}

noise <- function(s) rnorm(1, 0, s)

varmethod <- function(results, wl) {
    "Calculates the demeaned standard deviation of a vector along a sliding window."
    windows <- matrix(
        c(seq(1, length(results) - wl), seq(wl + 1, length(results))),
        byrow = FALSE, ncol = 2
    )
    
    sd_results <- apply(windows, 1, function(w) {
        vec <- results[w[1]:w[2]]
        vec.mean <- mean(vec)
        resid <- vec - vec.mean
        sd(resid)
        ##sd(vec)
    })

    sd_results
}

acmethod <- function(results, wl) {
    "Calculates the correlation of a vector x[t] with x[t+1] along a sliding window."
    windows <- matrix(
        c(seq(1, length(results) - wl), seq(wl + 1, length(results))),
        byrow = FALSE, ncol = 2
    )

    ac_results <- apply(windows[-nrow(windows), ], 1, function(w) {
        vec <- y[w[1]:w[2]]
        nextvec <- y[(w[1] + 1):(w[2] + 1)]
        cor(vec, nextvec, method = "pearson")
    })

    ac_results
}
