double_well <- function(x, r1, r2, r3, dt) {
    "Calculates the next step in a double well simulation assuming full determinism and no coupling. The state variable is x and r1, r2, and r3 are parameters."
    deltax <- (-(x - r1)*(x - r2)*(x - r3))*dt
    nextx <- x + deltax
    nextx
}

## double_well_coupled <- function(xi, r1, r2, r3, D, A, dt) {
##     "Calculates the next step in a double well simulation assuming full determinism and network-based coupling. Variables are as for `double_well`, with D as the coupling strength, A as the adjacency matrix, and `xj` a neighbor of `xi`."
##     deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*sum())*dt
## }


double_well_gao <- function(x, r1, r2, r3, D, A, beta_eff, dt) {# beta_eff has to be calculated from A
    "Calculates the next step in a double well simulation assuming full determinism and Gao method coupling. Variables are as for `double_well`, with D as the coupling strength, `beta_eff` and `x_eff` summarize the degree distribution of the network in different ways."
    x_eff <- mean(rowSums(A)*t(x))/mean(rowSums(A))
    deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*beta_eff*x_eff)*dt
    nextx <- x + deltax
    nextx
}

double_well_noisy <- function(x, r1, r2, r3, dt, noise = rnorm(1)) {
    "Calculates the next step in a double well simulation assuming random noise and no coupling. The state variable is x and r1, r2, and r3 are parameters."
    deltax <- (-(x - r1)*(x - r2)*(x - r3) + noise)*dt
    nextx <- x + deltax
    nextx
}

double_well_coupled <- function(x, r1, r2, r3, D, A, dt) {
    "Calculates the next step in a double well simulation assuming full determinism and network-based coupling. Variables are as for `double_well`, with `x` a row vector of current states, D as the coupling strength, and A as the adjacency matrix."
    ## x is now a row vector
    deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*colSums(A*x))*dt
    nextx <- x + deltax
    nextx
}

double_well_coupled_noisy <- function(x, r1, r2, r3, D, A, dt, noise = rnorm(length(x))) {
    "Calculates the next step in a double well simulation assuming full determinism and network-based coupling. Variables are as for `double_well`, with `x` a row vector of current states, D as the coupling strength, and A as the adjacency matrix."
    ## x is now a row vector
    deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*colSums(A*x) + noise)*dt
    nextx <- x + deltax
    nextx
}
