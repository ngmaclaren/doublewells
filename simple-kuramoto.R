library(igraph)
library(doublewells)

kuramoto <- function(theta, omega, A, K, N, dt = 0.01, noise = rnorm(length(theta), sd = 1)) {
    theta_j <- matrix(rep(theta, N), byrow = FALSE, ncol = N)
    theta_i <- matrix(rep(theta, N), byrow = TRUE, nrow = N)
    delta_theta <- (omega + (K/N)*colSums(A*sin(theta_j - theta_i)))*dt + noise

    next_theta <- theta + delta_theta
    return(next_theta)
}

very_simple <- TRUE # FALSE
save_plots <- FALSE # TRUE

### Network
if(very_simple) {
    N <- 4
    g <- sample_gnm(N, 4)
    plot(g)
} else {
    N <- 100
    ##g <- sample_gnp(N, .1)
    g <- sample_fitness_pl(no.of.nodes = N, no.of.edges = N*3, exponent.out = 2)
    g <- get_gcc(g)
    N <- vcount(g)
    plot(g, vertex.label = "", vertex.size = 3)
}
A <- as_adj(g, type = "both", sparse = FALSE)

### Model
initial_theta <- ((2*pi)/N)*(1:N)
omega <- rnorm(N, 1, 0.1)

K <- 0.1

s <- 0#0.1

τ <- 50
dt <- 0.01
T <- τ/dt
s <- s*sqrt(dt)

theta <- initial_theta
Theta <- matrix(0, nrow = T, ncol = N)

for(t in 1:T) {
    Theta[t, ] <- theta
    theta <- kuramoto(theta, omega, A, K, N, dt, noise = noise(N, s))
}

matplot(1:T, Theta, type = "l", lwd = 1, lty = 1, col = "olivedrab", xaxt = "n",
        ylab = expression(theta[i]), xlab = "t")
axis(side = 1, at = pretty(c(1, T)), labels = pretty(c(1, T))*dt)
