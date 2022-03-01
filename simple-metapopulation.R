library(igraph)
library(doublewells)

## metapopulation <- function(x, K, h, A, D, dt = 0.01, noise = NULL) {
##     "From Weinans et al 2019. x_i is the abundance at node i, K_i the carrying capacity, h_i the maximum harvesting rate (c in Weinans), d_ij is migration between i and j (symmetric). In Weinans et al, the noise process σ_xi*dW_i is a Wiener process (i.e., Gaussian noise) with mean 0 and variance σ. Weinans et al parameter settings assumed three patches with K_i = [10 13 8], c_i = [3, 2, 2.3], d = [[0, .2, .08] [.2 .08 0] [.08 .08 0] ]. x (N in Weinans), K, and c are vectors, A a symmetric adjacency matrix weighted by a connection parameter D."
##     N <- length(x)

##     x_i <- matrix(rep(x, nrow(A)), byrow = TRUE, nrow = nrow(A))
##     x_j <- matrix(rep(x, ncol(A)), byrow = FALSE, ncol = ncol(A))

##     ## Current <- A*matrix(rep(x, times = N), byrow = FALSE, nrow = N)
##     ## Change <- A*matrix(rep(x, times = N), byrow = TRUE, ncol = N)
##     ## Dispersal <- colSums(D*(Current - Change))
##     ## deltax <- (x*(1 - (x/K)) - ((h*(x^2))/(1 + (x^2))) + Dispersal)*dt + noise
##     deltax <- (x*(1 - (x/K)) - ((h*(x^2))/(1 + (x^2))) + D*colSums(A*(x_j - x_i)))*dt + noise
##     nextx <- x + deltax
##     nextx
## }

very_simple <- FALSE # TRUE
save_plots <- FALSE # TRUE

### Network
if(very_simple) {
    N <- 4
    g <- sample_gnm(N, 4)
} else {
    N <- 100
    ##g <- sample_gnp(N, .1)
    g <- sample_fitness_pl(no.of.nodes = N, no.of.edges = N*3, exponent.out = 2)
    g <- get_gcc(g)
    N <- vcount(g)
    plot(g, vertex.label = "", vertex.size = 3)
}
A <- as_adj(g, type = "both", sparse = FALSE)
##plot(g)

### Parameters for model
## Note these are chosen from Weinan et al's settings but all the same to facilitate easy scaling in N
K <- 10
h <- 1.79 # this is the control parameter
D <- .1
s <- 0.02
### Simulation settings
τ <- 100
dt <- 0.01
T <- τ/dt
s <- s*sqrt(dt)
initialx <- rep(.1, N) + noise(N, s)


## Main Loop
x <- initialx
X <- matrix(0, nrow = T, ncol = N)
for(t in 1:T) {
    X[t, ] <- x
    x <- metapopulation(x, K, h, A, D, dt, noise(N, s))
}

ht = 7
wd = 7
plotsamples <- seq(1, T, by = τ)
if(save_plots) {
    pdf("./img/example-metapopulation.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
matplot(plotsamples, X[plotsamples, ],
        type = "l", col = "slateblue", lwd = .25, lty = 1,
        xlab = "t (dt scale)", ylab = expression(x[i]))
if(save_plots) dev.off()
    
