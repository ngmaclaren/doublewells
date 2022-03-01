library(igraph)
library(doublewells)

## mutualistic <- function(x, A, D,
##                         model_params = list(B = 0.1, K = 5, C = 1, E = 5, H = 0.9, I = 0.1),
##                         dt = 0.01, s = 0.01, noise = rnorm(length(x), mean = 0, sd = s*sqrt(dt))) {
##     B <- model_params$B
##     K <- model_params$K
##     C <- model_params$C
##     E <- model_params$E
##     H <- model_params$H
##     I <- model_params$I

##     x_i <- matrix(rep(x, nrow(A)), byrow = TRUE, nrow = nrow(A))
##     x_j <- matrix(rep(x, ncol(A)), byrow = FALSE, ncol = ncol(A))
##     xi_xj <- (x_i*x_j)/(E + (H*x_i) + (I*x_j))
        
##     deltax <- (B + x*((1 - (x/K))*((x/C) - 1)) + D*colSums(A*xi_xj))*dt + noise
##     nextx <- x + deltax
##     nextx
## }

very_simple <- FALSE # TRUE
save_plots <- FALSE # TRUE

### Network
if(very_simple) {
    N <- 4
    g <- sample_gnm(N, 4)
    ##plot(g)
} else {
    N <- 100
    ##g <- sample_gnp(N, .1)
    g <- sample_fitness_pl(no.of.nodes = N, no.of.edges = N*3, exponent.out = 2)
    g <- get_gcc(g)
    N <- vcount(g)
    plot(g, vertex.label = "", vertex.size = 3)
}
A <- as_adj(g, type = "both", sparse = FALSE)


### Parameters for model
## these are from Kundu et al 2022, included as defaults above
D <- 0.3

### Parameters for the simulation
τ <- 500
s <- .05
dt <- 0.01

T <- τ/dt
s <- s*sqrt(dt)
initialx <- rep(.1, N) + noise(N, s)

x <- initialx
X <- matrix(0, nrow = T, ncol = N)
for(t in 1:T) {
    X[t, ] <- x
    x <- mutualistic(x, A, D, s = s, noise = noise(N, s))
}

ht = 7
wd = 7
plotsamples <- seq(1, T, by = τ)
if(save_plots) {
    pdf("./img/example-mutualistic.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
matplot(plotsamples, X[plotsamples, ],
        type = "l", lwd = .25, lty = 1, col = "slategray", xlab = "t", ylab = expression(x[i]))
if(save_plots) dev.off()
