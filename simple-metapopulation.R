library(igraph)
library(doublewells)

### Network
N <- 100
g <- sample_gnp(N, .06)
A <- as_adj(g, type = "both", sparse = FALSE)

### Parameters for model
## Note these are chosen from Weinan et al's settings but all the same to facilitate easy scaling in N
K <- 10
h <- 1.7 # this is the control parameter
D <- 0.08
s <- 0.02

### Simulation settings
initialx <- rep(1, N) + noise(N, s)
dt <- 0.01
T <- 5000

## Main Loop
x <- initialx
X <- matrix(0, nrow = T, ncol = N)
for(t in 1:T) {
    X[t, ] <- x
    x <- metapopulation(x, K, h, A, D, dt, noise(N, s))
}

matplot(1:T, X, type = "l", col = "black", lwd = .75, lty = 1, xlab = "t", ylab = expression(x[i]))
    
