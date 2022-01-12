library(igraph)
library(doublewells)

## deSolve does the same thing that Prosenjit is doing

## This is the model network
A <- as.matrix(read.csv("./data/matrix.csv", header = FALSE))
colnames(A) <- NULL
g <- graph_from_adjacency_matrix(A, mode = "undirected")

## Implement in some empirical networks

r1 <- 1
r2 <- 3
r3 <- 5
D <- .1
s <- 0
u <- 1.5
dt <- 0.01
T <- 5000

nnodes <-  vcount(g)
initialx <- rep(0.0001, nnodes)
bias <- rep(u, nnodes)
noise <- function(n, s) rnorm(n, 0, s)

results <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx
for(t in 1:T) {
    results[t, ] <- x

    x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(nnodes, s), bias)
}

matplot(
    1:T, results, type = "l",
    lty = 1,
    lwd = .2,
    col = "black",
    xlab = "t", ylab = expression(x[i]), main = "Node States"
)
