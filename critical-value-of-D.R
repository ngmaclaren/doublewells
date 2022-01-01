library(doublewells)

A <- matrix(1, nrow = 3, ncol = 3)
diag(A) <- 0

r1 <- 1
r2 <- 2
r3 <- 5

dt <- 0.001
T <- 20000

initialx <- rep(r1, 3)
noise_level <- 10
noise <- function(n, noise_level) rnorm(n = n, mean = 0, sd = noise_level)
##Ds <- round(seq(.2, .3, length.out = 10), 2)
D <- .26 # try several values here
niter <- 30

##results <- vector(mode = "list", length = length(Ds))
results <- vector(mode = "list", length = niter)
##names(results) <- Ds
results <- lapply(results, function(m) m <- matrix(0, ncol = length(initialx), nrow = T))

##for(i in 1:length(Ds)) {
    ##D <- Ds[[i]]
for(i in 1:niter) {
    x <- initialx
    for(t in 1:T) {
        results[[i]][t, ] <- x
        x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
    }
}

pal <- rainbow(niter)

dev.new(width = 20, height = 7)
plot(NULL, xlim = c(0, T), ylim = c(r1*.75, r3*1.25), main = "Undirected Triad",
     xlab = expression(t), ylab = expression(x[i]))
legend("topleft", title = paste0("D = ", D, "; ", niter, " Runs"),
       ##legend = c(as.character(Ds), c(expression(x[1]), expression(x[2]), expression(x[3]))),
       legend = c(as.character(1:niter), c(expression(x[1]), expression(x[2]), expression(x[3]))), 
       lty = c(rep(1, niter), 1, 2, 3), lwd = 2, col = c(pal, rep(pal[1], 3)))

for(i in 1:niter) {
    color <- pal[i]
    matlines(results[[i]], col = color, lwd = 2)
}
