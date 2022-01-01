library(doublewells)

varmethod <- function(wl, results) {
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

##set.seed(12345)

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
D <- .2
niter <- 1 # probably want 50 or 100 for production runs

##results <- vector(mode = "list", length = niter)
##results <- lapply(results, function(m) m <- matrix(0, ncol = length(initialx), nrow = T))
results <- matrix(0, ncol = length(initialx), nrow = T)
for(i in 1:niter) {
    x <- initialx
    for(t in 1:T) {
        ##results[[i]][t, ] <- x
        results[t, ] <- x
        x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
        ##x <- double_well(x, r1, r2, r3, dt, noise(length(x), noise_level))
    }
}

##apply(results, 2, summary)
##matplot(results, type = "l", lty = 1, col = 1:3)

wl <- 5000
sdresults <- apply(results, 2, function(column) varmethod(wl, column))

dev.new(width = 14, height = 7)
thislayout <- layout(matrix(c(1, 2), byrow = TRUE, ncol = 2, nrow = 1))
##layout.show(thislayout)

xcolors <- paste0("gray", c(30, 50, 70))
plot(NULL, xlim = c(0, T), ylim = c(0, 6), xlab = "t", ylab = expression(x[i]))
matlines(results, col = xcolors, lty = 1, lwd = 1)

plot(NULL, xlim = c(0, T), ylim = c(0, max(unlist(sdresults))*1.1),
     xlab = "t", ylab = expression(sd(x[i])))
sdcolors <- paste0("darkolivegreen", 2:4)
matlines((wl+1):T, sdresults, lty = 1, lwd = 1, col = sdcolors)
