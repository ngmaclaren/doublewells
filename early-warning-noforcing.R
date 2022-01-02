## Code by Neil G MacLaren, November 2021
## Updated January 2, 2022, to reflect doublewells package updates.
## This simulation tests whether looking at early warning signs measured on the individual nodes can
## predict the transition to the upper equilbirum.
## Currently runs for both "lagmethod" and "sdmethod".

library(doublewells)

##set.seed(12345)

use <- "lagmethod" # "sdmethod"

A <- matrix(1, nrow = 3, ncol = 3)
diag(A) <- 0

r1 <- 1
r2 <- 2
r3 <- 5
dt <- 0.001
T <- 10000
initialx <- rep(r1, 3)
noise_level <- 1
noise <- function(n, noise_level) rnorm(n = n, mean = 0, sd = noise_level)
D <- .32
niter <- 1 # probably want 50 or 100 for production runs

results <- matrix(0, ncol = length(initialx), nrow = T)
for(i in 1:niter) {
    x <- initialx
    for(t in 1:T) {
        results[t, ] <- x
        x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
    }
}

wl <- floor(T/4)
if(use == "sdmethod") {
    ewresults <- apply(results, 2, function(column) windowed_sdmethod(column, wl))
    ylab <- expression(sd(x[i]))
    xvals <- (wl + 1):T
} else if(use == "lagmethod") {
    lag <- 100
    ewresults <- apply(results, 2, function(column) windowed_lagmethod(column, wl, lag))
    ylab <- expression(cor(x[i], lag(x[i])))
    xvals <- (wl + lag + 1):T
}

dev.new(width = 14, height = 7)
thislayout <- layout(matrix(c(1, 2), byrow = TRUE, ncol = 2, nrow = 1))
##layout.show(thislayout)
legendtext <- c(expression(x[1]), expression(x[2]), expression(x[3]))

xcolors <- paste0("gray", c(30, 50, 70))
plot(NULL, xlim = c(0, T), ylim = c(0, 6), xlab = "t", ylab = expression(x[i]))
matlines(results, col = xcolors, lty = 1, lwd = 1)
legend("topleft", legend = legendtext, lty = 1, lwd = 2, col = xcolors, bty = "n")

plot(NULL, xlim = c(0, T), ylim = c(0, max(unlist(ewresults))*1.1),
     xlab = "t", ylab = ylab)
ewcolors <- paste0("darkolivegreen", 2:4)

matlines(xvals, ewresults, lty = 1, lwd = 1, col = ewcolors)
legend("topleft", legend = legendtext, lty = 1, lwd = 2, col = ewcolors, bty = "n")
