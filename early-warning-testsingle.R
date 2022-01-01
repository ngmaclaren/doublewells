## To be more principled about this, I want to get simple slopes of the e.g. sd_results data for sliding windows before the cutoff and after the pig has passed through the snake.

## for now, try switching to the single variable version?
## single variable shows no effect. I think I should save these results (maybe make a few figures?), then try to implement some of the multivariate versions.
## Before I do that, though, I want to be very certain that I'm doing this right.

library(doublewells)

##set.seed(12345)

##A <- matrix(1, nrow = 3, ncol = 3)
##diag(A) <- 0

r1 <- 1
r2 <- 2
r3 <- 5

dt <- 0.001
T <- 20000

##initialx <- rep(r1, 3)
initialx <- 1
##noise_level <- 10
noise_level <- 20
noise <- function(n, noise_level) rnorm(n = n, mean = 0, sd = noise_level)
D <- .2
niter <- 1 # probably want 50 or 100 for production runs

##results <- vector(mode = "list", length = niter)
##results <- lapply(results, function(m) m <- matrix(0, ncol = length(initialx), nrow = T))
results <- numeric(T)

for(i in 1:niter) {
    x <- initialx
    for(t in 1:T) {
        ##results[[i]][t, ] <- x
        results[t] <- x
        ##x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
        x <- double_well(x, r1, r2, r3, dt, noise(length(x), noise_level))
    }
}

##y <- results[[3]][, 1] # this has a breakout
y <- results
summary(y)

thresh <- 2.5
mark <- which(y > thresh)[1]
wl <- floor(.5 * mark)
##wl <- 200 # needs to be a sliding window but this doesn't work because not long enough
windows <- matrix(c(seq(1, length(y) - wl), seq(wl + 1, length(y))), byrow = FALSE, ncol = 2)

## Variance/sd method -- this one works but the others don't seem to.
sd_results <- apply(windows, 1, function(w) {
    vec <- y[w[1]:w[2]]
    vec.mean <- mean(vec)
    resid <- vec - vec.mean
    sd(resid)
    ##sd(vec)
})

breakout <- which(y > thresh)[1]
dev.new()
layout(matrix(c(rep(1, 4), rep(2, 2)), byrow = TRUE, nrow = 3))
plot(1:T, y, type = "l", lwd = 1, lty = 1, col = "black",
     xlim = c(0, length(y)), ylim = c(.75, 6),
     xlab = "t", ylab = "x", main = paste0("Variance, Window Length: ", wl))
abline(h = thresh, lty = 1, lwd = 1, col = "green")
##abline(v = which(y == min(y[14000:15000])), col = "purple")
abline(v = breakout, col = "purple")
plot((wl+1):T, sd_results, type = "l", lwd = 1, lty = 1, col = "red",
     xlim = c(0, length(y)), ylim = c(min(sd_results)*.75, max(sd_results*1.25)),
     xlab = "t", ylab = "sd(x)")
##abline(v = which(y == min(y[14000:15000])), col = "purple")
abline(v = breakout, col = "purple")

### do the simple slopes analysis here
breakout_wl <- 1000
breakout_windows <- matrix(
    c(seq(1, length(sd_results) - breakout_wl),
      seq(breakout_wl + 1, length(sd_results))),
    byrow = FALSE, ncol = 2
)
lms_sd_results <- apply(breakout_windows, 1, function(w) {
    dat <- data.frame(x = w[1]:w[2], y = sd_results[w[1]:w[2]])
    m <- lm(y ~ x, data = dat)
    m$coefficients[2]
})

plot(1:length(lms_sd_results), lms_sd_results, type = "l")
abline(v = (breakout - wl), col = "purple")



## Autocorrelation method; should be correlation of x[t] with x[t-1]
ac_results <- apply(windows[-nrow(windows), ], 1, function(w) {
    vec <- y[w[1]:w[2]]
    nextvec <- y[(w[1] + 1):(w[2] + 1)]
    cor(vec, nextvec, method = "pearson")
})

dev.new()
layout(matrix(c(rep(1, 4), rep(2, 2)), byrow = TRUE, nrow = 3))
plot(1:T, y, type = "l", lwd = 1, lty = 1, col = "black",
     xlim = c(0, length(y)), ylim = c(.75, 6),
     xlab = "t", ylab = "x", main = paste0("Autocorrelation, Window Length: ", wl))
abline(h = thresh, lty = 1, lwd = 1, col = "green")
abline(v = breakout, col = "purple")
##abline(v = which(y == min(y[14000:15000])), col = "purple")
plot((wl+1):(T-1), ac_results, type = "l", lwd = 1, lty = 1, col = "red",
     xlim = c(0, length(y)), ylim = c(min(ar_results)*.75, max(ar_results)*1.25),
     xlab = "t", ylab = "AR1 Coef")
##abline(v = which(y == min(y[14000:15000])), col = "purple")
abline(v = breakout, col = "purple")


## AR1 method
ar_results <- apply(
    windows, 1, function(w) ar(y[w[1]:w[2]], order.max = 1, method = "ols", demean = TRUE)$ar[1]
)

dev.new()
layout(matrix(c(rep(1, 4), rep(2, 2)), byrow = TRUE, nrow = 3))
plot(1:T, y, type = "l", lwd = 1, lty = 1, col = "black",
     xlim = c(0, length(y)), ylim = c(.75, 6),
     xlab = "t", ylab = "x", main = paste0("AR1 Coefficient, Window Length: ", wl))
abline(h = thresh, lty = 1, lwd = 1, col = "green")
abline(v = breakout, col = "purple")
##abline(v = which(y == min(y[14000:15000])), col = "purple")
plot((wl+1):T, ar_results, type = "l", lwd = 1, lty = 1, col = "red",
     xlim = c(0, length(y)), ylim = c(min(ar_results)*.75, max(ar_results)*1.25),
     xlab = "t", ylab = "AR1 Coef")
##abline(v = which(y == min(y[14000:15000])), col = "purple")
abline(v = breakout, col = "purple")


## Huh. I don't think I'm doing this right. Maybe try the autocorrelation and/or variance methods.
