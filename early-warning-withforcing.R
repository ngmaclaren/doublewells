## Code by Neil G MacLaren, November 2021
## Updated January 2, 2022, to reflect doublewells package updates.
## This simulation tests whether looking at early warning signs measured on the individual nodes can
## predict the transition to the upper equilbirum.
## Not clear that I have implemented the multivariate methods appropriately. Moving that testing to a different file. 

library(doublewells)

windowed_eigenmethod <- function(results, wl) {
    "For multivariate systems only. Find the dominant eigenvector of the covariance matrix for a multivariate system along a sliding window."
    windows <- matrix(
        c(seq(1, nrow(results) - wl), seq(wl + 1, nrow(results))),
        byrow = FALSE, ncol = 2
    )

    eig_results <- apply(windows, 1, function(w) {
        dat <- results[w[1]:w[2], ]
        covar <- cov(dat, method = "pearson")#cor(dat, method = "pearson")
        eigen(covar, symmetric = TRUE, only.values = TRUE)$values[1]
    })

    eig_results
}

windowed_pcamethod <- function(results, wl) {
    windows <- matrix(
        c(seq(1, nrow(results) - wl), seq(wl + 1, nrow(results))),
        byrow = FALSE, ncol = 2
    )

    pca_results <- apply(windows, 1, function(w) {
        dat <- results[w[1]:w[2], ]
        pc1 <- prcomp(dat)$x[, 1]
        var(pc1)
    })

    pca_results
}
        
##set.seed(12345)

A <- matrix(1, nrow = 3, ncol = 3)
diag(A) <- 0

r1 <- 1
r2 <- 2
r3 <- 5
dt <- 0.01
T <- 10000
initialx <- rep(r1, 3)
noise_level <- 2.5
noise <- function(n, noise_level) rnorm(n = n, mean = 0, sd = noise_level)

nsteps <- 4
sl <- floor(T/nsteps)
Ds <- seq(.1, .3, length.out = nsteps)
niter <- 1
p <- 1
D <- Ds[p]

results <- matrix(0, ncol = length(initialx), nrow = T)
for(i in 1:niter) {
    x <- initialx
    for(t in 1:T) {
        if(t %% sl == 0 & t != T) {
            p <- p + 1
            D <- Ds[p]
        }
        
        results[t, ] <- x
        x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
    }
}

wl <- floor(T/(nsteps*2))
lag <- 1/dt
sdresults <- apply(results, 2, function(column) windowed_sdmethod(column, wl))
acresults <- apply(results, 2, function(column) windowed_lagmethod(column, wl, lag))
eigenresults <- windowed_eigenmethod(results, wl)
pcaresults <- windowed_pcamethod(results, wl)

dev.new()
##pdf("./img/four-indicators.pdf", width = 20, height = 10)
##pdf("./img/four-indicators-run5.pdf", width = 20, height = 10)
layout(matrix(c(rep(1, 4), 2, 3, 4, 5), byrow = FALSE, ncol = 4, nrow = 2))
xicolors <- c("firebrick", "darkslateblue", "gold")
xcolor <- "dodgerblue"
steps <- seq(0, T - sl, by = sl)
matplot(1:T, results, type = "l", lty = 1, col = xicolors,
        xlim = c(0, T), ylim = c(0, 6),
        xlab = "t", ylab = expression(x[i]), main = "Triangle of Double-Well Systems")
text(x = steps, y = 0, srt = 90,
     labels = paste0("D = ", round(Ds, 4)),
     adj = c(0, .5), col = "black", cex = 1)
matplot((wl+1):T, sdresults, type = "l", lty = 1, col = xicolors,
        xlim = c(0, T), ylim = c(0, max(unlist(sdresults))*1.1),
        xlab = "t", ylab = expression(sd(x[i])), main = "Variance Method")
plot((wl+1):T, pcaresults, type = "l", col = xcolor,
     xlim = c(0, T), ylim = c(0, max(pcaresults)*1.1),
     xlab = "t", ylab = "var(PC1)", main = "Variance Method, First Principal Component")
matplot((wl+lag+1):T, acresults, type = "l", lty = 1, col = xicolors,
        xlim = c(0, T), ylim = c(-1, 1),
        xlab = "t", ylab = expression(AC1(x[i])), main = "Autocorrelation Method (Lag '1')")
plot((wl+1):T, eigenresults, type = "l", col = xcolor,
     xlim = c(0, T), ylim = c(0, max(eigenresults)*1.1),#c(1, 3),
     xlab = "t", ylab = expression(lambda), main = "Dominant Eigenvalue of the Covariance Matrix")
##dev.off()
