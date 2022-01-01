library(doublewells)

varmethod <- function(results, wl) {
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

acmethod <- function(results, wl, lag) {
    "Calculates the correlation of a vector x[t] with x[t+1] along a sliding window."
    windows <- matrix(
        c(seq(1, length(results) - wl - lag), seq(1 + wl + lag, length(results))),
        byrow = FALSE, ncol = 2
    )

    ac_results <- apply(windows[-nrow(windows), ], 1, function(w) {
        vec <- results[w[1]:w[2]]
        nextvec <- results[(w[1] + lag):(w[2] + lag)]
        cor(vec, nextvec, method = "pearson")
    })

    ac_results
}

eigenmethod <- function(results, wl) {
    "For multivariate systems only. Find the dominant eigenvector of the correlation matrix for a multivariate system along a sliding window."
    windows <- matrix(
        c(seq(1, nrow(results) - wl), seq(wl + 1, nrow(results))),
        byrow = FALSE, ncol = 2
    )

    eig_results <- apply(windows, 1, function(w) {
        dat <- results[w[1]:w[2], ]
        corr <- cov(dat, method = "pearson")#cor(dat, method = "pearson")
        eigen(corr, symmetric = TRUE, only.values = TRUE)$values[1]
    })

    eig_results
}

pcamethod <- function(results, wl) {
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
dt <- 0.001
T <- 2.5e5#20000
initialx <- rep(r1, 3)
noise_level <- 2.5
noise <- function(n, noise_level) rnorm(n = n, mean = 0, sd = noise_level)

## here make a vector like the one in the replication attempt
sl <- T/20#2500#1000
##Ds <- c(rep(.1, 5), seq(.1, .7, length.out = 15))
Ds <- seq(.25, .35, length.out = 20)
niter <- 1 # probably want 50 or 100 for production runs
p <- 1
D <- Ds[p]

##results <- vector(mode = "list", length = niter)
##results <- lapply(results, function(m) m <- matrix(0, ncol = length(initialx), nrow = T))
results <- matrix(0, ncol = length(initialx), nrow = T)
for(i in 1:niter) {
    x <- initialx
    for(t in 1:T) {
        if(t %% sl == 0 & t != T) {
            p <- p + 1
            D <- Ds[p]
        }
        
        ##results[[i]][t, ] <- x
        results[t, ] <- x
        x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
        ##x <- double_well(x, r1, r2, r3, dt, noise(length(x), noise_level))
    }
}

##apply(results, 2, summary)
##matplot(results, type = "l", lty = 1, col = 1:3)

wl <- sl/2#1000
lag <- 1000 # get back to 1 from dt
sdresults <- apply(results, 2, function(column) varmethod(column, wl))
## testing
acresults <- apply(results, 2, function(column) acmethod(column, wl, lag))
eigenresults <- eigenmethod(results, wl)
pcaresults <- pcamethod(results, wl)

##dev.new()
pdf("./img/four-indicators.pdf", width = 20, height = 10)
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
matplot((wl+1):(T-lag-1), acresults, type = "l", lty = 1, col = xicolors,
        xlim = c(0, T), ylim = c(-1, 1),
        xlab = "t", ylab = expression(AC1(x[i])), main = "Autocorrelation Method (Lag '1')")
plot((wl+1):T, eigenresults, type = "l", col = xcolor,
     xlim = c(0, T), ylim = c(0, max(eigenresults)*1.1),#c(1, 3),
     xlab = "t", ylab = expression(lambda), main = "Dominant Eigenvalue of the Covariance Matrix")
dev.off()

## dev.new(width = 14, height = 7)
## ##pdf("./img/doublewell-withforcing.pdf", width = 14, height = 7)
## ##thislayout <- layout(matrix(c(1, 2), byrow = TRUE, ncol = 2, nrow = 1))
## thislayout <- layout(matrix(c(rep(1, 4), rep(rep(c(2, 3), each = 2), 4)), byrow = TRUE, ncol = 4))
## ##layout.show(thislayout)

## oldmar <- par()$mar

## plot.new()
## par(mar = oldmar/2)
## text(.5, .5, labels = "Double-Well Triad With Forcing on D", cex = 2)

## par(mar = c(5.1, 4.1, 0, 2.1))
## steps <- seq(0, T - sl, by = 1000)

## xcolors <- c("firebrick", "dodgerblue", "darkorange")#paste0("gray", c(30, 50, 70))
## plot(NULL, xlim = c(0, T), ylim = c(0, 6), xlab = "t", ylab = expression(x[i]))
## text(x = steps, y = 0, srt = 90,
##      labels = paste0("D = ", round(Ds, 4)),
##      adj = c(0, .5), col = "black", cex = .75)
## matlines(results, col = xcolors, lty = 1, lwd = 1)

## ## plot(NULL, xlim = c(0, T), ylim = c(-.25, max(unlist(sdresults))*1.1),
## ##      xlab = "t", ylab = expression(sd(x[i])))
## ## sdcolors <- xcolors##paste0("darkolivegreen", 2:4)
## ## text(x = steps, y = -.25, srt = 90,
## ##      labels = paste0("D = ", round(Ds, 4)),
## ##      adj = c(0, .5), col = "black", cex = .75)
## ## matlines((wl+1):T, sdresults, lty = 1, lwd = 1, col = sdcolors)

## plot(NULL, xlim = c(0, T), ylim = c(-.25, max(eigenresults)*1.1),
##      xlab = "t", ylab = "Dominant Eigenvalue of cor(x)")
## text(x = steps, y = -.25, srt = 90,
##      labels = paste0("D = ", round(Ds, 4)),
##      adj = c(0, .5), col = "black", cex = .75)
## lines((wl+1):T, eigenresults, lty = 1, lwd = 1, col = "chartreuse")


## par(mar = oldmar)
## ##dev.off()
