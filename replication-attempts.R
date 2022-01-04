library(doublewells)
library(earlywarnings)

noise <- function(s) rnorm(1, 0, s)

use <- "lagmethod" # "sdmethod"

## Model parameters
r <- 1
h <- 1
K <- 10
s <- 0.27#0.03
## Simulation parameters
#steplength <- 5000 # for the forcing on c
cps <- seq(1, 2.6771, length.out = 1000) ## a sequence of "forcing" steps on the harvesting parameter
T <- length(cps)#*steplength
dt <- 1
wl <- floor(T/2)#floor(steplength/2) # for calculating the windowed versions of the early warning signals
## Simulation
initialx <- 9
nruns <- 1
results <- vector(mode = "list", length = nruns)
results <- lapply(results, function(l) l <- numeric(T))
for(i in 1:nruns) {
    p <- 1
    cp <- cps[p]
    x <- initialx
    for(t in 1:T) {
        ## if(t %% steplength == 0 & t != T) {
        ##     p <- p + 1
        ##     cp <- cps[p]
        ## }
        results[[i]][[t]] <- x
        x <- harvestmodel(x = x, r = r, K = K, cp = cp, h = h, dt = dt, s = s, noise = TRUE)
        if(x <= 0) x <- 0 ## seems like not the best solution
        p <- p + 1
        cp <- cps[p]
    }
}
gews <- generic_ews(results[[1]], detrending = "gaussian", bandwidth = 5)

## for testing
plot(NULL, xlab = "t", ylab = "x", xlim = c(0, T), ylim = c(0, 10))
for(i in 1:nruns) lines(1:T, y = results[[i]], lwd = .5, col = pal[i])



## Calculate windowed early warning indicator
if(use == "sdmethod") {
    ewresults <- vector(mode = "list", length = nruns)
    ewresults <- lapply(results, function(x) windowed_sdmethod(x, wl))
} else if(use == "lagmethod") {
    ewresults <- vector(mode = "list", length = nruns)
    lag <- 100
    ewresults <- lapply(results, function(x) windowed_lagmethod(x, wl, lag))
}

## Plot results
dev.new(width = 14, height = 8)
##pdf("./img/harvestmodel-withforcing.pdf", width = 14, height = 8)
thislayout <- layout(matrix(c(rep(1, 4), rep(rep(c(2, 3), each = 2), 4)), byrow = TRUE, ncol = 4))
##layout.show(thislayout)
pal <- rainbow(10)

oldmar <- par()$mar

plot.new()
par(mar = oldmar/2, cex.lab = 1.5)
text(.5, .5, labels = "Harvest Model With Forcing on Grazing Rate (c), 10 trials", cex = 2)

par(mar = c(5.1, 4.1, 0, 2.1))
plot(NULL, xlab = "t", ylab = "x", xlim = c(0, T), ylim = c(0, 10))
## steps <- seq(0, T - steplength, by = steplength)
## text(x = steps, y = 0, srt = 90,
##      labels = paste0("c = ", round(cps, 4)),
##      adj = c(0, .5), col = "black", cex = 1)
for(i in 1:nruns) lines(1:T, y = results[[i]], lwd = .5, col = pal[i])


##plot(NULL, xlim = c(0, T), ylim = c(0, max(unlist(sdresults))*1.1), xlab = "t", ylab = "sd(x)")
plot(NULL, xlim = c(0, T),
     ylim = c(min(unlist(ewresults))*.9, max(unlist(ewresults))*1.1),
     xlab = "t", ylab = "cor(x, lag(x))")
## text(##x = steps, y = 0, srt = 90,
##     x = steps, y = .8, srt = 90,
##     labels = paste0("c = ", round(cps, 4)),
##     adj = c(0, .5), col = "black", cex = 1)
if(use == "sdmethod") {
    for(i in 1:nruns) lines((wl+1):T, ewresults[[i]], lwd = .75, col = pal[i])
} else if(use == "lagmethod") {
    for(i in 1:nruns) lines((wl + lag + 1):T, ewresults[[i]], lwd = .75, col = pal[i])
}
par(mar = oldmar)
##dev.off()
