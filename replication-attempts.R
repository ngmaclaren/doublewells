## touch up these functions and move them to the package

harvestmodel <- function(x, r, K, cp, h, dt, noise = NULL) {
    ## Harvest model used by Dakos et al. 2012
    ## Deterministic system should bifurcate at cp == 2.604
    if(is.null(noise)) {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt
    } else {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt + noise
    }

    nextx <- x + diff
    nextx
}

windowed_sdmethod <- function(results, wl) {
    windows <- matrix(
        c(seq(1, length(results) - wl), seq(wl + 1, length(results))),
        byrow = FALSE, ncol = 2
    )
    
    sd_results <- apply(windows, 1, function(w) {## why is this so slow?
        window <- w[1]:w[2]
        vec <- results[window]
        samplefit <- lm(vec ~ window)
        resid <- samplefit$residuals
        sd(resid)
    })

    sd_results
}

windowed_lagmethod <- function(results, wl, lag) {
    wstarts <- seq(1, length(results) - wl - lag)
    wstops <- wstarts + wl
    windows <- matrix(c(wstarts, wstops), byrow = FALSE, ncol = 2)

    lag_results <- apply(windows, 1, function(w) {
        vec <- results[w[1]:w[2]]
        nextvec <- results[(w[1] + lag):(w[2] + lag)]
        cor(vec, nextvec, method = "pearson")
    })

    lag_results
}

noise <- function(s) rnorm(1, 0, s)

## Model parameters
r <- 1
h <- 1
K <- 10
s <- 0.01

## Simulation parameters
steplength <- 5000 # for the forcing on c
wl <- floor(T/4) # for calculating the windowed versions of the early warning signals

## Control parameter
cps <- seq(1, 2.65, length.out = 5) ## a sequence of "forcing" steps on the harvesting parameter
T <- length(cps)*steplength
dt <- .01

## Simulation
initialx <- 9
nruns <- 10
results <- vector(mode = "list", length = nruns)
results <- lapply(results, function(l) l <- numeric(T))
for(i in 1:nruns) {
    p <- 1
    cp <- cps[p]

    x <- initialx
    for(t in 1:T) {
        if(t %% steplength == 0 & t != T) {
            p <- p + 1
            cp <- cps[p]
        }

        results[[i]][[t]] <- x
        x <- harvestmodel(x, r, K, cp, h, dt, noise(s))
        if(x <= 0) x <- 0 ## seems like not the best solution
    }
}

## Calculate windowed early warning indicator
##sdresults <- vector(mode = "list", length = nruns)
##sdresults <- lapply(results, function(l) windowed_sdmethod(l, wl))

lagresults <- vector(mode = "list", length = nruns)
lag <- 100
lagresults <- lapply(results, function(l) windowed_lagmethod(l, wl, lag))

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
steps <- seq(0, T - steplength, by = steplength)
text(x = steps, y = 0, srt = 90,
     labels = paste0("c = ", round(cps, 4)),
     adj = c(0, .5), col = "black", cex = 1)
for(i in 1:nruns) lines(1:T, y = results[[i]], lwd = .5, col = pal[i])


##plot(NULL, xlim = c(0, T), ylim = c(0, max(unlist(sdresults))*1.1), xlab = "t", ylab = "sd(x)")
plot(NULL, xlim = c(0, T),
     ylim = c(min(unlist(lagresults))*.9, max(unlist(lagresults))*1.1),
     xlab = "t", ylab = "cor(x, lag(x))")
text(##x = steps, y = 0, srt = 90,
    x = steps, y = .8, srt = 90,
     labels = paste0("c = ", round(cps, 4)),
     adj = c(0, .5), col = "black", cex = 1)
##for(i in 1:nruns) lines((wl+1):T, sdresults[[i]], lwd = .75, col = pal[i])
for(i in 1:nruns) lines((wl + lag + 1):T, lagresults[[i]], lwd = .75, col = pal[i])
par(mar = oldmar)
##dev.off()
