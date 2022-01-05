library(earlywarnings)

harvestmodel <- function(x, r, K, cp, h, dt, s, noise = FALSE) {
    "Calculates the next step in the harvest model system used by Dakos et al. (2012)."
    whitenoise <- function(s) rnorm(1, 0, s)
    
    if(noise) {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt + whitenoise(s)
    } else {
        diff <- ((r*x*(1 - (x/K))) - (cp*(x^2/(x^2 + h^2))))*dt
    }

    nextx <- x + diff
    nextx
}

## Model parameters
r <- 1
h <- 1
K <- 10
s <- 0.27#0.03

## Simulation parameters
cps <- seq(1, 2.6771, length.out = 1000)
T <- length(cps)
dt <- 1
wl <- floor(T/2)

## Simulation
initialx <- 9
nruns <- 1 # 
results <- vector(mode = "list", length = nruns)
results <- lapply(results, function(l) l <- numeric(T))
for(i in 1:nruns) {
    p <- 1
    cp <- cps[p]
    x <- initialx
    for(t in 1:T) {
        results[[i]][[t]] <- x
        x <- harvestmodel(x = x, r = r, K = K, cp = cp, h = h, dt = dt, s = s, noise = TRUE)

        if(x <= 0) {
            results[[i]] <- results[[i]][1:t]
            break
        }
        
        p <- p + 1
        cp <- cps[p]
    }
}

if(nruns == 1) {
    gews <- generic_ews(results[[1]], detrending = "gaussian", bandwidth = 5, powerspectrum = TRUE)
} else {
    plot(NULL, xlab = "t", ylab = "x", xlim = c(0, T), ylim = c(0, 10))
    for(i in 1:nruns) lines(1:T, y = results[[i]], lwd = .5, col = pal[i])
}
