## Code by Neil MacLaren, 27 Jan 2022

library(igraph)
library(doublewells)

noise <- function(n, s) rnorm(n, 0, s)

Iadj <- function(X, A, t = NULL, nsteps = NULL, times = NULL) {
    ## A is a 2D array of connection weights (all w_ij ∈ {0, 1} for now)
    ## t is the chosen time
    ## X is a 2D array where the rows are t_t and the columns are x_i
    ##stopifnot((length(t) == 1 & length(nsteps) == 1) | length(times) > 1)
    if(length(times) == 0) times <- seq(t - nsteps + 1, t)
    if(is.null(nsteps)) nsteps <- length(times)

    X <- X[times, ]
    N <- nrow(A)
    W <- sum(A)

    muI <- colMeans(X)
    deltaX <- apply(X, 2, function(x) x[nsteps] - mean(x))
    
    numerator <- 0
    for(i in 1:nrow(A)) {
        for(j in 1:ncol(A)) {
            if(j >= i) next
            y <- A[i, j]*deltaX[i]*deltaX[j]
            numerator <- numerator + y
        }
    }
    denomenator <- sum((X[nrow(X), ] - muI)^2)

    (N/W)*(numerator/denomenator)
}

## Graph
nnodes <- 100
cprob <- 0.06
g <- sample_gnp(nnodes, cprob)
A <- as_adj(g, type = "both", sparse = FALSE)
hdnodes <- V(g)[which(degree(g) >= 10)]
mdnode <- V(g)[which(degree(g) == max(degree(g)))]

## Node Systems
r <- c(1, 4, 7)
s <- 0.005
D <- 0.565
u <- 3*s
stress <- rep(0, nnodes)

dt <- 0.01
T <- 1000

## Parameters for the stress impulse
impulse_start <- 500
impulse_dur <- 10
impulse_stop <- impulse_start + impulse_dur

## Set/Re-Set X state
initialx <- rep(1, nnodes) + noise(nnodes, s)
X <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx

## Simulation
for(t in 1:T) {
    X[t, ] <- x

    x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), stress)
    if(t >= impulse_start & t < (impulse_start + impulse_dur)) x <- x + u
}

### Analysis
## Recovery
baseline <- X[350:450, mdnode]
bl_mean <- mean(baseline)
bl_sd <- sd(baseline)
stressedstate <- X[impulse_stop, mdnode]
time_to_recover <- which(X[impulse_stop:T, mdnode] <= (bl_mean + 1.96*bl_sd))[1]
recovered <- time_to_recover + impulse_stop
if(is.na(time_to_recover)) time_to_recover <- "Failed to Recover"

## Early Warning
## select five time points that happen before the impulse.
samples <- seq(100, impulse_start, by = 100)
I. <- Iadj(X, A, times = samples)

matplot(1:T, X, type = "l", lty = 1, lwd = .5, col = "steelblue",
        xlab = "t", ylab = expression(x[i]),
        main = paste0("D = ", D, "; Time to Recover: ", time_to_recover))
mtext(paste0("I' = ", round(I., 3)), side = 3, line = 0)
abline(h = bl_mean + 1.96*bl_sd, lty = 1, lwd = 1, col = "slategray")
abline(v = recovered, lty = 1, lwd = 1, col = "slategray")

## matplot(1:T, X[, hdnodes], type = "l", lty = 1, lwd = .5, col = "slategray",#"turquoise",
##         xlab = "t", ylab = expression(x[i]))
## abline(v = impulse_start, col = "turquoise", lty = 2, lwd = 1)
## abline(v = impulse_start + impulse_dur, col = "firebrick", lty = 2, lwd = 1)


## Set/Re-Set X state
Ds <- c(.2, .3, .4, .5, .525, .55, .5525, .555, .5575, .56, .561, .562, .563, .564, .565)
nsims <- 10
results <- vector("list", length(Ds)*nsims)
## D, nsim, I., rtime
sim_n <- 1
for(i in 1:length(Ds)) {
    for(j in 1:nsims) {
        D <- Ds[i]
        
        initialx <- rep(1, nnodes) + noise(nnodes, s)
        X <- matrix(0, nrow = T, ncol = nnodes)
        x <- initialx

        ## Simulation
        for(t in 1:T) {
            X[t, ] <- x

            x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), stress)
            if(t >= impulse_start & t < (impulse_start + impulse_dur)) x <- x + u
        }

        baseline <- X[350:450, mdnode]
        bl_mean <- mean(baseline)
        bl_sd <- sd(baseline)
        stressedstate <- X[impulse_stop, mdnode]
        time_to_recover <- which(X[impulse_stop:T, mdnode] <= (bl_mean + 1.96*bl_sd))[1]

        samples <- seq(100, impulse_start, by = 100)
        I. <- Iadj(X, A, times = samples)

        result <- data.frame(D = D, nsim = j, I. = I., rtime = time_to_recover)
        results[[sim_n]] <- result
        sim_n <- sim_n + 1
    }
}

results <- do.call(rbind, results)
