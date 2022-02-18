library(doublewells)

save_plots <- TRUE

## System Setup

                                        # Parameters of the double-well system
r <- c(1, 4, 7)
                                        # Strength of the noise process (sd)
                                        # If dt = dt/2, s = s/sqrt(2)
s <- 0.005
                                        # Perturbation strength, if using
p <- if(s > 0) 3*s else 0.015
                                        # Stress
u <- 0
                                        # Duration of the simulation in time units
τ <- 75
                                        # Step size for integration
dt <- 0.01
                                        # Duration of the simulation in integration step units
T <- τ/dt
                                        # Initial state of x
initialx <- 1 + noise(1, s)

                                        # Number of samples for early warning indicators
nsamples <- 250
                                        # Sample spacing, τ scale
sample_spacing <- .1
                                        # Sample spacing, Δt scale
sample_spacing <- sample_spacing/dt
                                        # Sample indices in X
samples <- seq(T, T - (nsamples*sample_spacing), by = -sample_spacing)
                                        # Lag for autocorrelation, τ scale
lag <- .1
                                        # Lag for autocorrelation, Δt scale
lag <- lag/dt

## Simulation Setup

x <- initialx
                                        # For incrementing the bifurcation parameter, u
                                        # (because there is no D for a single node system)
stepsize <- 1e-1
                                        # Code to stop the simulation
stopval <- if(r[2] == 3) {
               3 - (2/sqrt(3))
           } else if (r[2] == 4 & r[3] == 7) {
               2.3 # get a better value if it becomes important
           } else if(r[2] == 4) {
               (10 - sqrt(13))/3
           } else if(r[2] == 2) {
               1.5 # get better value, but this is close
           } 

                                        # Storage vectors
acs <- numeric()
sds <- numeric()
us <- numeric()

## Main Loop
i <- 1
while(x < stopval) {
                                        # Debugging
    ## print(x)
                                        # Storage vector for value of x at each Δt
    X <- numeric(T)

                                        # Integration
    for(t in 1:T) {
        X[t] <- x
        
        x <- double_well(x, r[1], r[2], r[3], dt, noise(1, s), u)
    }

                                        # Early Warning Indicators
    acs[i] <- cor(X[samples], X[samples - lag])
    sds[i] <- sd(X[samples])
    us[i] <- u

                                        # Iterate u and i
    u <- u + stepsize
    i <- i + 1
    
}

## Plotting

wd <- 14
ht <- 7
if(save_plots) {
    pdf("./img/single-system-ewi.pdf", width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mfrow = c(1, 2))
plot(us, acs, main = "Autocorrelation", xlab = "u", ylab = expression(ac[1](x)),
     pch = 1, col = "steelblue", lwd = 2)
plot(us, sds, main = "Variance", xlab = "u", ylab = expression(sd(x)),
     pch = 2, col = "darkseagreen", lwd = 2)
if(save_plots) dev.off()
