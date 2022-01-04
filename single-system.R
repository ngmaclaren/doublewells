library(earlywarnings)
library(doublewells)

noise <- function(s) rnorm(1, 0, s)

r1 <- 1
r2 <- 2
r3 <- 5
dt <- 0.01
T <- 1000

initialx <- 1

noise_level <- .01
usteps <- T
                                        # same simulation method as Dakos et al (2012)
U <- seq(0, .85, length.out = usteps)

results <- numeric(T)
for(u in U) {
    x <- initialx
    
    for(t in 1:T) {
        results[t] <- x
        x <- double_well(x, r1, r2, r3, dt, noise(noise_level), u)
    }
}

## These methods do not work in this system---at least not with Gaussian smoothing. Try the others.
gews <- generic_ews(results, winsize = 25, detrending = "gaussian", bandwidth = 5)

## Try the other package functions here. 
