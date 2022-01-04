library(earlywarnings)
library(doublewells)

noise <- function(s) rnorm(1, 0, s)

r1 <- 1
r2 <- 2
r3 <- 5
dt <- 0.01
T <- 1000

initialx <- 1

noise_level <- .5
usteps <- T
                                        # same simulation method as Dakos et al (2012)
U <- seq(0, 1.1, length.out = usteps)

results <- numeric(T)

x <- initialx
    
for(t in 1:T) {
    results[t] <- x
    
    u <- U[t]
    x <- double_well(x, r1, r2, r3, dt, noise(noise_level), u)

    if(x > 2.5) {
        results <- results[1:t]
        break
    }
        
}

gews <- generic_ews(
    results,
    winsize = 25,
    ## detrending = "no",
    detrending = "gaussian", # this seems most appropriate b/c I am adding Gaussian noise
    ## detrending = "loess",
    ## detrending = "linear",
    ## detrending = "first-diff",
    bandwidth = 5
)
