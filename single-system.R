library(earlywarnings)
library(doublewells)

noise <- function(s) rnorm(1, 0, s)

r1 <- 1
r2 <- 3
r3 <- 5
dt <- 0.01
T <- 5000

stopval <- if(r2 == 3) {
               3 - (2/sqrt(3))
           } else if(r2 == 4) {
               (10 - sqrt(13))/3
           } else if(r2 == 2) {
               1.5 # get better value, but this is close
           }

initialx <- 1

noise_level <- .005
usteps <- T
##usteps <- T/2
                                        # same simulation method as Dakos et al (2012)
U <- seq(0, 10, length.out = usteps) # if r2 = 2, U 0:1.1
##U <- rep(seq(0, 1.1, length.out = usteps), each = 2)

results <- numeric(T)

x <- initialx
    
for(t in 1:T) {
    results[t] <- x
    
    u <- U[t]
    x <- double_well(x, r1, r2, r3, dt, noise(noise_level), u)

    if(x > stopval) {
        results <- results[1:t]
        break
    }
        
}

## Some additional early warning signals
## bews <- bdstest_ews(cbind(1:T, results)) # funny error message
## chews <- ch_ews(cbind(1:T, results))
## ddj <- ddjnonparam_ews(cbind(1:T, results))
## liv <- livpotential_ews(results) # this one is weird
## qda <- qda_ews(results, detrending = "gaussian", bandwidth = 5) ## Error: "non-stationary AR part from CSS"
## sens <- sensitivity_ews(results, indicator = "ar1", detrending = "gaussian") ## seems to suggest a larger window with a filtering bandwidth of 800 is about right? There was also an "island" of almost-as-good at a sliding window of 40% and bandwidth of ~10%, which seems much more reasonable.

## The generic early warning signals appear to "work" in some simulations and not others. It would be possible to investigate systematically how often (over simulation runs) these early warning signals "worked", but for now I'm going to try some of the model-based signals.
gews <- generic_ews(
    results,
    winsize = 40,
    detrending = "no",
    ## detrending = "gaussian", # this seems most appropriate b/c I am adding Gaussian noise
    ## detrending = "loess",
    ## detrending = "linear",
    ## detrending = "first-diff",
    bandwidth = 10,
    AR_n = FALSE,#TRUE,
    powerspectrum = FALSE#TRUE
)
