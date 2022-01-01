library(doublewells)

noise <- function(s) rnorm(1, 0, s)

r1 <- 1
r2 <- 2
r3 <- 5
dt <- 0.01
T <- 1000

initialx <- 1

noise_level <- .1
usteps <- 1000
                                        #usteplength <- floor(T/usteps)
U <- seq(.878, .882, length.out = usteps)
                                        #U <- rep(.89, usteps)

sdresults <- data.frame(u = U, sd = NA)

for(u in U) {
    x <- initialx
    i <- 1
                                        #u <- U[i]
    results <- numeric(T)
    
    for(t in 1:T) {
        ## if(t %% usteplength == 0 & t != T) {
        ##     i <- i + 1
        ##     u <- U[i]
        ## }

        results[t] <- x
        x <- double_well(x, r1, r2, r3, dt, noise(noise_level), u)
    }
    ## what is sd(x'), where x' is detrended x, during timesteps 250 to 750?
    markers <- floor(.5*T):T
    xsample <- results[markers]
    samplefit <- lm(xsample ~ markers)
    xp <- samplefit$residuals

    ew <- sd(xp)
    sdresults$sd[sdresults$u == u] <- ew

    ac <- cor(results[markers - (1/dt)], results[markers])
    sdresults$ac <- ac
    
    ##dev.new()
    ##plot(1:T, results, main = paste0("u = ", round(u, 2), " sd(x') = ", round(ew, 5)))

}

dev.new(); par(mfrow = c(1, 2))
plot(sd ~ u, data = sdresults)
plot(ac ~ u, data = sdresults)



## Simple sim to check u
u <- .880
x <- initialx
longT <- T*10
results <- numeric(longT)
for(t in 1:longT) {
    results[t] <- x
    x <- double_well(x, r1, r2, r3, dt, noise(noise_level), u)
}
dev.new()
plot(1:longT, results, main = paste0("u = ", u))
abline(v = c(markers[1], markers[length(markers)]), col = "red", lwd = 2)
