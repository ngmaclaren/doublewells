source("double-well-functions.R")

dt <- 0.001
r1 <- 1
r2 <- 2
r3 <- 5
D <- 0.1
β <- 7
## x <- 2

nsteps <- 4000
xs <- numeric(nsteps)
T <- 1:nsteps

## for(t in T) {
##     xs[t] <- x
##     x <- double_well(x, r1, r2, r3, dt)
## }
## for(t in T) {
##     xs[t] <- x
##     x <- double_well_gao(x, r1, r2, r3, D, β)
## }
## for(t in T) {
##     xs[t] <- x
##     x <- double_well_noisy(x, r1, r2, r3, dt, noise = rnorm(1, 0, .0025))
## }

pal <- list(# underworld_palette
    "Ghostpipe" = "#eff8f6", "WhiteToadstool" = "#e4e6d5", "Cupshroom" = "#d0c597",
    "Lichen" = "#95a696", "ShroomBloom" = "#5a5046", "FeatherFern" = "#93c54a",
    "YellowToadstool" = "#e9d466", "FeatherShroom" = "#daa72d", "FalseChant" = "#ec9813",
    "CrackingTable" = "#d67534", "TinyToadstool" = "#db7050", "RedButton" = "#d15935",
    "PurpleButton" = "#805350"
)

initialx <- 2
nruns <- 15

for(run in 1:nruns) {
    x <- initialx

    for(t in T) {
        xs[t] <- x
        x <- double_well_noisy(x, r1, r2, r3, dt, noise = rnorm(1, 0, 10))
    }

    if(run == 1) {
        plot(T, xs, pch = ".", col = pal$Lichen,
             xlab = "t", ylab = "x", main = "Noisy Double Well",
             ylim = c(0.5, 5.5))
    } else {
        points(T, xs, pch = ".", col = pal$Lichen)
    }
}

