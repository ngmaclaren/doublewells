## Code by Neil MacLaren 6/7/2022

                                        # The number of species (nodes) in the network
nspecies <- 10 # 20 # 40 #
                                        # The number of time units in the simulation
TU <- 20
                                        # The standard deviation of the noise process
δ <- 0.05
                                        # Δt
dt <- 0.001
                                        # The noise function, to be used inside the diffeq
ξ <- function(n, ...) {rnorm(n, 0, δ*sqrt(dt))}
                                        # carrying capacities
if(nspecies == 10) {
    K <- runif(nspecies, min = 5, max = 6)
} else if(nspecies == 20) {
    K <- runif(nspecies, min = 7.63, max = 9.15)
} else if(nspecies == 40) {
    K <- runif(nspecies, min = 12.89, max = 15.47)
}
                                        # Equilibrium abundances (without intervention)
Nhat <- runif(nspecies, min = 1.5, max = 2.5)
                                        # Competition parameters
C <- matrix(
    data = runif(nspecies^2, min = 0.04, max = 0.16),
    nrow = nspecies, ncol = nspecies
); diag(C) <- 1
                                        # Mortality rate
d <- runif(nspecies, min = 0.15, max = 0.25)
                                        # Critical abundance (Allee constant)
A <- rbeta(nspecies, shape1 = 5, shape2 = 1) * 1.5
                                        # Relative facilitative benefits; sum to 1 for each row
                                        # scaling parameter that feeds into γ_ij
fill_theta <- function(nspecies, prob = 0.75) {
    ## Making some assumptions.
    ## First, I think I want a n x n square matrix where n is the number of species
    ## Second, I think I remove the 75% of edges first and then do the random draws so that the random draws sum to 1
    sel_prob <- runif(nspecies)
    θ <- ifelse(sel_prob <= prob, 0, 1)
    l <- sum(θ)
    x <- MCMCpack::rdirichlet(1, rep(1, l)) 
    θ[which(θ == 1)] <- x
    return(θ)
}
θ <- matrix(nrow = nspecies, ncol = nspecies) 
for(i in 1:nspecies) {
    θ[i, ] <- fill_theta(nspecies)
}
                                        # Facilitation parameters
                                        #
                                        # initial facilitation parameters
γ_0 <- matrix(
    runif(nspecies^2, min = 0.2, max = 1.8),
    nrow = nspecies, ncol = nspecies
); diag(γ_0) <- 1
                                        # final facilitation parameters
fill_gamma <- function(nspecies, γ_i, Nhat) {
    Γ <- sum(γ_i*Nhat)
    Γ/Nhat
}
γ_final <- matrix(nrow = nspecies, ncol = nspecies)
for(i in 1:nspecies) {# making an assumption about the use of γ_0 here
    γ_final[i, ] <- fill_gamma(nspecies, γ_0[i, ], Nhat) * θ[i, ]
}; diag(γ_final) <- 1
# within-iteration γ
update_γ <- function(γ_0, γ_final, M) {# these need to be recalculated in every iteration
    γ_0 + (γ_final - γ_0)*M
}

Ms <- seq(0, 1, length.out = 10) # conditions; 100 is arbitrary. ###Question### follow their scaling method or do the equilibrium analysis? Both just to check?

                                        # intrinsic growth rate
r <- sapply(1:nspecies, function(i) {
    (d[i]*Nhat[i]*K[i]) / ((sum(γ_0[i, ]*Nhat) - A[i])*(K[i] - sum(C[i, ]*Nhat))*Nhat[i])
})

model <- function(noise = NULL, ...) {
    γ_ <- sum(γ*n) # check this
    c_ <- sum(C*n) # check this
    n + (r*n*((γ_/A) - 1)*(1 - (c_/K)) - d*n)*dt + noise
} # abundance of i

totalT <- TU/dt

for(i in 1:length(Ms)) {
    ##i <- 50
    M <- Ms[i]
    γ <- update_γ(γ_0, γ_final, M)
    N <- matrix(nrow = totalT, ncol = nspecies)
    n <- Nhat
    for(t in 1:totalT) {
        N[t, ] <- n
        n <- model(n, dt = dt, r = r, γ = γ, C = C, A = A, K = K, d = d,
                   noise = ξ(nspecies, δ = δ, dt = dt))
    }

                                        # for now, check the integration
    dev.new()
    matplot(1:totalT, N, type = "l", lwd = .75, col = "black", lty = 1,
            xlab = "T", ylab = "N", main = paste("M =", M))
}

