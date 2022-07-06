                                        # Use this file interactively
library(igraph)
library(doublewells)

                                        # Calculate Aparicio et al.'s ρ for the two example networks
choices <- c("pref_attach", "dolphins")
data(list = choices)
set.seed(123)

                                        # BA network first
g <- pref_attach
                                        # The calculations for ρ are in the doublewells package
                                        #
                                        # Run the simulation with standard settings for nodes starting
                                        # in the lower state
lowtohigh <- simulation(g, calc_rho = TRUE)
                                        # For nodes starting in the upper state, use the same alternate
                                        # settings as in examples-sims.R
hightolow <- simulation(g, from_upper = TRUE, D.init = 1, D.stop = 0, stepsize = -5e-3,
                        u = rep(-15, vcount(g)), calc_rho = TRUE)
                                        # collect...
rhos <- c(lowtohigh = median(colMeans(lowtohigh)), hightolow = median(colMeans(hightolow)))
                                        # and print to stdout
rhos
                                        # Then the dolphins network
g <- dolphins
lowtohigh <- simulation(g, calc_rho = TRUE)
hightolow <- simulation(g, from_upper = TRUE, D.init = 1.2, D.stop = 0, stepsize = -5e-3,
                        u = rep(-15, vcount(g)), calc_rho = TRUE)
rhos <- c(lowtohigh = median(colMeans(lowtohigh)), hightolow = median(colMeans(hightolow)))
rhos
