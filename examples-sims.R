## Code by Neil MacLarenn 5/31/2022

                                        # This file generates and saves simulations on the example
                                        # networks
library(igraph)
library(doublewells)

outfile_lower <- "./data/examples-lower.rda"
outfile_upper <- "./data/examples-upper.rda"

choices <- c("pref_attach", "dolphins")

data(list = choices)
set.seed(123)

examples_lower <- list()
examples_upper <- list()
for(i in 1:length(choices)) {
    print(choices[i])
    
    g <- get(choices[i])
                                        # Standard settings going from low to high, but calculate all
                                        # early warning signals
    print("lower")
    examples_lower[[i]] <- simulation(g, check_alts = TRUE)
                                        # Alternate settings are needed to go from high to low.
                                        # In particular, we step D in a negative direction from 1 to 0
                                        # and use u = -15 for all i
    print("upper")
    examples_upper[[i]] <- simulation(g,  from_upper = TRUE, D.init = 1, D.stop = 0, stepsize = -5e-3, u = rep(-15, vcount(g)), check_alts = TRUE)
}

save(examples_lower, file = outfile_lower)

save(examples_upper, file = outfile_upper)
