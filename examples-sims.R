## Code by Neil MacLaren 5/31/2022
## R1 update 12/14/2022

                                        # This file generates and saves simulations on the example
                                        # networks
library(igraph)
library(doublewells)

## outfile_lower <- "./data/examples-lower-r1.rda" # -r1 to differentiate from original 
## outfile_upper <- "./data/examples-upper-r1.rda" # -r1 to differentiate from original 
outfile_lower <- "./data/examples-lower-r1alt.rda" # -r1alt with assessment samples
outfile_upper <- "./data/examples-upper-r1alt.rda" # -r1alt with assessment samples
## outfile_lower <- "./data/examples-lower-r1-withX.rda # if a copy of X is needed

choices <- c("powerlaw", "dolphins")

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
    examples_lower[[i]] <- simulation(
        g, check_alts = TRUE, return_histories = TRUE
        ##, return_X = TRUE # if saving a copy of X for each simulation
      , assessment_samples = 1:10
    )
                                        # Alternate settings are needed to go from high to low.
                                        # In particular, we step D in a negative direction from 1 to 0
                                        # and use u = 15 for all i (there's a `+` in the `simulation`
                                        # code to accommodate bias in either direction; here we use
                                        # negative bias).
    print("upper")
    examples_upper[[i]] <- simulation(
        g,  from_upper = TRUE, D.init = 1, D.stop = 0, stepsize = -5e-3, u = rep(-15, vcount(g)),
        check_alts = TRUE, return_histories = TRUE
      , assessment_samples = 1:10
    )
}

save(examples_lower, file = outfile_lower)

save(examples_upper, file = outfile_upper)
