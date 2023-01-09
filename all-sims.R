## Code by Neil MacLaren 12/16/2022

                                        # This file generates and saves simulations on all of the
                                        # networks in the package. The powerlaw and dolphins networks
                                        # are excluded because simulations are already available in
                                        # the examples files.
                                        #
                                        # Running time was around 10 hours. Parallelize if re-running.

library(parallel)
library(igraph)
library(doublewells)

networks <- c(
    "erdos_renyi", "er_islands", "barabasi_albert", "LFR", "powerlaw", "fitness",
    empiricals
)
networks <- networks[-which(networks %in% c("powerlaw", "dolphins"))]

outfile_lower <- "./data/allnets-lower-r1.rda"
outfile_upper <- "./data/allnets-upper-r1.rda"

nthreads <- detectCores() - 1

threads <- makeCluster(nthreads)
clusterExport(threads, varlist = list("networks"))

allnets_lower <- clusterApply(
    cl = threads, x = networks,
    fun = function(x) {
        require(igraph)
        require(doublewells)
        data(list = networks)
        set.seed(123)
        g <- get(x)
        simulation(g, check_alts = TRUE, return_histories = TRUE)#, assessment_samples = 1:10)
    }
)

allnets_upper <- clusterApply(
    cl = threads, x = networks,
    fun = function(x) {
        require(igraph)
        require(doublewells)
        data(list = networks)
        set.seed(123)
        g <- get(x)
        simulation(
            g, from_upper = TRUE, D.init = 1, D.stop = 0, stepsize = -5e-3,
            u = rep(-15, vcount(g)), check_alts = TRUE, return_histories = TRUE
            ## , assessment_samples = 1:10
        )}
)
save(allnets_lower, file = outfile_lower)
save(allnets_upper, file = outfile_upper)

stopCluster(threads)
