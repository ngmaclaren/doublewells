## Code by Neil MacLaren 3/23/2022

library(parallel)
library(igraph)
library(doublewells)

                                        # Select and make available the networks on which to simulate.
                                        # `empiricals` is a list of networks available in "doublewells"
                                        # originally drawn from the "networkdata" package.
networks <- c(
    "max_entropy", "me_islands", "pref_attach", "scale_free_2", "LFR",
    empiricals
)
data(list = networks)

                                        # Store results in a "sims" dir, creating subdirs as necessary
outdir <- "./data/sims/"
if(!dir.exists("./data")) dir.create("data")
if(!dir.exists(outdir)) dir.create(outdir)

                                        # How many simulation sets should be run?
                                        # Each set runs through the basic simulation on 22 networks.
                                        # With three worker processes, each run takes about an hour
                                        # on a typical office computer.
nruns <- 50

                                        # For reproducibility, set random seeds ahead of time
set.seed(456)
seeds <- sample.int(10000, size = nruns)

                                        # Choose the number of worker processes. Using one less than
                                        # the number available should keep the computer available for
                                        # other use while the simulations are running.
nthreads <- detectCores() - 1

                                        # Package the simulation tasks in a function that can be
                                        # exported to the workers.
run_sim <- function(network, seed) {
    require(igraph)
    require(doublewells)
                                        # Take the seed from the previously generated vector
    set.seed(seed)
                                        # Load the network on the worker
    data(list = c(network))
    g <- get(network)
                                        # And run the simulation with default settings.
    simulation(g)
}

                                        # Main loop
for(i in 1:nruns) {
                                        # For checking work time during the runs
    time1 <- Sys.time()
    print(paste("Sim:", i))
                                        # File name for sim results
    outfile <- paste0(outdir, "/net-var-results-", i, ".rda")
                                        # Make the cluster for this run
    threads <- makeCluster(nthreads)
                                        # Choose the seed for this run
    seed <- seeds[i]
                                        # Make the necessary variables available on the worker
    clusterExport(threads, varlist = list("networks", "seed"))
                                        # Apply the simulation to the networks on the workers
    results <- clusterApply(
        cl = threads,
        x = networks,
        fun = run_sim,
        seed = seed
    )

                                        # Mark metadata
    for(j in 1:length(results)) {
        results[[j]]$results$run <- i
    }
                                        # Name the variable to match the results file, and save it
    assign(paste0("net_var_results_", i), results)
    save(list = paste0("net_var_results_", i), file = outfile)
                                        # Clean up
    stopCluster(threads)
                                        # Print the running time for this set of networks
    time2 <- Sys.time()
    print(time2 - time1)
}
