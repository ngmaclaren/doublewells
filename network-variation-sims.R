## Code by Neil MacLaren 3/23/2022

library(parallel)
library(igraph)
library(doublewells)

## To run the simulation in clusterApply or parLapply we need a function
run_sim <- function(network, seed) {
    require(igraph)
    require(doublewells)
    
    set.seed(seed)
    
    ##net <- networks[i]
    ##print(net)
    data(list = c(network))
    g <- get(network)
    
    simulation(g) ## The main simulation is already available in the package
}

## Select and make available the networks on which to simulate
networks <- c(
    "max_entropy", "me_islands", "pref_attach", "scale_free_2", "LFR",
    empiricals
)
data(list = networks)

## To store results
outdir <- "./data/sims/"
if(!dir.exists("./data")) dir.create("data")
if(!dir.exists(outdir)) dir.create(outdir)

## This is a loop to keep things organized
## Will output each results file separately (so, 22 networks once through)
nruns <- 50

## For reproducibility, set random seeds ahead of time
set.seed(456)
seeds <- sample.int(10000, size = nruns)
nthreads <- detectCores() - 1

for(i in 1:nruns) {
    time1 <- Sys.time()
    
    outfile <- paste0(outdir, "/net-var-results-", i, ".rda")
    print(paste("Sim:", i))

    threads <- makeCluster(nthreads)
    seed <- seeds[i]
    clusterExport(threads, varlist = list("networks", "seed"))
    
    results <- clusterApply(
        cl = threads,
        x = networks,
        fun = run_sim,
        seed = seed
    )

    ##results$results$run <- i
    for(j in 1:length(results)) {
        results[[j]]$results$run <- i
    }

    assign(paste0("net_var_results_", i), results)
    save(list = paste0("net_var_results_", i), file = outfile)

    stopCluster(threads)
    
    time2 <- Sys.time()
    print(time2 - time1)
}
