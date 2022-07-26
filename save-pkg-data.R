library(igraph)
                                        # The doublewells package is not imported because this file
                                        # builds data used for that package.
                                        #
                                        # Empirical networks come from the "networkdata" package, but
                                        # that package doesn't need to be called directly for this
                                        # code to run. poweRlaw is also a dependency (see the
                                        # `sample_powerlaw` function.
##library(networkdata) 

## NB: netsci and others have weights, but downstream analysis will ignore the weights.
## Other networks in this list are also formed from data that originally could be considered weighted,
## but are dichotomized in the network.

                                        # For reproducibility
set.seed(123)
                                        # Should files be saved?
save_files <- TRUE # FALSE
                                        # Should the networks be inspected?
check_networks <- FALSE # TRUE

## Testing and Collection Functions
get_gcc <- function(g) {
    "Return the largest connected component of a graph."
    require(igraph)
                                        # Get all of the components
    comps <- components(g)
                                        # and find the largest.
    gcc_id <- which.max(comps$csize)
                                        # Get the nodes that are members of that component
    vids <- V(g)[comps$membership == gcc_id]
                                        # Make an induced subgraph with just those nodes
    g <- induced_subgraph(g, vids)
                                        # and return that subgraph
    g
}
                                        # Make all networks uniformly undirected and simple, and only
                                        # use the GCC.
convert <- function(g) as.undirected(simplify(get_gcc(g)))
                                        # A function to generate a network from a power law degree
                                        # distribution.
sample_powerlaw <- function(N, alpha = 2, kmin = 1, max.iter = 10) {
    require(igraph)
    require(poweRlaw)
    if(alpha > 2.5) stop("Power law exponent too large: unlikely to make acceptable graph.")
    if(alpha < 1.8) stop("Power law exponent too small: unlikely to make acceptable graph.")
    if(alpha < 1.9) warning("Power law exponent may be too small to make acceptale graph. Try increasing `max.iter` or `alpha`.")
    
    test <- FALSE
    iter <- 1

    while(!test) {
        k <- rpldis(n = N, xmin = kmin, alpha = alpha)
        if(max(k) < N) {
            test <- tryCatch(
                sample_degseq(
                    k,
                    method = "vl",
                ),
                error = function(cond) return(FALSE)
            )
        }

        if(is.igraph(test)) {
            g <- test
            test <- TRUE
            return(g)
        } else {
            iter <- iter + 1
            if(iter == max.iter) {
                test <- TRUE
                message("Max iter reached without acceptable degree sequence.")
            }
        }
    }
}


                                        # Where should files be saved?
fileloc <- "./doublewells/data/"
                                        # Generate (or load) the model networks in this analysis.
model_networks <- c("erdos_renyi", "er_islands", "barabasi_albert", "LFR", "powerlaw", "fitness")
                                        # Some parameters for the Erdos-Renyi family networks
nnodes <- 100; cprob <- 0.05; nislands <- 5; nbridges <- 1
                                        # The fundamental Erdos-Renyi model
erdos_renyi <- convert(
    sample_gnp(nnodes, cprob)
)
                                        # An igraph adjustment to the ER model, enforcing community
                                        # structure.
er_islands <- convert(
    sample_islands(
        islands.n = nislands, islands.size = nnodes/nislands,
        islands.pin = cprob*nislands*(100 - nbridges)/100,
        n.inter = nbridges
    )
)
                                        # A fundamental Barabasi-Albert model
barabasi_albert <- convert(
    sample_pa(nnodes, directed = FALSE, m = 2, start.graph = make_full_graph(2))
)
                                        # LFR model, using output from a Python/NetworkX model
LFR <- convert(
    graph_from_adjacency_matrix(
        as.matrix(read.table("./data/LFR-NetworkX.txt")),
        mode = "undirected"
    )
)
                                        # A configuration model network with a powerlaw degree
                                        # distribution.
powerlaw <- convert(
    sample_powerlaw(N = nnodes, alpha = 2, kmin = 1, max.iter = 50)
)
                                        # A Goh et al. 2001  node fitness network
fitness <- convert(
    sample_fitness_pl(no.of.nodes = nnodes, no.of.edges = ecount(barabasi_albert), exponent.out = 2)
)
                                        # For looking at the network drawings and some summary
                                        # statistics.
if(check_networks) {
    for(graphname in model_networks) {
        g <- get(graphname)
        dev.new(); plot(g, vertex.size = 3, vertex.label = "", main = graphname)
        print(graphname); print(vcount(g)); print(ecount(g)); print(edge_density(g))
    }
}
                                        # If files are to be saved, do so.
if(save_files) {
    for(i in 1:length(model_networks)) {
        filename <- paste0(fileloc, model_networks[i], ".rda")
        save(list = c(model_networks[i]), file = filename)
    }
}

## Empirical Networks
                                        # A data frame (for convenience) of the chosen empirical
                                        # networks. All networks were taken from the networkdata
                                        # package.
empiricals <- data.frame( 
    dataselect = c(
        "karate", "covert_38", "netsci", "jazz", "covert_17", "hall", "highschool_boys", "surfersb",
        "animal_10", "dolphins_2", "animal_11", "animal_36", "animal_18", "animal_23", "animal_35",
        "animal_6", "animal_33"
    ),
    listno = c(NA, 1, NA, NA, NA, NA, NA, NA, 17, NA, 3, 1, 1, 1, 2, 1, 88),
    pkgname = c("karate", "pira", "netsci", "jazz", "drugusers", "hall", "highschoolboys", "surfers",
        "weaverbirds", "dolphins", "nestbox", "lizards", "bats", "elephantseals", "tortoises",
        "housefinches", "voles"
    )
)
                                        # Load the chosen networks
data(list = empiricals$dataselect, package = "networkdata")
                                        # For each network, choose it, simplify it, and save it under
                                        # a memorable name.
for(i in 1:nrow(empiricals)) {
    g <- get(empiricals$dataselect[i])
    if(is.list(g) & !is.igraph(g)) g <- g[[empiricals$listno[i]]]
    g <- convert(g)
    assign(empiricals$pkgname[i], g)
    if(save_files) {
        filename <- paste0(fileloc, empiricals$pkgname[i], ".rda")
        save(list = c(empiricals$pkgname[i]), file = filename)
    }
}
