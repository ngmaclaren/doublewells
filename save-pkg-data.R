library(igraph)
##library(networkdata)

## NB: sn_auth and netsci have weights, but downstream analysis will ignore the weights. Other networks in this list are also formed from data that originally could be considered weighted, but are dichotomized in the network.

set.seed(123)
save_files <- TRUE

## Testing and Collection Functions
get_gcc <- function(g) {
    require(igraph)
    comps <- components(g)
    gcc_id <- which.max(comps$csize)
    vids <- V(g)[comps$membership == gcc_id]
    g <- induced_subgraph(g, vids)
    g
}

CVk <- function(g) {
    k <- degree(g)
    sd(k)/mean(k)
}

checks <- function(g) {
    cat("Directed: ", is_directed(g), "\n",
        "Weighted: ", is_weighted(g), "\n",
        "Simple: ", is_simple(g), "\n",
        "Bipartite: ", is_bipartite(g), "\n",
        "N Nodes: ", vcount(g), "\n",
        "N Edges: ", ecount(g), "\n",
        "CV_k: ", CVk(g), "\n",
        "Size of GCC: ", vcount(get_gcc(g)), "\n")
}

checkplot <- function(g) plot(g, vertex.label = "", vertex.size = 3)

## Random Networks
generate_network <-
    function(choice = c("random_regular", "max_entropy", "sphere_surface",
                        "me_islands", "pref_attach", "small_world"),
             nnodes = 100, cprob = .06,
             rr.k = 6,
             sph.dim = 3, sph.radius = .279,
             nislands = 5, nbridges = 1,
             pa.power = 1.5,
             pa.outdist = c(0, 10, (8+1/3), (6+2/3), 5, (3+1/3), (1+2/3), (5/6))*cprob,
             sw.dim = 1, sw.nei = 3, sw.p = .1)
{
    require(igraph)
    
    choice <- match.arg(choice)
    
    switch(choice,
           random_regular = sample_k_regular(nnodes, rr.k),
           max_entropy = sample_gnp(nnodes, cprob),
           sphere_surface = sample_dot_product(
               sample_sphere_surface(dim = sph.dim, n = nnodes, radius = sph.radius)),
           me_islands = sample_islands(
               islands.n = nislands, islands.size = nnodes/nislands,
               islands.pin = cprob*nislands * (100 - nbridges)/100,
               n.inter = nbridges),
           pref_attach = sample_pa(nnodes, power = pa.power, out.dist = pa.outdist, directed = FALSE),
           small_world = sample_smallworld(dim = sw.dim, size = nnodes, nei = sw.nei, p = sw.p)
           )
}

fileloc <- "./doublewells/data/"

### Random Graphs

## Standard igraph models
graphs <- c("random_regular", "max_entropy", "sphere_surface",
            "me_islands", "pref_attach", "small_world")

for(i in 1:length(graphs)) {
    assign(graphs[i], generate_network(graphs[i]))
    filename <- paste0(fileloc, graphs[i], ".rda")
    if(save_files) save(list = c(graphs[i]), file = filename)
}

## LFR model, using output from a Python/NetworkX model
LFR <- as.matrix(read.table("./data/LFR-NetworkX.txt"))
LFR <- graph_from_adjacency_matrix(LFR, mode = "undirected")
if(save_files) save(LFR, file = paste0(fileloc, "LFR.rda"))

## A sequence of node fitness networks with different power law exponents
nnodes <- vcount(pref_attach)
nedges <- ecount(pref_attach)
exponents <- c(2, 2.1, 2.2, 2.4, 2.6, 2.8, 3, 5, 10, 100)
##graphlist <- list()
graphnames <- paste0("scale_free_", exponents)
outfilenames <- paste0(fileloc, graphnames, ".rda")
for(i in 1:length(exponents)) {
    assign(
        graphnames[i], sample_fitness_pl(no.of.nodes = nnodes, no.of.edges = nedges,
                                         exponent.out = exponents[i])
    )
    if(save_files) save(list = graphnames[i], file = outfilenames[i])
}

## Empirical Networks

empiricals <- c("mine", "jpr", "sn_auth", "covert_38", "karate", "jazz",
                "animal_10", "covert_27", "taro", "wiring",
                "dolphins_2", "netsci")
                                        # choose from animal_10 at random
                                        # choose [[1]] from covert_27

data(list = empiricals, package = "networkdata")

## Choose network at random from list objects
## animal_10

## Get GCCs
sn_auth <- get_gcc(sn_auth)
##petster <- get_gcc(petster) # too big
pira <- get_gcc(covert_38[[1]])
weaverbirds <- get_gcc(animal_10[[4]]) # a larger network from the list
linux <- get_gcc(covert_27[[1]]) # the smaller of the two networks
dolphins <- get_gcc(dolphins_2)
netsci <- get_gcc(netsci)
empiricals[which(empiricals == "covert_38")] <- "pira"
empiricals[which(empiricals == "animal_10")] <- "weaverbirds"
empiricals[which(empiricals == "covert_27")] <- "linux"
empiricals[which(empiricals == "dolphins_2")] <- "dolphins"

for(i in 1:length(empiricals)) {
    filename <- paste0(fileloc, empiricals[i], ".rda")
    if(save_files) save(list = c(empiricals[i]), file = filename)
}

## other_empirical <- c("dolphins", "netscience-lcc")
## for(choice in other_empirical) {
##     infilename <- paste0("./data/", choice, ".mat")
##     edgelist <- as.matrix(read.table(file = infilename, header = FALSE, sep = " ", skip = 1)[, 1:2])
##     if(choice == "netscience-lcc") choice <- "NNS"
##     assign(choice, graph_from_edgelist(edgelist, directed = FALSE))
##     outfilename <- paste0(fileloc, choice, ".rda")
##     save(list = c(choice), file = outfilename)
## }
