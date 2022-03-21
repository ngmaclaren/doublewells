library(igraph)
##library(networkdata)

## NB: sn_auth, netsci, and others have weights, but downstream analysis will ignore the weights. Other networks in this list are also formed from data that originally could be considered weighted, but are dichotomized in the network.

set.seed(123)
save_files <- TRUE # FALSE

## Testing and Collection Functions
get_gcc <- function(g) {
    require(igraph)
    comps <- components(g)
    gcc_id <- which.max(comps$csize)
    vids <- V(g)[comps$membership == gcc_id]
    g <- induced_subgraph(g, vids)
    g
}

convert <- function(g) as.undirected(simplify(get_gcc(g)))

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

check_sizes <- function(graphlist) {
    sapply(graphlist, function(g) vcount(get_gcc(g)))
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
    assign(graphs[i], get_gcc(generate_network(graphs[i])))
    filename <- paste0(fileloc, graphs[i], ".rda")
    if(save_files) save(list = c(graphs[i]), file = filename)
}

## LFR model, using output from a Python/NetworkX model
LFR <- as.matrix(read.table("./data/LFR-NetworkX.txt"))
LFR <- get_gcc(graph_from_adjacency_matrix(LFR, mode = "undirected"))
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
    assign(graphnames[i], get_gcc(get(graphnames[i])))
    if(save_files) save(list = graphnames[i], file = outfilenames[i])
}

## Empirical Networks

## Thursday, March 3, 2022: take out jpr and sn_auth, probably the small ones, too?
## Add: animal_11[[1]], animal_36[[1]], covert_17 (GCC)
## Maybe: animal_18[[1]], animal_23[[1]] (GCC), animal_33[[21]] (GCC), animal_35[[4]], animal_6[[1]], hall (as.undirected), highschool_boys (as.undirected, simplify), surfersb (simplify)
empiricals <- c(
    ## Used for Prosenjit's paper but not EWI
    "mine", "jpr", "sn_auth", "taro", "wiring", "covert_27",
    ## EWI; these all have a reference
    "covert_38", "karate", "jazz", "animal_10", "dolphins_2", "netsci",
    "animal_11", "animal_36", "covert_17", "animal_18", "animal_23", "animal_33", "animal_35",
    "animal_6", "hall", "highschool_boys", "surfersb"
)
                                        # choose from animal_10 at random
                                        # choose [[1]] from covert_27

data(list = empiricals, package = "networkdata")

## Choose network at random from list objects
## animal_10

## Get GCCs
sn_auth <- convert(sn_auth)
##petster <- get_gcc(petster) # too big
pira <- convert(covert_38[[1]]) # close to 100, less tree-like than [[2]]
weaverbirds <- convert(animal_10[[17]]) # the largest network from the list, was 4
linux <- convert(covert_27[[1]]) # the smaller of the two networks
dolphins <- convert(dolphins_2)
netsci <- convert(netsci)
empiricals[which(empiricals == "covert_38")] <- "pira"
empiricals[which(empiricals == "animal_10")] <- "weaverbirds"
empiricals[which(empiricals == "covert_27")] <- "linux"
empiricals[which(empiricals == "dolphins_2")] <- "dolphins"

nestbox <- convert(animal_11[[3]]) # has more bridging connections, near 100
empiricals[which(empiricals == "animal_11")] <- "nestbox"
lizards <- convert(animal_36[[1]])
empiricals[which(empiricals == "animal_36")] <- "lizards"
drugusers <- convert(covert_17)
empiricals[which(empiricals == "covert_17")] <- "drugusers"
bats <- convert(animal_18[[1]])
empiricals[which(empiricals == "animal_18")] <- "bats"
elephantseals <- convert(animal_23[[1]])
empiricals[which(empiricals == "animal_23")] <- "elephantseals"
tortoises <- convert(animal_35[[2]])
empiricals[which(empiricals == "animal_35")] <- "tortoises"
housefinches <- convert(animal_6[[1]])
empiricals[which(empiricals == "animal_6")] <- "housefinches"
voles <- convert(animal_33[[88]])
empiricals[which(empiricals == "animal_33")] <- "voles"
hall <- convert(hall)
highschool_boys <- convert(highschool_boys)
surfersb <- convert(surfersb)


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


## for(var in empiricals) {
##     obj <- get(var)
##     if(class(obj) == "igraph") {
##         dev.new()
##         checkplot(obj)
##     }
## }
