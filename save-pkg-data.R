library(igraph)
library(networkdata)

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

## Random Graphs

graphs <- c("random_regular", "max_entropy", "sphere_surface",
            "me_islands", "pref_attach", "small_world")

fileloc <- "./doublewells/data/"

for(i in 1:length(graphs)) {
    assign(graphs[i], generate_network(graphs[i]))
    filename <- paste0(fileloc, graphs[i], ".rda")
    save(list = c(graphs[i]), file = filename)
}

## Empirical Networks

empiricals <- c("mine", "jpr", "surfersb")
data(list = empiricals)

for(i in 1:length(empiricals)) {
    filename <- paste0(fileloc, empiricals[i], ".rda")
    save(list = c(empiricals[i]), file = filename)
}


