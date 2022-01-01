library(igraph)
source("double-well-functions.R")

nnodes <- 500
ndegree <- 8

rrg <- sample_k_regular(nnodes, ndegree)
rrA <- as_adj(rrg, type = "both", sparse = FALSE)
## plot(rrg, vertex.size = 2, vertex.label = NA)

sfg <- sample_fitness_pl(nnodes, nnodes * ndegree, 3.5)
sfA <- as_adj(sfg, type = "both", sparse = FALSE)
## plot(sfg, vertex.size = 2, vertex.label = NA)

r1 <- 1
r2 <- 2
r3 <- 5
D <- 0.1
dt <- 0.001
T <- 4000

dev.new()
plot(NULL,  xlab = expression(t), ylab = expression(x[i]), xlim = c(0, T), ylim = c(0, 10))
abline(h = 1, col = "lightgray")
abline(h = 2, col = "lightgray")
abline(h = 5, col = "lightgray")
text(T, y = 1, labels = expression(r[1]), adj = c(1, 1), col = "lightgray")
text(T, y = 2, labels = expression(r[2]), adj = c(1, 1), col = "lightgray")
text(T, y = 5, labels = expression(r[3]), adj = c(1, 1), col = "lightgray")
colors <- c("#db7050", "#daa72d")
legend("topright", bty = "n", legend = c("Random Regular", "Scale-Free"), col = colors, pch = 1)

##color1 <- "#db7050" ## random regular
##color2 <- "#daa72d" ## scale-free

## x_i == 0.01, 10 ∀ i
As <- list(rrA, sfA)
xs <- list(rep(0.01, nnodes), rep(10, nnodes))
for(i in 1:2) {
    A <- As[[i]]
    color <- colors[i]
    for(j in 1:2) {
        x <- xs[[j]]
        sto <- matrix(0, ncol = length(x), nrow = T)

        for(t in 1:T) {
            sto[t, ] <- x
            x <- double_well_coupled(x, r1, r2, r3, D, A, dt)
        }

        points(1:T, rowMeans(sto), col = color, pch = 1)
    }
}

beta_eff <- function(g) {
    s.in <- degree(g, mode = "in")
    s.out <- degree(g, mode = "out")

    mean(s.in * s.out)/mean(s.out)
}

beta_eff(rrg) # 8
beta_eff(sfg) # 22.1825

### How do I adjust β_eff on the graphs themselves?
rredges <- sample(E(rrg), 1500)
rrg2 <- delete_edges(rrg, rredges)
beta_eff(rrg2)
rrA2 <- as_adj(rrg2, type = "both", sparse = FALSE)

sfedges <- sample(E(sfg), 3500)
sfg2 <- delete_edges(sfg, sfedges)
beta_eff(sfg2)
sfA2 <- as_adj(sfg2, type = "both", sparse = FALSE)

### deleting edges works. So I think that's how they did the sim.
### these are kind of like network death bifurcation diagrams
### but shouldn't there be lots of ways to get the same β_eff?

dev.new()
plot(NULL,  xlab = expression(t), ylab = expression(x[i]), xlim = c(0, T), ylim = c(0, 10))
abline(h = 1, col = "lightgray")
abline(h = 2, col = "lightgray")
abline(h = 5, col = "lightgray")
text(T, y = 1, labels = expression(r[1]), adj = c(1, 1), col = "lightgray")
text(T, y = 2, labels = expression(r[2]), adj = c(1, 1), col = "lightgray")
text(T, y = 5, labels = expression(r[3]), adj = c(1, 1), col = "lightgray")
colors <- c("#db7050", "#daa72d")
legend("topright", bty = "n", legend = c("Random Regular", "Scale-Free"), col = colors, pch = 1)

##color1 <- "#db7050" ## random regular
##color2 <- "#daa72d" ## scale-free

## x_i == 0.01, 10 ∀ i
As <- list(rrA2, sfA2)
xs <- list(rep(0.01, nnodes), rep(10, nnodes))
for(i in 1:2) {
    A <- As[[i]]
    color <- colors[i]
    for(j in 1:2) {
        x <- xs[[j]]
        sto <- matrix(0, ncol = length(x), nrow = T)

        for(t in 1:T) {
            sto[t, ] <- x
            x <- double_well_coupled(x, r1, r2, r3, D, A, dt)
        }

        points(1:T, rowMeans(sto), col = color, pch = 1)
    }
}
