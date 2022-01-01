library(igraph)
source("double-well-functions.R")

## Two Node
## A <- matrix(
##     c(0, 1,
##       1, 0),
##     byrow = TRUE, nrow = 2
## )

## Three Node
A <- matrix(
    c(0, 1, 1,
      0, 0, 1,
      0, 0, 0),
    byrow = TRUE, nrow = 3
)

## Four Node
## A <- matrix(
##     c(0, 1, 1, 1,
##       1, 0, 1, 1,
##       1, 1, 0, 1,
##       1, 1, 1, 0),
##     byrow = TRUE, nrow = 4
## )

check_graph_structure <- FALSE

if(check_graph_structure) {
    dev.new()
    g <- graph_from_adjacency_matrix(A, mode = "directed")
    plot(g)
}

dt <- 0.001
r1 <- 1
r2 <- 2
r3 <- 5
D <- 0.1
noise_level <- 10
initialx <- rep(2, nrow(A))

x <- initialx
T <- 2000
xs <- matrix(0, ncol = length(x), nrow = T)

for(t in 1:T) {
    xs[t, ] <- x
    x <- double_well_coupled_noisy(x, r1, r2, r3, D, A, dt, noise = rnorm(length(x), 0, noise_level))
}

pal <- list(# underworld_palette
    "Ghostpipe" = "#eff8f6", "WhiteToadstool" = "#e4e6d5", "Cupshroom" = "#d0c597",
    "Lichen" = "#95a696", "ShroomBloom" = "#5a5046", "FeatherFern" = "#93c54a",
    "YellowToadstool" = "#e9d466", "FeatherShroom" = "#daa72d", "FalseChant" = "#ec9813",
    "CrackingTable" = "#d67534", "TinyToadstool" = "#db7050", "RedButton" = "#d15935",
    "PurpleButton" = "#805350"
)
##colors <- unlist(rev(pal)[1:length(x)])#c(pal$Lichen, pal$Cupshroom)
shift <- 2
colors <- unlist(pal[seq(1 + shift, length(x)*2 + shift, by = 2)])

dev.new(width = 13, height = 7)
plot(NULL, xlim = c(0, T), ylim = c(0, 6),
     xlab = expression(t), ylab = expression(x[i]), main = "Coupled Double Wells")
legend("bottomright", bty = "n",
       ##legend = sapply(1:length(x), paste0("x_", i)),
       legend = 1:length(x), title = expression(x[i]),
       pch = 19, col = colors)
for(i in 1:length(x)) points(1:T, xs[, i], col = colors[i], pch = 1)
