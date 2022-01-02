## Code by Neil G. MacLaren, October, 2021
## This file shows the effect of varying the connection strength parameter, D, on selected triads.

library(igraph)
library(doublewells)

##set.seed(12345)

matrices <- list(
    complete = complete,
    feedback = feedback,
    feedforward = feedforward,
    fromthemiddleout = fromthemiddleout,
    linemotif = linemotif,
    tothemiddle = tothemiddle
)

r1 <- 1
r2 <- 2
r3 <- 5

dt <- 0.001
T <- 5000

initialx <- rep(r1, 3)
noise_level <- 10
noise <- function(n, noise_level) rnorm(n = n, mean = 0, sd = noise_level)
Ds <- seq(.1, 1, by = .1)

results <- vector(mode = "list", length = length(matrices))
names(results) <- names(matrices)
results <- lapply(
    results, function(m) {
        m <- vector(mode = "list", length = length(Ds))
        names(m) <- Ds
        lapply(m, function(n) n <- matrix(0, ncol = length(initialx), nrow = T))
    }
)

for(i in 1:length(matrices)) {
    A <- matrices[[i]]

    for(j in 1:length(Ds)) {
        D <- Ds[[j]]

        x <- initialx
        for(t in 1:T) {
            results[[i]][[j]][t, ] <- x
            x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise(length(x), noise_level))
        }
    }
}

### Plotting
pal <- rainbow(10)

dev.new(width = 20, height = 15)
##pdf("./img/varying-D.pdf", width = 20, height = 15)

nrows <- length(matrices)
layoutmatrix <- matrix(
    c(1, 2, 3, 7, 8, 9, 7, 8, 9, 4, 5, 6, 10, 11, 12, 10, 11, 12),
    byrow = FALSE, ncol = 6
)
figlayout <- layout(layoutmatrix)
##layout.show(figlayout)

linecolor <- "black"
fillcolor <- "gainsboro"
for(i in 1:nrows) {
    A <- matrices[[i]]
    g <- graph_from_adjacency_matrix(A, mode = "directed")
    plot(
        g,
        vertex.color = fillcolor,
        vertex.frame.color = linecolor,
        vertex.label.color = linecolor,
        vertex.size = 25,
        vertex.label.cex = 2,
        edge.color = linecolor
    )
}

for(i in 1:nrows) {
    plot(
        NULL,
        xlim = c(0, T), ylim = c(r1*.75, r3*1.25),
        main = "",
        xlab = expression(t), ylab = expression(x[i])
    )

    if(i == 4) legend("topleft", legend = as.character(Ds), lty = 1, lwd = 4, col = pal)
    if(i == 5) legend("topleft",
                      legend = c(expression(x[1]), expression(x[2]), expression(x[3])),
                      lty = 1:3, lwd = 3, col = pal[1])

    for(j in 1:length(Ds)) {
        color <- pal[j]
        matlines(results[[i]][[j]], col = color, lwd = 2)
    }
}
##dev.off()
