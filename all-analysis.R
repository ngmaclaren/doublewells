## Code by Neil MacLaren 12/17/2022

library(igraph)
library(doublewells)

save_plots <- TRUE # FALSE

networks <- c(
    "erdos_renyi", "er_islands", "barabasi_albert", "LFR", "powerlaw", "fitness",
    empiricals
)
networks <- networks[-which(networks %in% c("powerlaw", "dolphins"))]

load("./data/allnets-lower-r1.rda")
load("./data/allnets-upper-r1.rda")

names(allnets_lower) <- names(allnets_upper) <- networks

load("./data/examples-lower-r1.rda")
load("./data/examples-upper-r1.rda")

allnets_lower[["powerlaw"]] <- examples_lower[[1]]
allnets_lower[["dolphins"]] <- examples_lower[[2]]
allnets_upper[["powerlaw"]] <- examples_upper[[1]]
allnets_upper[["dolphins"]] <- examples_upper[[2]]
networks <- c(networks, "powerlaw", "dolphins")
data(list = networks)

plotthis <- function(agg, marker = 1) {
    points(
        agg[, "n_tr"], agg[, "error"],
        pch = marker, col = marker, cex = 1.25, lwd = 1.25
    )
    lines(
        agg[, "n_tr"], agg[, "error"],
        col = marker, lty = 1, lwd = 1.25
    )
}

aggthis <- function(sp) {
    allsp <- do.call(rbind, sp)
    n_tr_bins <- sort(unique(allsp$n_tr))
    count <- sapply(n_tr_bins, function(x) sum(allsp$n_tr == x))
    errors <- sapply(n_tr_bins, function(x) sum(allsp$error[allsp$n_tr == x]))
    p_errors <- errors/count
    df <- data.frame(n_tr = n_tr_bins, error = errors, count = count)
    return(df)
}
    
nodesets <- c("sentinels", "random", "largecorr", "largesd")
plotdata <- lapply(nodesets, function(x) mod_sp(allnets_lower, x))
names(plotdata) <- nodesets
for(i in seq_along(plotdata)) plotdata[[i]] <- aggthis(plotdata[[i]])
yrange <- range(unlist(lapply(plotdata, function(x) x$count)))
xrange <- range(unlist(lapply(plotdata, function(x) x$n_tr)))

ht <- 7; wd <- 7
if(save_plots) {
    cairo_pdf(
        "./img/counting-errors2.pdf", height = ht, width = wd, family = "Bitsream Vera Sans"
    )
} else dev.new(height = ht, width = wd)
par(mar = c(4, 4, 1, 1) + .4)
plot(
    NULL,
    xlim = xrange, ylim = yrange,
    xlab = "# Nodes Transitioning", ylab = "# Errors", log = "x",
    cex.lab = 1.75, cex.axis = 1.75
)
markers <- c(3, 6, 7, 8)
for(i in seq_along(plotdata)) plotthis(plotdata[[i]], marker = markers[i])
points(plotdata[[1]]$n_tr, plotdata[[1]]$count, col = 5, pch = 1, cex = 1.25, lwd = 1.25)
lines(plotdata[[1]]$n_tr, plotdata[[1]]$count, col = 5, lty = 1, lwd = 1.25)
legend(
    "topright", bty = "n",
    legend = c("High Input", "Random", "Large Corr", "Large SD", "# Transitions"),
    pch = c(markers, 1), col = c(markers, 5),
    pt.lwd = 1.5, cex = 1.25
)
if(save_plots) dev.off()
