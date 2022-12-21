## Code by Neil MacLaren 12/17/2022

library(igraph)
library(doublewells)

save_plots <- FALSE # TRUE
palette("R4")

networks <- c(
    "erdos_renyi", "er_islands", "barabasi_albert", "LFR", "powerlaw", "fitness",
    empiricals
)
networks <- networks[-which(networks %in% c("powerlaw", "dolphins"))]

load("./data/allnets-lower-r1.rda")

names(allnets_lower) <- names(allnets_upper) <- networks

load("./data/examples-lower-r1.rda")

allnets_lower[["powerlaw"]] <- examples_lower[[1]]
allnets_lower[["dolphins"]] <- examples_lower[[2]]
networks <- c(networks, "powerlaw", "dolphins")
data(list = networks)
samples <- rev(seq(from = 75/.01, by = -.1/.01, length.out = 240))

all_cvs <- lapply(allnets_lower, function(dl) rev(find_critical_values(dl)))

## applied to the above list, conduct the simulations at just those values of D, returning all histories (including X, which is what I need).
## Requires modification to `simulation` to allow passing specified values of D.

simresults <- list()
for(i in seq_along(all_cvs)) {
    g <- get(networks[i])
    simresults[[i]] <- simulation(
        g, check_alts = TRUE, return_histories = TRUE, assessment_samples = 1:10,
        return_X = TRUE,
        specific_values = TRUE, Dvec = all_cvs[[i]]
    )
}
names(simresults) <- networks

## need: pl_kcs, pl_rs
## need a big palette (with at least 23 colors)
kcs <- list()
rs <- list()
for(i in seq_along(all_cvs)) {
    g <- get(networks[i])
    k <- degree(g)
    cvs <- all_cvs[[i]]
    kcs[[i]] <- Kendall_correlations(allnets_lower[[i]]$df)$kendalls[, "avgac_largecorr"]
    rs[[i]] <- sapply(cvs, function(D) cor(k, est_degree(simresults[[i]], D), method = "spearman"))
}
names(kcs) <- networks
names(rs) <- networks

## cor(unlist(kcs), unlist(rs))
## model <- lm(unlist(rs) ~ unlist(kcs))

library(nlme)

y <- unlist(kcs)
x <- unlist(rs)
vm <- gls(y ~ x, weights = varFixed(~ x))
## plot(vm) # evidence for changing variance with x is not as strong with rs as x var; diagnostic plot suggests this model is not appropriate
lsm <- lm(y ~ x)
## plot(lsm) # diagnostics look generally ok. Decreased sample size at larger values of x is potentially behind the apparent decrease in variance with x.

ht <- 7; wd <- 7
if(save_plots) {
    cairo_pdf("./img/est-degree-correlations2.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(4, 4, 1, 1))
plot(NULL, xlim = c(-.25, 1), ylim = c(-.25, 1), ylab = "Large Corr (τ)", xlab = "ρ(k, k')")
abline(v = 0, lwd = .5); abline(h = 0, lwd = .5)
for(i in seq_along(all_cvs[-which(names(all_cvs) %in% c("powerlaw", "dolphins"))]))
    points(rs[[i]], kcs[[i]], col = 1, pch = 19, cex = 1.5)
points(rs[["powerlaw"]], kcs[["powerlaw"]], pch = 19, col = 3, cex = 1.5)
points(rs[["dolphins"]], kcs[["dolphins"]], pch = 19, col = 4, cex = 1.5)
## lines(unlist(kcs), predict(model), col = "black")
lines(x, predict(lsm), col = 2, lwd = 2)
## lines(x, predict(vm), col = 2, lwd = 2)
text(
    .6, .1, col = 2,
    paste0("y = ", round(coefficients(lsm)[1], 3), " + ", round(coefficients(lsm)[2], 3), "x")
)
legend(
    "bottomrigh", bty = "n", col = c(3, 4, 1), pch = 19,
    legend = c("Power-Law", "Dolphins", "Other Networks")
)
if(save_plots) dev.off()
