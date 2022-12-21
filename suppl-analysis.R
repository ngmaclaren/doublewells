library(igraph)
library(doublewells)

save_plots <- FALSE # TRUE

load("./data/examples-lower-r1-withX.rda")
data(list = c("powerlaw", "dolphins"))

pl <- examples_lower[[1]]
dol <- examples_lower[[2]]

## pl_A <- as_adj(powerlaw, type = "both", sparse = FALSE)
## dol_A <- as_adj(dolphins, type = "both", sparse = FALSE)

pl_cv <- rev(find_critical_values(pl))
dol_cv <- rev(find_critical_values(dol))

pl_k <- degree(powerlaw)
dol_k <- degree(dolphins)

samples <- rev(seq(from = 75/.01, by = -.1/.01, length.out = 240))

pl_r <- est_degree(pl, pl_cv[1])
dol_r <- est_degree(dol, dol_cv[1])
cor(pl_k, pl_r, method = "spearman") # "pearson"
cor(dol_k, dol_r, method = "spearman") # "pearson"

pl_rs <- sapply(pl_cv, function(D) cor(pl_k, est_degree(pl, D), method = "spearman"))
dol_rs <- sapply(dol_cv, function(D) cor(dol_k, est_degree(dol, D), method = "spearman"))

pl_kcs <- Kendall_correlations(pl$df)$kendalls[, "avgac_largecorr"]
dol_kcs <- Kendall_correlations(dol$df)$kendalls[, "avgac_largecorr"]

pl_sizes <- find_n_nodes_tr(pl_cv, pl)
dol_sizes <- find_n_nodes_tr(dol_cv, dol)

ht <- 7; wd <- 7
if(save_plots) {
    cairo_pdf("./img/est-degree-correlations.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(4, 4, 1, 1))
plot(NULL, xlim = c(-.5, 1), ylim = c(-.5, 1), xlab = "Large Corr (τ)", ylab = "ρ(k, k')")
points(pl_kcs, pl_rs, col = 2, pch = 19, lwd = 2, cex = sqrt(pl_sizes))
## text(pl_kcs, pl_rs, adj = 2)
points(dol_kcs, dol_rs, col = 3, pch = 19, lwd = 2, cex = sqrt(dol_sizes))
## text(dol_kcs, dol_rs, adj = 2)
abline(v = 0, lty = 1, lwd = .5)
abline(h = 0, lty = 1, lwd = .5)
legend(
    "bottomleft", bty = "n", col = c(2, 3), pch = 19,
    legend = c("Power-Law", "Dolphins")
)
if(save_plots) dev.off()

ht <- 7; wd <- 14
dev.new(height = ht, width = wd)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
pl_tdf <- data.frame(
    k = pl_k,
    est_k = est_degree(pl, pl_cv[2]),
    state = colMeans(pl$histories$X[[which(pl$df$Ds == pl_cv[2])]][samples, ])
)
pl_tdf$color <- ifelse(pl_tdf$state > 5, "orange", "blue")
plot(est_k ~ k, data = pl_tdf, type = "p", pch = 1, lwd = 1, cex = 1, col = pl_tdf$color,
     main = "Power-Law", xlab = "Degree", ylab = "Estimated Degree")
dol_tdf <- data.frame(
    k = dol_k,
    est_k = est_degree(dol, dol_cv[1]),
    state = colMeans(dol$histories$X[[which(dol$df$Ds == dol_cv[1])]][samples, ])
)
dol_tdf$color <- ifelse(dol_tdf$state > 5, "orange", "blue")
plot(est_k ~ k, data = dol_tdf, type = "p", pch = 1, lwd = 1, cex = 1, col = dol_tdf$color,
     main = "Power-Law", xlab = "Degree", ylab = "Estimated Degree")
