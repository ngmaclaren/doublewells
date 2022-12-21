## Code by Neil MacLaren 12/17/2022

library(igraph)
library(doublewells)

save_plots <- FALSE # TRUE

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

## The below definition of an error may be too strict.
## What about, instead of binary, we put a %?

sp <- mod_sp(allnets_lower, nodeset = "sentinels")
## srs <- unlist(lapply(sp, function(x) x$stable_range))
## pes <- unlist(lapply(sp, function(x) x$p_error))
## plot(srs, pes, xlab = "Magnitude of Stable Range", ylab = "Error (% of missed nodes)")


### try error/not vs. number of nodes transitioned; this would be some kind of bar chart
### or, try number of errors for length of stable range



n_tr_bins <- sort(unique((unlist(lapply(sp, function(x) x$n_tr)))))
count_errors <- sapply(n_tr_bins, function(x) {
    sum( unlist( lapply(sp, function(df) df$error[df$n_tr == x]) ) )
})
count_totals <- sapply(n_tr_bins, function(x) {
    sum( unlist( lapply(sp, function(df) length(df$error[df$n_tr == x])) ) )
})
## barplot(
##     count_errors,
##     ##legend.text = c("Error", "No Error"),
##     names.arg = n_tr_bins
## )

## hist(unlist(lapply(sp, function(x) x$n_tr)),
##      breaks = 20, main = "", xlab = "# Nodes in a Transition")

plotthis <- function(sp, filter = NULL, marker = 1) {
    ##sp <- mod_sp(allnets_lower, "sentinels")
    allsp <- do.call(rbind, sp)
    if(is.null(filter)) {
        agg <- aggregate(error ~ n_tr, data = allsp, FUN = sum)
    } else {
        agg <- aggregate(error ~ n_tr, data = allsp[allsp$stable_range > filter, ], FUN = sum)
    }
    points(error ~ n_tr, data = agg, pch = marker, col = marker, cex = 1.25, lwd = 1.25)
    lines(error ~ n_tr, data = agg, col = marker, lty = 1)
}

                                        # This one
ht <- 7; wd <- 7
if(save_plots) {
    pdf("./img/counting-errors2.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(4, 4, 1, 1))
plot(NULL, xlim = c(1, 200), ylim = c(0, 160), xlab = "# Nodes Transitioning", ylab = "# Errors", log = "x")
nodesets <- c("sentinels", "random", "largecorr", "largesd")
markers <- c(3, 6, 7, 8)
for(i in 1:length(nodesets)) plotthis(mod_sp(allnets_lower, nodesets[i]), marker = markers[i])
legend(
    "topright", bty = "n",
    legend = c("High Input", "Random", "Large Corr", "Large SD"),
    pch = markers, col = markers, pt.lwd = 1.5
)
if(save_plots) dev.off()

## plotthis(mod_sp(allnets_lower, "sentinels"), 0.1, 3)
## plotthis(mod_sp(allnets_lower, "random"), 0.1, 6)
## plotthis(mod_sp(allnets_lower, "largecorr"), 0.1, 7)
## plotthis(mod_sp(allnets_lower, "largesd"), 0.1, 8)


res_lower <- calc_results(allnets_lower)
x <- res_lower$n_tr
e_highinput <- res_lower$errors
e_random <- calc_results(allnets_lower, nodeset = "random")$errors
e_largecorr <- calc_results(allnets_lower, nodeset = "largecorr")$errors
e_largesd <- calc_results(allnets_lower, nodeset = "largesd")$errors
y <- list(e_highinput, e_random, e_largecorr, e_largesd)
marker <- c(3, 6, 7, 8)

ht <- 7; wd <- 7
if(save_plots) {
    pdf("./img/counting-errors.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(4, 4, 1, 1), fig = c(0, 1, 0, 1)) # x1, x2, y1, y2
plot(
    NULL, xlab = "Number of Transition Events", ylab = "Number of Errors",# log = "x",
    xlim = range(x), ylim = range(unlist(y))
)
for(i in 1:length(y)) points(x, y[[i]], pch = marker[i], col = marker[i], lwd = 2, cex = 2)
legend(
    "toplef", bty = "n",
    legend = c("High Input", "Random", "Large Corr", "Large SD"),
    pch = marker, col = marker, pt.lwd = 2
)
par(mar = c(1, 1, 1, 1), fig = c(.15, .5, .5, .85), new = TRUE)
plot(
    NULL, xlab = "", ylab = "", xlim = c(0, 5), ylim = c(0, 5), cex.axis = .75#, axes = FALSE
)
## box()
## axis(1, labels = FALSE)
## axis(2, labels = FALSE)
for(i in 1:length(y)) points(x, y[[i]], pch = marker[i], col = marker[i], lwd = 1.5, cex = 1.5)
if(save_plots) dev.off()

## res_upper <- calc_results(allnets_upper, fromlower = FALSE, nodeset)




### Errors are defined as a transition (not necessarily after a stable range) for which not all sentinels transitioned, but more nodes transitioned than sentinels did
## res_lower[order(res_lower$errors, decreasing = TRUE), ]
## res_upper[order(res_upper$errors, decreasing = TRUE), ]



kendalls <- lapply(allnets_lower, function(dl) {
    ##for(i in 1:length(allnets_lower)) {
    df <- dl$df
    rdf <- df[, -c(grep("lowinput", colnames(df)))]
    dat <- as.data.frame(Kendall_correlations(rdf)$means)
    dat$fullname <- rownames(dat)
    dat <- cbind(dat, do.call(rbind, strsplit(dat$fullname, "_")))
    dat <- dat[, -grep("fullname", colnames(dat))]
    colnames(dat) <- c("tau", "signal", "set")
    rownames(dat) <- NULL
    dat
})

dat <- lapply(kendalls, function(dat) {
    dat[dat$signal == "avgac" & dat$set %in% c("sentinel", "random", "largecorr", "largesd"), ]
})
for(i in 1:length(dat)) dat[[i]]$marker <- c(3, 6, 7, 8)
for(i in 1:length(dat)) dat[[i]]$ypos <- 4:1
for(i in 1:length(dat)) dat[[i]]$xpos <- 1:4

tau_highinput <- sapply(dat, function(x) x$tau[x$set == "sentinel"])
tau_random <- sapply(dat, function(x) x$tau[x$set == "random"])
tau_largecorr <- sapply(dat, function(x) x$tau[x$set == "largecorr"])
tau_largesd <- sapply(dat, function(x) x$tau[x$set == "largesd"])
                                        # and this one
ht <- 7; wd <- 7
if(save_plots) {
    pdf("./img/errors-vs-tau.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(4, 4, 1, 1))
plot(NULL, xlab = expression(tau), ylab = "Number of Errors", xlim = c(0, 1), ylim = range(unlist(y)))
points(tau_highinput, e_highinput, pch = 3, col = 3, lwd = 1)
points(tau_random, e_random, pch = 6, col = 6, lwd = 1)
points(tau_largecorr, e_largecorr, pch = 7, col = 7, lwd = 1)
points(tau_largesd, e_largesd, pch = 8, col = 8, lwd = 1)
points(mean(tau_highinput), mean(e_highinput), pch = 3, col = 3, lwd = 3, cex = 4)
points(mean(tau_random), mean(e_random), pch = 6, col = 6, lwd = 3, cex = 4)
points(mean(tau_largecorr), mean(e_largecorr), pch = 7, col = 7, lwd = 3, cex = 4)
points(mean(tau_largesd), mean(e_largesd), pch = 8, col = 8, lwd = 3, cex = 4)
points(tau_highinput[22:23], e_highinput[22:23], pch = 1, col = 2, lwd = 1, cex = 2.5)
points(tau_random[22:23], e_random[22:23], pch = 1, col = 2, lwd = 1, cex = 2.5)
points(tau_largecorr[22:23], e_largecorr[22:23], pch = 1, col = 2, lwd = 1, cex = 2.5)
points(tau_largesd[22:23], e_largesd[22:23], pch = 1, col = 2, lwd = 1, cex = 2.5)
legend(
    "topright", bty = "n",
    legend = c("High Input", "Random", "Large Corr", "Large SD"),
    pch = markers, col = markers, pt.lwd = 1.5
)
if(save_plots) dev.off()

ht <- 7; wd <- 7
dev.new(height = ht, width = wd)
par(mar = c(3, 4, 1, 1))
plot(NULL, xlim = c(.5, 4.5), ylim = c(0, 1), ylab = expression(tau), xlab = "", axes = FALSE)
box()
axis(2)
axis(1, at = 1:4, labels = c("High Input", "Random", "Large Corr", "Large SD"), las = 1)
for(i in 1:length(dat)) {
    points(jitter(dat[[i]]$xpos), dat[[i]]$tau, col = dat[[i]]$marker, pch = dat[[i]]$marker,
           lwd = 2, cex = 2)
}
    
##     dat$signal <- factor(dat$signal, levels = rev(signals))
##     dat$set <- factor(dat$set, levels = nodesets)
##     dat$signal_ <- as.numeric(dat$signal) + adjust[i]
##     dat$set_ <- as.numeric(dat$set)
##     points(signal_ ~ tau, data = dat, col = dat$set_, pch = dat$set_, lwd = pt.lwd, cex = pt.cex)
##     print(summary(dat$tau[dat$set %in% c("all", "lower", "sentinel")]))
## }
