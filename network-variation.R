## Need to add another layer to the multilevel model: `run`
## There should be 30--50 runs, the results of which should each be saved to a separate .rda
## There should be a subdirectory under data for all those .rda files
## I also want to split this file into two:
### network-variation-sims.R that conducts the sims and saves the results
### network-variation-analysis.R that loads the results, runs the models, and makes the figure

## Code by Neil MacLaren 3/1/2022
library(nlme)
library(igraph)
library(doublewells)

set.seed(123) # v1: 123, v2: 1234, v3: 12345

run_sims <- TRUE # FALSE
save_plots <- FALSE # TRUE
save_results <- TRUE # FALSE

networks <- c(
    "max_entropy", "me_islands", "pref_attach", "scale_free_2", "LFR",
    empiricals
)

if(run_sims) {
    data(list = networks)

    results <- list()
    for(i in 1:length(networks)) {
        net <- networks[i]
        print(net)
        g <- get(net)
        results[[i]] <- simulation(g, TU = 125, check_equil = TRUE)
    }

    network_variation_results <- results
    save(network_variation_results, file = "./data/network-variation-results-test.rda")
} else {
    load("./data/network-variation-results.rda")
    ##load("./data/network-variation-results-v1.rda")
    ##load("./data/network-variation-results-v2.rda")
    ##load("./data/network-variation-results-v3.rda")
    results <- network_variation_results
    names(results) <- networks
}

corr_results <- lapply(results, function(x) Kendall_correlations(x[[1]])$means)
df <- as.data.frame(t(do.call(cbind, corr_results)))
varcols <- colnames(df)
df$network <- networks

rdf <- reshape(
    df,
    varying = varcols,
    v.names = "tau",
    timevar = "ewi_",
    times = varcols,
    new.row.names = 1:10000,
    direction = "long"
)

newcols <-  strsplit(rdf$ewi_, "_")
newcols <- as.data.frame(do.call(rbind, newcols))
colnames(newcols) <- c("ewi", "nodeset")
rdf <- cbind(rdf, newcols)
rdf$nodeset <- factor(rdf$nodeset, levels = c("all", "lower", "sentinel"), ordered = FALSE)

## This is the anova way
##m1 <- aov(tau ~ ewi_ + Error(network/ewi_), data = rdf[grep("avgac", rdf$ewi_), ])

## The mixed effects way is more useful as we get a coefficient for eeach level of ewi
## ~ 1 | network/nodeset means the intercept is random over nodesets within networks
avgsd <- lme(tau ~ nodeset, random = ~ 1 | network/nodeset, data = rdf[rdf$ewi == "avgsd", ])
avgac <- lme(tau ~ nodeset, random = ~ 1 | network/nodeset, data = rdf[rdf$ewi == "avgac", ])

if(save_results) {
    if(!dir.exists("./r-out")) dir.create("./r-out")
    sink("./r-out/network-variation-test.txt", append = FALSE, type = "output", split = TRUE)

    print(summary(avgsd))
    print(summary(avgac))

    sink()
}


get_values <- function(model) {
    ints <- intervals(model, which = "fixed")
    ests <- ints$fixed[, 2]
    ests[2:3] <- ests[1] + ests[2:3]
    lowers <- ints$fixed[, 1]
    lowers[2:3] <- ests[1] + lowers[2:3]
    uppers <- ints$fixed[, 3]
    uppers[2:3] <- ests[1] + uppers[2:3]

    df <- as.data.frame(t(rbind(lowers, ests, uppers)))
    ## df <- df[c(2, 1, 3), ] # needed if lower is the ref in the model
    rownames(df) <- c("all", "lower", "sentinel")
    df
}

plotvals <- list(
    avgsd = get_values(avgsd),
    avgac = get_values(avgac)
)

colors <- c(
    "sienna", # all
    "navy", # lower
    "olivedrab" # sentinel
)

yticklabels <- c("Average Standard Deviation", "Average Autocorrelation")
legendlabels <- c("All Nodes", "Lower State Nodes", "Sentinel Nodes")
ht <- 7
wd <- 8
if(save_plots) {
    pdf("./img/kendall-corr-figure.pdf", width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mar = c(5, 12, 1, 1) + 0.1)
ylims <- c(.5, 2.5)
plot(NULL, ylim = ylims, xlim = c(0.45, .95), #xlim = c(min(lowers)*.9, max(uppers)*1.1),
     ylab = "", xlab = expression(tau), yaxt = "n")
axis(side = 2, tick = FALSE, labels = yticklabels, at = 2:1, las = 1)
ypos <- 2 + c((1/3), 0, -(1/3))
points(plotvals[[1]]$ests, ypos, pch = 16:18, col = colors, cex = 2)
segments(x0 = plotvals[[1]]$lowers, x1 = plotvals[[1]]$uppers, y0 = ypos,
         lty = 1, lwd = 2, col = colors)
abline(h = 1.5, col = "gray", lwd = .5, lty = 1)
ypos <- 1 + c((1/3), 0, -(1/3))
points(plotvals[[2]]$ests, ypos, pch = 16:18, col = colors, cex = 2)
segments(x0 = plotvals[[2]]$lowers, x1 = plotvals[[2]]$uppers, y0 = ypos,
         lty = 1, lwd = 2, col = colors)
legend("bottomleft", legend = legendlabels,
       col = colors, bty = "n", pch = 16:18, lty = 1, lwd = 1)#, pt.cex = 2)
if(save_plots) dev.off()




names(results) <- networks
for(i in 1:length(results)) {
    dl <- results[[i]]
    dev.new()
    plot(NULL, xlim = c(0, 11), ylim = c(0, 11), xlab = "At Equil", ylab = "At End",
         main = networks[i])
    for(df in dl$equil) points(df)
}
