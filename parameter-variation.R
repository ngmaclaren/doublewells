## Code by Neil MacLaren 3/1/2022

library(igraph)
library(doublewells)

save_plots <- FALSE
save_results <- TRUE # FALSE
run_sims <- TRUE #  FALSE

which_network <- "dolphins" # "pref_attach" # 
filename <- paste0("./data/parameter-variation-results-", which_network, ".rda")
outfile <- paste0("./r-out/parameter-variation-", which_network, ".txt")

intensities <- c(.01, .05, .1, .5)
sample_sizes <- rev(c(250, 150, 50, 25))
rs <- list(c(1, 3, 5), c(1, 2.5, 7), c(1, 4, 7), c(1, 5.5, 7))
TUs <- c(50, 75, 90, 100)
stopifnot(length(intensities) == length(sample_sizes) & length(sample_sizes) == length(rs))
length.each <- length(intensities)

if(run_sims) {
    ## data(pref_attach)
    ## g <- pref_attach
    data(list = which_network)
    g <- get(which_network)

    results <- list()

    for(i in 1:length.each) {
        print(intensities[i])
        set.seed(123)
        results[[i]] <- simulation(g, s = intensities[i], check_equil = TRUE)
    }

    for(i in 1:length.each) {
        print(sample_sizes[i])
        set.seed(123)
        results[[i + length.each]] <- simulation(g, nsamples = sample_sizes[i], check_equil = TRUE)
    }

    for(i in 1:length.each) {
        print(rs[[i]])
        set.seed(123)
        results[[i + 2*length.each]] <- simulation(g, r = rs[[i]], check_equil = TRUE)
    }

    for(i in 1:length.each) {
        print(TUs[i])
        set.seed(123)
        results[[i + 3*length.each]] <- simulation(g, TU = TUs[i], check_equil = TRUE)
    }

    parameter_variation_results <- results
    save(parameter_variation_results, file = filename)#"./data/parameter-variation-results.rda")
} else {
    load(filename)#"./data/parameter-variation-results.rda")
    results <- parameter_variation_results
}

if(save_results) {
    if(!dir.exists("./r-out")) dir.create("./r-out")
    sink(outfile, append = FALSE, type = "output", split = TRUE)
}

for(i in 1:length(results)) {
    if(i <= length.each) {
        print(paste0("Intensity: ", intensities[i]))
    } else if(i > length.each & i <= 2*length.each) {
        print(paste0("N Samples: ", sample_sizes[i - length.each]))
    } else if(i > 2*length.each & i <= 3*length.each) {
        print(paste0("r = ", rs[[i - 2*length.each]]))
    } else print(paste0("TU = ", TUs[i - 3*length.each]))
    corr_results <- Kendall_correlations(results[[i]]$results)
    print("Means")
    print(corr_results$means)
    ##print("SDs")
    ##print(corr_results$sds)
}

if(save_results) sink()

## Plot testing below here ##

## corr_results <- lapply(results, function(x) Kendall_correlations(x$results)$means)
## corrs <- do.call(cbind, corr_results)
## image(corrs)


## ## colorRamp testing
## plotdat <- round(t(corrs[nrow(corrs):1, ]), 3)

## yticklabels <- paste0(
##     c("Max. Eigenvalue", "... ", "... ",
##       "Max. SD", "... ", "... ",
##       "Avg. SD", "... ", "... ",
##       "Max. AC", "... ", "... ",
##       "Avg. AC", "... ", "... "),
##     c(": All", ": Lower", ": Sentinel")
## )
## xticklabels <- c("Intensity", "# Samples", "r", "TU")
## yticks <- 1:nrow(corrs)
## xticks <- 1:ncol(corrs)

## color <- colorRamp(c("white", "navy"))
## colorref <- round(seq(0, 1, length.out = 1000), 3)
## colors <- apply(color(colorref), 1, function(x) rgb(t(x), maxColorValue = 255))
## col <- matrix(colors[match(plotdat, colorref)], byrow = FALSE, nrow = nrow(plotdat))
## test <- matrix(colorref[match(plotdat, colorref)], byrow = FALSE, nrow = nrow(plotdat))
## col <- colors[match(plotdat, colorref)]

## ht <- 7
## wd <- 8
## if(save_plots) {
##     pdf("./img/parameter-variation.pdf", height = ht, width = wd)
## } else dev.new(height = ht, width = wd)
## par(mar = c(1, 10, 5, 1) + 0.1)
## image(y = yticks, x = xticks, z = plotdat,
##       col = col,#hcl.colors(10, "blues", rev = TRUE),
##       axes = FALSE, xlab = "", ylab = "")
## axis(side = 2, at = rev(yticks), labels = yticklabels, tick = FALSE, las = 1)
## axis(side = 3, at = seq(from = 2.5, by = 4, length.out = 4), labels = xticklabels,
##      tick = FALSE, adj = .5)
## mtext("Preferential Attachment", font = 2, line = 3, adj = -.5)
## if(save_plots) dev.off()



## library(lattice)

## ncuts <- 20
## levelplot(
##     plotdat,
##     cuts = ncuts,
##     col.regions = hcl.colors(ncuts+1, "blues", rev = TRUE),
##     xlab = "", ylab = "",
##     scales = list(y = list(at = rev(yticks), label = yticklabels),
##                   x = list(at = xticks, label = xticklabels))
## )
