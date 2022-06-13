## Code by Neil MacLaren 3/1/2022, edited 6/3/2022

library(igraph)
library(doublewells)

save_plots <- TRUE # FALSE
save_results <- FALSE # TRUE
run_sims <- FALSE # TRUE
filename <- "./data/parameter-variation-results.rda"

## which_network <- "pref_attach" # "dolphins" # 
## filename <- paste0("./data/parameter-variation-results-", which_network, ".rda")
## outfile <- paste0("./r-out/parameter-variation-", which_network, ".txt")
choices <- c("pref_attach", "dolphins")

intensities <- c(.01, .1, .5) # , .05
sample_sizes <- rev(c(150, 50, 25)) # 250, 
rs <- list(c(1, 3, 5), c(1, 2.5, 7), c(1, 5.5, 7)) # , c(1, 4, 7)
TUs <- c(50, 90, 100) # 75
stopifnot(length(intensities) == length(sample_sizes) & length(sample_sizes) == length(rs))
length.each <- length(intensities)
n.param <- 4

if(run_sims) {
    ## data(pref_attach)
    ## g <- pref_attach
    data(list = choices)#which_network)
    ##g <- get(which_network)

    results <- vector(mode = "list", length = n.param * length.each)

    for(i in 1:length.each) {
        print(intensities[i])
        set.seed(123)
        results[[i]] <- vector(mode = "list", length = length(choices))
        for(j in 1:length(choices)) {
            g <- get(choices[j])
            results[[i]][[j]] <- simulation(g, s = intensities[i], check_alts = TRUE, TU = 50)
                                        #, check_equil = TRUE)
        }
    }

    for(i in 1:length.each) {
        print(sample_sizes[i])
        set.seed(123)
        adjust <- length.each
        results[[i + adjust]] <- vector(mode = "list", length = length(choices))
        for(j in 1:length(choices)) {
            g <- get(choices[j])
            results[[i + adjust]][[j]] <- simulation(g, nsamples = sample_sizes[i],
                                                        check_alts = TRUE, TU = 50)
                                        #, check_equil = TRUE)
        }
    }

    for(i in 1:length.each) {
        print(rs[[i]])
        set.seed(123)
        adjust <- 2*length.each
        results[[i + adjust]] <- vector(mode = "list", length = length(choices))
        for(j in 1:length(choices)) {
            g <- get(choices[j])
            results[[i + adjust]][[j]] <- simulation(g, r = rs[[i]], check_alts = TRUE)
                                        #, check_equil = TRUE)
        }
    }

    for(i in 1:length.each) {
        print(TUs[i])
        set.seed(123)
        adjust <- 3*length.each
        results[[i + adjust]] <- vector(mode = "list", length = length(choices))
        for(j in 1:length(choices)) {
            g <- get(choices[j])
            results[[i + adjust]][[j]] <- simulation(g, TU = TUs[i], check_alts = TRUE)
                                        #, check_equil = TRUE)
        }
    }

    parameter_variation_results <- results
    save(parameter_variation_results, file = filename)
} else {
    load(filename)
    results <- parameter_variation_results
}

## Plotting
pal <- c(
    ## https://personal.sron.nl/~pault/#sec:qualitative, the muted qualitative colour scheme
    "#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE",
    "#882255", "#44AA99", "#999933", "#AA4499", "#DDDDDD"
    
)
palette(pal)

## Loop through all the data frames, producing one plot, code below, for each one
nodesets <- c("all", "lower", "sentinel", "upper", "lowrank", "antisentinel", "altsentinel")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
pch <- colors <- 1:length(nodesets)
pt.lwd <- 3
pt.cex <- 3
adjust <- c(0.25, -0.25)
ht <- 7
wd <- 14#7
saveorder <- c("s", "M", "r", "TU")
filenames <- as.character(t(outer(saveorder, 1:length.each, FUN = paste, sep = "-")))
for(i in 1:length(results)) {
    dl <- results[[i]]
    imgfile <- paste0("./img/param-variation-", filenames[i], ".pdf")
    if(save_plots) {
        pdf(imgfile, height = ht, width = wd)
    } else dev.new(height = ht, width = wd)
    par(mar = c(4, 7, 0, 6.1) + 0.3, xpd = TRUE)
    yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
    ylim <- c(0.5, length(signals) + 0.5)
    xlim <- c(-1, 1)
    plot(NULL, xlab = expression(italic(tau)), ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
         bty = "n", cex.lab = 1.75, cex.axis = 1.75)
    axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1, cex.lab = 1.75, cex.axis = 1.75)
    text(x = -1, y = c(5.25, 4.75), labels = c("BA", "Dolphins"), adj = c(0, .5), col = "black", font = 3)
    legend("bottomright", inset = c(-.1, 0), cex = 1.25,
           legend = c("All", "Lower State", "Sentinel", "Upper State", "Low Rank", "Random", "Alternate"),
           pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
           bty = "n")
    segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
    segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .5, col = "gray")
    for(j in 1:length(dl)) {
        df <- dl[[j]]
        rdf <- df[, -c(grep("revdir", colnames(df)))]
        dat <- as.data.frame(Kendall_correlations(rdf)$means)
        dat$fullname <- rownames(dat)
        dat <- cbind(dat, do.call(rbind, strsplit(dat$fullname, "_")))
        dat <- dat[, -grep("fullname", colnames(dat))]
        colnames(dat) <- c("tau", "signal", "set")
        rownames(dat) <- NULL
        dat$signal <- factor(dat$signal, levels = rev(signals))
        dat$set <- factor(dat$set, levels = nodesets)
        dat$signal_ <- as.numeric(dat$signal) + adjust[j]
        dat$set_ <- as.numeric(dat$set)
        points(signal_ ~ tau, data = dat, col = dat$set_, pch = dat$set_, lwd = pt.lwd, cex = pt.cex)
    }
    if(save_plots) dev.off()
}



## Uncomment and fix to print out results from the new format.
## if(save_results) {
##     if(!dir.exists("./r-out")) dir.create("./r-out")
##     sink(outfile, append = FALSE, type = "output", split = TRUE)
## }
## for(i in 1:length(results)) {
##     if(i <= length.each) {
##         print(paste0("Intensity: ", intensities[i]))
##     } else if(i > length.each & i <= 2*length.each) {
##         print(paste0("N Samples: ", sample_sizes[i - length.each]))
##     } else if(i > 2*length.each & i <= 3*length.each) {
##         print(paste0("r = ", rs[[i - 2*length.each]]))
##     } else print(paste0("TU = ", TUs[i - 3*length.each]))
##     corr_results <- Kendall_correlations(results[[i]]$results)
##     print("Means")
##     print(corr_results$means)
##     ##print("SDs")
##     ##print(corr_results$sds)
## }
## if(save_results) sink()
