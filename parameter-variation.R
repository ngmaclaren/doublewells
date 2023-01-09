## Code by Neil MacLaren 3/1/2022, edited 6/3/2022

library(igraph)
library(doublewells)

save_plots <- TRUE # FALSE
run_sims <- FALSE # TRUE
filename <- "./data/parameter-variation-results-r1.rda"

choices <- c("powerlaw", "dolphins")

                                        # Typical values commented below
intensities <- c(.01, .1, .5) # , .05
sample_sizes <- c(25, 50, 150) # 250, 
rs <- list(c(1, 3, 5), c(1, 2.5, 7), c(1, 5.5, 7)) # , c(1, 4, 7)
TUs <- c(50, 100, 125) # 75
                                        # All parameter vectors should be the same length to support
                                        # the loops below
stopifnot(length(intensities) == length(sample_sizes) & length(sample_sizes) == length(rs))
length.each <- length(intensities)
saveorder <- c("s", "M", "r", "TU")
n.param <- 4
                                        # Main loop.
                                        # For each of the two example networks, run the standard
                                        # simulations with the exception of the parameter being
                                        # adjusted.
if(run_sims) {
    data(list = choices)

    results <- vector(mode = "list", length = n.param * length.each)

    for(i in 1:length.each) {
        print(intensities[i])
        set.seed(123)
        results[[i]] <- vector(mode = "list", length = length(choices))
        for(j in 1:length(choices)) {
            g <- get(choices[j])
            results[[i]][[j]] <- simulation(g, s = intensities[i], check_alts = TRUE, TU = 50, assessment_samples = 1:10)
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
                                                        check_alts = TRUE, TU = 50, assessment_samples = 1:10)
        }
    }

    for(i in 1:length.each) {
        print(rs[[i]])
        set.seed(123)
        adjust <- 2*length.each
        results[[i + adjust]] <- vector(mode = "list", length = length(choices))
        for(j in 1:length(choices)) {
            g <- get(choices[j])
            results[[i + adjust]][[j]] <- simulation(g, r = rs[[i]], check_alts = TRUE, assessment_samples = 1:10)
        }
    }

    for(i in 1:length.each) {
        print(TUs[i])
        set.seed(123)
        adjust <- 3*length.each
        results[[i + adjust]] <- vector(mode = "list", length = length(choices))
        for(j in 1:length(choices)) {
            g <- get(choices[j])
            results[[i + adjust]][[j]] <- simulation(g, TU = TUs[i], check_alts = TRUE, assessment_samples = 1:10)
        }
    }

    names(results) <- as.character(t(outer(saveorder, 1:length.each, FUN = paste, sep = "_")))
    parameter_variation_results <- results
    save(parameter_variation_results, file = filename)
} else {
    load(filename)
    results <- parameter_variation_results
}

## Plotting
palette("R4")

## Loop through all the data frames, producing one plot, code below, for each one
nodesets <- c("all", "lower", "sentinel", "upper", "lowrank", "random", "largecorr", "largesd")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
pch <- colors <- 1:length(nodesets)
pt.lwd <- 3
pt.cex <- 3
adjust <- c(0.25, -0.25)
ht <- .75*(7*length.each)
wd <- 7*2
filenames <- paste0("./img/param-variation-", saveorder, ".pdf")
for(i in 1:n.param) {
    imgfile <- filenames[i]
    if(save_plots) {
        cairo_pdf(imgfile, height = ht, width = wd, family = "Bitstream Vera Sans")
    } else dev.new(height = ht, width = wd)
    par(mfrow = c(3, 1), mar = c(4, 8, 0, 2) + 0.5, xpd = TRUE) # 6.1
    yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
    ylim <- c(0.5, length(signals) + 0.5)
    xlim <- c(-1, 1)
    for(j in 1:length.each) {
        idx <- (i - 1)*length.each + j
        dl <- results[[idx]]
        plot(NULL, xlab = "", ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
             bty = "n", cex.axis = 2.5)
        title(xlab = "τ", cex.lab = 2.5, font.lab = 3)
        mtext(LETTERS[j], font = 2, line = -3, adj = -0.06, cex = 3)
        axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1,
             cex.lab = 1.75, cex.axis = 1.75)
        text(x = -1, y = c(5.25, 4.75), labels = c("Power-law network", "Dolphin network"),
             adj = c(0, .5), col = "black", cex = 1.5)
        segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
        segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5,
                 lty = 2, lwd = .5, col = "gray")
        for(k in 1:length(dl)) {
            df <- dl[[k]]
            rdf <- df[, -c(grep("lowinput", colnames(df)))]
            dat <- as.data.frame(Kendall_correlations(rdf)$means)
            dat$fullname <- rownames(dat)
            dat <- cbind(dat, do.call(rbind, strsplit(dat$fullname, "_")))
            dat <- dat[, -grep("fullname", colnames(dat))]
            colnames(dat) <- c("tau", "signal", "set")
            rownames(dat) <- NULL
            dat$signal <- factor(dat$signal, levels = rev(signals))
            dat$set <- factor(dat$set, levels = nodesets)
            dat$signal_ <- as.numeric(dat$signal) + adjust[k]
            dat$set_ <- as.numeric(dat$set)
            points(signal_ ~ tau, data = dat, col = dat$set_, pch = dat$set_,
                   lwd = pt.lwd, cex = pt.cex)
        }
        if(j == length.each) {
            legend("bottomleft", cex = 1.5, #inset = c(-0.01, 0),
                   legend = c("All", "Lower State", "High Input",
                              "Upper State", "Lower Half", "Random", "Large Corr.", "Large SD"),
                   pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
                   bty = "n")
        }
    }
    if(save_plots) dev.off()
}
