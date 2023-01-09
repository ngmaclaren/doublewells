## Code by Neil MacLaren 12/15/2022

library(igraph)
library(doublewells)

save_plots <- TRUE # FALSE
run_sims <- FALSE # TRUE
choices <- c("powerlaw", "dolphins")
palette("R4")

data(list = choices)

outfile_lower <- "./data/examples-lower-u.rda"
outfile_upper <- "./data/examples-upper-u.rda"

if(isTRUE(run_sims)) {
    set.seed(123)

    examples_lower <- list()
    examples_upper <- list()

    for(i in 1:length(choices)) {
        print(choices[i])

        g <- get(choices[i])

        print("lower")
        examples_lower[[i]] <- simulation(
            g,
            check_alts = TRUE,
            bifurcation_parameter = "u",
            D.init = 0.1,
            stepsize = 0.05,
            return_histories = TRUE
        )

        print("upper")
        examples_upper[[i]] <- simulation(
            g,
            from_upper = TRUE,
            check_alts = TRUE,
            bifurcation_parameter = "u",
            D.init = 0.1,
            stepsize = -0.05,
            return_histories = TRUE
        )
    }

    save(examples_lower, file = outfile_lower)
    save(examples_upper, file = outfile_upper)
} else {
    load(outfile_lower)
    load(outfile_upper)
}

nodesets <- c("all", "lower", "sentinel")
signals <- c("avgac")
columns <- paste(signals, nodesets, sep = "_")
d_imgfile <- "./img/examples-multistage-transition-u.pdf"
wd <- 9
ht <- 8
colors <- 1:length(nodesets)
ylims <- list(
    c(0, 110),
    c(0, 70)
)
if(save_plots) {
    cairo_pdf(d_imgfile, width = wd, height = ht, family = "Bitsream Vera Sans")
} else dev.new(width = wd, height = ht)
par(mfrow = c(2, 1), mar = c(3.5, 4, 1, 4)+.4, xpd = TRUE)
for(i in 1:length(examples_lower)) {
    df <- examples_lower[[i]]$df
    df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
    df$p_upperstate <- round( (1 - (df$n_lowerstate/max(df$n_lowerstate)))*100 )
    plot(df$us, df$p_upperstate,
         type = "l", lwd = 3, lty = 1,
         col = palette()[length(palette())],
         xlim = range(df$us), ylim = c(0, 100),
         xlab = expression(italic(u)), ylab = "% Nodes in Upper State", yaxt = "n",
         cex.lab = 1.75, cex.axis = 1.75)
    axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100), cex.axis = 1.75)
    cvs <- rev(find_critical_values(examples_lower[[i]], bifurcation_parameter = "u"))
    arrows(
        x0 = cvs - .1, x1 = cvs,
        y0 = df[df$us %in% cvs, "p_upperstate"] + 10, y1 = df[df$us %in% cvs, "p_upperstate"] + 1,
        length = .1, angle = 20
    )
    if(i == 1) {
        legend("topright", inset = c(0, -.15),
               ncol = 4, bty = "n",
               lty = 1, lwd = 4,
               col = c(length(palette()), colors),
               cex = 1.25,
               legend = c("State", "All", "Lower", "High Input")
               )
    }
    mtext(LETTERS[i], adj = 0.01, line = -2.5, cex = 3, font = 2)
    par(new = TRUE)
    plotsamples <- 1:nrow(df)
    matplot(df$us[plotsamples], df[plotsamples, columns],
            type = "l", lwd = 2, lty = 1,
            col = colors,
            xlim = range(df$us), ylim = c(0, 1), axes = FALSE,
            xlab = "", ylab = "")
    axis(4, at = pretty(c(0, 1)), cex.axis = 1.75)
    mtext("Average Autocorrelation",side = 4, cex = 1.75, font = 1, line = 3)
}
if(save_plots) dev.off()

nodesets <- c("all", "lower", "sentinel", "upper", "lowrank", "random", "largecorr", "largesd")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
tau_imgfile <- "./img/examples-taus-u.pdf"
pch <- colors <- 1:length(nodesets)
pt.lwd <- 3
pt.cex <- 3
adjust <- c(0.25, -0.25)
ht <- 7
wd <- 14
if(save_plots) {
    cairo_pdf(tau_imgfile, height = ht, width = wd, family = "Bitsream Vera Sans")
} else dev.new(height = ht, width = wd)
par(mar = c(4, 7, 0, 6.1) + 0.5, xpd = TRUE)
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
ylim <- c(0.5, length(signals) + 0.5)
xlim <- c(-1, 1)
plot(NULL, xlab = "", ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
     bty = "n", cex.axis = 2.5)
title(xlab = "τ", cex.lab = 2.5, font.lab = 3)
axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1, cex.lab = 1.75, cex.axis = 1.75)
text(x = -1, y = c(5.25, 4.75), labels = c("Power-law network", "Dolphin network"),
     adj = c(0, .5), col = "black", cex = 1.5)
legend("bottomright", cex = 1.25, inset = c(-.1, 0),
       legend = c("All", "Lower State", "High Input",
                  "Upper State", "Lower Half", "Random", "Large Corr.", "Large SD"),
       pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
       bty = "n")
segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .5, col = "gray")
for(i in 1:length(examples_lower)) {
    df <- examples_lower[[i]]$df
    rdf <- df[, !grepl("lowinput", colnames(df))]
    dat <- as.data.frame(Kendall_correlations(rdf, b_param = "us")$means)
    dat$fullname <- rownames(dat)
    dat <- cbind(dat, do.call(rbind, strsplit(dat$fullname, "_")))
    dat <- dat[, -grep("fullname", colnames(dat))]
    colnames(dat) <- c("tau", "signal", "set")
    rownames(dat) <- NULL
    dat$signal <- factor(dat$signal, levels = rev(signals))
    dat$set <- factor(dat$set, levels = nodesets)
    dat$signal_ <- as.numeric(dat$signal) + adjust[i]
    dat$set_ <- as.numeric(dat$set)
    points(signal_ ~ tau, data = dat, col = dat$set_, pch = dat$set_, lwd = pt.lwd, cex = pt.cex)
    print(summary(dat$tau[dat$set %in% c("all", "lower", "sentinel")]))
}
if(save_plots) dev.off()

nodesets <- c("all", "lower", "sentinel", "upper", "lowinput_sentinel", "random")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
tau_imgfile_u <- "./img/examples-upper-taus-u.pdf"
pch <- colors <- 1:length(nodesets)
pt.lwd <- 3
pt.cex <- 3
adjust <- c(0.25, -0.25)
ht <- 7
wd <- 14
if(save_plots) {
    cairo_pdf(tau_imgfile_u, height = ht, width = wd, family = "Bitsream Vera Sans")
} else dev.new(height = ht, width = wd)
par(mar = c(4, 7, 0, 6.1) + 0.5, xpd = TRUE)
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
ylim <- c(0.5, length(signals) + 0.5)
xlim <- c(-1, 1)
plot(NULL, xlab = "", ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
     bty = "n", cex.axis = 2.5)
title(xlab = "τ", cex.lab = 2.5, font.lab = 3)
axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1, cex.lab = 1.75, cex.axis = 1.75)
text(x = .85, y = c(5.25, 4.75), labels = c("Power-law network", "Dolphin network"),
     adj = c(0, .5), col = "black", cex = 1.5)
legend("bottomright", cex = 1.25, inset = c(-.12, 0),
       legend = c("All", "Lower State", "High Input", "Upper State", "Low Input", "Random"),
       pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
       bty = "n")
segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .5, col = "gray")
for(i in 1:length(examples_upper)) {
    df <- examples_upper[[i]]$df
    rdf <- df[, !grepl("lowrank", colnames(df))]
    dat <- as.data.frame(Kendall_correlations(rdf, b_param = "us")$means)
    dat$fullname <- rownames(dat)
    dat$fullname <- gsub("lowinput_sentinel", "lowinputSentinel", dat$fullname)
    dat <- cbind(dat, do.call(rbind, strsplit(dat$fullname, "_")))
    dat <- dat[, -grep("fullname", colnames(dat))]
    colnames(dat) <- c("tau", "signal", "set")
    dat$set <- gsub("lowinputSentinel", "lowinput_sentinel", dat$set)
    rownames(dat) <- NULL
    dat$signal <- factor(dat$signal, levels = rev(signals))
    dat$set <- factor(dat$set, levels = nodesets)
    dat$signal_ <- as.numeric(dat$signal) + adjust[i]
    dat$set_ <- as.numeric(dat$set)
    points(signal_ ~ tau, data = dat, col = dat$set_, pch = dat$set_, lwd = pt.lwd, cex = pt.cex)
}
if(save_plots) dev.off()

### Figure 4(SI). use examples_upper
nodesets <- c("all", "upper", "lowinput_sentinel")
signals <- c("avgac")
columns <- paste(signals, nodesets, sep = "_")
u_imgfile <- "./img/examples-multistage-transition-upper-u.pdf"
wd <- 9
ht <- 8
colors <- c(1, 4, 5)
ylims <- list(
    c(0, 110),
    c(0, 70)
)
if(save_plots) {
    cairo_pdf(u_imgfile, width = wd, height = ht, family = "Bitsream Vera Sans")
} else dev.new(width = wd, height = ht)
par(mfrow = c(2, 1), mar = c(3.5, 4, 1, 4)+.4, xpd = TRUE)
for(i in 1:length(examples_upper)) {
    df <- examples_upper[[i]]$df
    df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
    df$p_upperstate <- round( (1 - (df$n_lowerstate/max(df$n_lowerstate)))*100 )
    plot(df$us, df$p_upperstate,
         type = "l", lwd = 3, lty = 1,
         col = palette()[length(palette())],
         xlim = rev(range(df$us)), ylim = c(0, 100),
         xlab = expression(italic(u)), ylab = "% Nodes in Upper State", yaxt = "n",
         cex.lab = 1.75, cex.axis = 1.75)
    axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100), cex.axis = 1.75)
    cvs <- find_critical_values(examples_upper[[i]], bifurcation_parameter = "u")
    cvs <- sapply(as.integer(names(cvs)), function(x) min(df$us[which(df$n_lowerstate == x)]))
    arrows(
        x0 = cvs - .1, x1 = cvs,
        y0 = df[df$us %in% cvs, "p_upperstate"] + 10, y1 = df[df$us %in% cvs, "p_upperstate"] + 1,
        length = .1, angle = 20
    )
    if(i == 1) {
        legend("topright", inset = c(0, -.15),
               ncol = 4, bty = "n",
               lty = 1, lwd = 4,
               col = c(length(palette()), colors),
               cex = 1.25,
               legend = c("State", "All", "Upper", "Low Input")
               )
    }
    mtext(LETTERS[i], adj = 0.01, line = -2.75, cex = 3, font = 2)
    par(new = TRUE)
    plotsamples <- 1:nrow(df)
    matplot(df$us[plotsamples], df[plotsamples, columns],
            type = "l", lwd = 2, lty = 1,
            col = colors,
            xlim = rev(range(df$us)), ylim = c(0, 1), axes = FALSE,
            xlab = "", ylab = "")
    axis(4, at = pretty(c(0, 1)), cex.axis = 1.75)
    mtext("Average Autocorrelation",side = 4, cex = 1.75, font = 1, line = 3)
}
if(save_plots) dev.off()
