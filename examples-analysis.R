## Code by Neil MacLaren 5/31/2022

### Update this file to produce the new early warning signals

library(igraph)
library(doublewells)

### This file produces the various examples figures

## options
save_plots <- TRUE # FALSE
palette("R4")

load("./data/examples-lower.rda")
load("./data/examples-upper.rda")

### Figure 1. use examples_lower
### x = D, y1 = prop lower state, y2 = avg ac
nodesets <- c("all", "lower", "sentinel")
signals <- c("avgac")
columns <- paste(signals, nodesets, sep = "_")
d_imgfile <- "./img/examples-multistage-transition.pdf"
wd <- 9
ht <- 8
colors <- 1:length(nodesets)
ylims <- list(
    c(0, 110),
    c(0, 70)
)
if(save_plots) {
    pdf(d_imgfile, width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mfrow = c(2, 1), mar = c(3.5, 4, 1, 4)+.4, xpd = TRUE)
for(i in 1:length(examples_lower)) {
    df <- examples_lower[[i]]
    df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
    df$p_upperstate <- round( (1 - (df$n_lowerstate/max(df$n_lowerstate)))*100 )
    plot(df$Ds, df$p_upperstate,#df$p_lowerstate,
         type = "l", lwd = 3, lty = 1,
         col = palette()[length(palette())],
         xlim = range(df$Ds), ylim = c(0, 100),
         xlab = expression(italic(D)), ylab = "% Nodes in Upper State", yaxt = "n",
         cex.lab = 1.75, cex.axis = 1.75)
    axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100), cex.axis = 1.75)
    if(i == 1) {
        legend("topright", inset = c(0, -.15),
               ncol = 4, bty = "n",
               lty = 1, lwd = 4,
               col = c(length(palette()), colors),
               cex = 1.25,
               legend = c("State", "All", "Lower", "High Input")
               )
    }
    mtext(LETTERS[i], adj = 0.01, line = -1.5, cex = 1.5, font = 2)
    par(new = TRUE)
    plotsamples <- 1:nrow(df)
    matplot(df$Ds[plotsamples], df[plotsamples, columns],
            type = "l", lwd = 2, lty = 1,
            col = colors,
            xlim = range(df$Ds), ylim = c(0, 1), axes = FALSE,
            xlab = "", ylab = "")
    axis(4, at = pretty(c(0, 1)), cex.axis = 1.75)
    mtext("Average Autocorrelation",side = 4, cex = 1.75, font = 1, line = 3)
}
if(save_plots) dev.off()

## Figure 2. use examples_lower
nodesets <- c("all", "lower", "sentinel", "upper", "lowrank", "random", "largecorr", "largesd")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
tau_imgfile <- "./img/examples-taus.pdf"
pch <- colors <- 1:length(nodesets)
pt.lwd <- 3
pt.cex <- 3
adjust <- c(0.25, -0.25)
ht <- 7
wd <- 14#7
if(save_plots) {
    pdf(tau_imgfile, height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(4, 7, 0, 6.1) + 0.3, xpd = TRUE)
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
ylim <- c(0.5, length(signals) + 0.5)
xlim <- c(-1, 1)
plot(NULL, xlab = expression(italic(tau)), ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
     bty = "n", cex.lab = 2.5, cex.axis = 2.5)#cex.lab = 1.75, cex.axis = 1.75)
axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1, cex.lab = 1.75, cex.axis = 1.75)
text(x = -1, y = c(5.25, 4.75), labels = c("BA", "Dolphins"),
     adj = c(0, .5), col = "black", font = 3, cex = 1.5)
legend("bottomright", cex = 1.25, inset = c(-.1, 0),#c(-.32, 0),
       legend = c("All", "Lower State", "High Input",
                  "Upper State", "Lower Half", "Random", "Large Corr.", "Large SD"),
       pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
       bty = "n")
segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .5, col = "gray")
for(i in 1:length(examples_lower)) {
    df <- examples_lower[[i]]
    rdf <- df[, -c(grep("lowinput", colnames(df)))]
    dat <- as.data.frame(Kendall_correlations(rdf)$means)
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
}
if(save_plots) dev.off()

## Figure 3. use examples_upper
nodesets <- c("all", "lower", "sentinel", "upper", "lowinput_sentinel", "random")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
tau_imgfile_u <- "./img/examples-upper-taus.pdf"
pch <- colors <- 1:length(nodesets)
pt.lwd <- 3
pt.cex <- 3
adjust <- c(0.25, -0.25)
ht <- 7
wd <- 14#7
if(save_plots) {
    pdf(tau_imgfile_u, height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(4, 7, 0, 6.1) + 0.3, xpd = TRUE)
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
ylim <- c(0.5, length(signals) + 0.5)
xlim <- c(-1, 1)
plot(NULL, xlab = expression(italic(tau)), ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
     bty = "n", cex.lab = 2.5, cex.axis = 2.5)
axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1, cex.lab = 1.75, cex.axis = 1.75)
text(x = .85, y = c(5.25, 4.75), labels = c("BA", "Dolphins"),
     adj = c(0, .5), col = "black", font = 3, cex = 1.5)
legend("bottomright", cex = 1.25, inset = c(-.12, 0),# c(-.32, 0),
       legend = c("All", "Lower State", "High Input", "Upper State", "Low Input", "Random"),
       pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
       bty = "n")
segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .5, col = "gray")
for(i in 1:length(examples_upper)) {
    df <- examples_upper[[i]]
    rdf <- df[, -c(grep("lowrank", colnames(df))
                 ##, grep("random", colnames(df))
                   )]
    dat <- as.data.frame(Kendall_correlations(rdf)$means)
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
u_imgfile <- "./img/examples-multistage-transition-upper.pdf"
wd <- 9
ht <- 8
colors <- c(1, 4, 5)#1:length(nodesets)
ylims <- list(
    c(0, 110),
    c(0, 70)
)
if(save_plots) {
    pdf(u_imgfile, width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mfrow = c(2, 1), mar = c(3.5, 4, 1, 4)+.4, xpd = TRUE)
for(i in 1:length(examples_upper)) {
    df <- examples_upper[[i]]
    df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
    df$p_upperstate <- round( (1 - (df$n_lowerstate/max(df$n_lowerstate)))*100 )
    plot(df$Ds, df$p_upperstate,#df$p_lowerstate,
         type = "l", lwd = 3, lty = 1,
         col = palette()[length(palette())],
         xlim = rev(range(df$Ds)), ylim = c(0, 100),
         xlab = expression(italic(D)), ylab = "% Nodes in Upper State", yaxt = "n",
         cex.lab = 1.75, cex.axis = 1.75)
    axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100), cex.axis = 1.75)
    if(i == 1) {
        legend("topright", inset = c(0, -.15),
               ncol = 4, bty = "n",
               lty = 1, lwd = 4,
               col = c(length(palette()), colors),
               cex = 1.25,
               legend = c("State", "All", "Upper", "Low Input")
               )
    }
    mtext(LETTERS[i], adj = 0.01, line = -1.5, cex = 1.5, font = 2)
    par(new = TRUE)
    plotsamples <- 1:nrow(df)
    matplot(df$Ds[plotsamples], df[plotsamples, columns],
            type = "l", lwd = 2, lty = 1,
            col = colors,
            xlim = rev(range(df$Ds)), ylim = c(0, 1), axes = FALSE,
            xlab = "", ylab = "")
    axis(4, at = pretty(c(0, 1)), cex.axis = 1.75)
    mtext("Average Autocorrelation",side = 4, cex = 1.75, font = 1, line = 3)
}
if(save_plots) dev.off()
