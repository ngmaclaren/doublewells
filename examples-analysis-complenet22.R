## Code by Neil MacLaren 5/31/2022

## Thicken lines
## new color palette
## save as SVG with no background
## split into components (i.e., not A/B figures)

library(igraph)
library(doublewells)

### This file produces the various examples figures

## options
save_plots <- FALSE # TRUE
                                        #palette("R4")
palette("Okabe-Ito")
newpal <- palette.colors()
newpal["yellow"] <- "#DBD03D"
palette(newpal)

load("./data/examples-lower.rda")
load("./data/examples-upper.rda")

### Figure 1. use examples_lower
### x = D, y1 = prop lower state, y2 = avg ac
nodesets <- c("all", "lower", "sentinel")
signals <- c("avgac")
columns <- paste(signals, nodesets, sep = "_")
d_imgfile <- "./img/examples-multistage-transition-complenet22.svg"
wd <- 12
ht <- 6
colors <- 2:(length(nodesets)+1)
lwd <- 3
fontscale <- 1.6
ylims <- c(0, 70)
if(save_plots) {
    ##pdf(d_imgfile, width = wd, height = ht)
    svg(d_imgfile, width = wd, height = ht, bg = "transparent")
} else dev.new(width = wd, height = ht)
par(mar = c(3.5, 4, 1, 4)+.4, xpd = TRUE)
df <- examples_lower[[2]]
df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
df$p_upperstate <- round( (1 - (df$n_lowerstate/max(df$n_lowerstate)))*100 )
plot(df$Ds, df$p_upperstate,#df$p_lowerstate,
     type = "l", lwd = lwd, lty = 1,
     col = 1,#palette()[length(palette())],
     xlim = range(df$Ds), ylim = c(0, 100),
     xlab = expression(italic(D)), ylab = "% Nodes in Upper State", yaxt = "n",
     cex.lab = fontscale, cex.axis = fontscale)
axis(side = 2, at = c(0, 25, 50, 75, 100), labels = c(0, 25, 50, 75, 100), cex.axis = fontscale)
legend("topleft", #"topright", inset = c(0, -.09),
       #ncol = 4,
       bty = "n",
       lty = 1, lwd = 4,
       col = c(1, colors),
       cex = 1.25,
       legend = c("State", "All", "Lower State", "High Input")
       )
par(new = TRUE)
plotsamples <- 1:nrow(df)
matplot(df$Ds[plotsamples], df[plotsamples, columns],
        type = "l", lwd = lwd, lty = 1,
        col = colors,
        xlim = range(df$Ds), ylim = c(0, 1), axes = FALSE,
        xlab = "", ylab = "")
axis(4, at = pretty(c(0, 1)), cex.axis = fontscale)
mtext("Average Autocorrelation",side = 4, cex = fontscale, font = 1, line = 3)
if(save_plots) dev.off()

## Figure 2. use examples_lower
##nodesets <- c("all", "lower", "sentinel", "upper", "lowrank", "random", "largecorr", "largesd")
nodesets <- c("all", "lower", "sentinel", "upper", "random")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
tau_imgfile <- "./img/examples-taus-complenet22.svg"
pch <- 1:length(nodesets)
colors <- 2:(length(nodesets) + 1)
pt.lwd <- 4
pt.cex <- 3
## adjust <- c(0.25, -0.25)
ht <- 6
wd <- 12#7
if(save_plots) {
    svg(tau_imgfile, height = ht, width = wd, bg = "transparent")#"#c2cdd9")#
} else dev.new(height = ht, width = wd)
par(mar = c(4, 7, 0, 1) + 0.3)#6.1, xpd = TRUE)
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
ylim <- c(0.5, length(signals) + 0.5)
xlim <- c(-1, 1)
plot(NULL, xlab = "Rank Correlation", #xlab = expression(italic(tau)),
     ylab = "", yaxt = "n", ylim = ylim, xlim = xlim, #bty = "n",
     cex.lab = fontscale, cex.axis = fontscale)#cex.lab = 1.75, cex.axis = 1.75)
axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1,
     cex.lab = fontscale, cex.axis = fontscale)
## text(x = -1, y = c(5.25, 4.75), labels = c("Power-Law network", "Dolphin network"),
##      adj = c(0, .5), col = "black", font = 3, cex = 1.5)
legend("bottomleft", cex = 1.25, #inset = c(-.1, 0),#c(-.32, 0),
       legend = c("All", "Lower State", "High Input", "Upper State", "Random"),
       pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
       bty = "n")
segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .9, col = "gray")
segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .9, col = "gray")
df <- examples_lower[[2]]
rdf <- df[, -c(grep("lowrank", colnames(df)), grep("large", colnames(df)),
               grep("lowinput", colnames(df)))]
dat <- as.data.frame(Kendall_correlations(rdf)$means)
dat$fullname <- rownames(dat)
dat <- cbind(dat, do.call(rbind, strsplit(dat$fullname, "_")))
dat <- dat[, -grep("fullname", colnames(dat))]
colnames(dat) <- c("tau", "signal", "set")
rownames(dat) <- NULL
dat$signal <- factor(dat$signal, levels = rev(signals))
dat$set <- factor(dat$set, levels = nodesets)
dat$signal_ <- as.numeric(dat$signal) + seq(-.5, .5, length.out = length(colors)+2)[colors]#adjust[2]
dat$set_ <- as.numeric(dat$set)
points(signal_ ~ tau, data = dat, col = dat$set_ + 1, pch = dat$set_, lwd = pt.lwd, cex = pt.cex)
print(summary(dat$tau[dat$set %in% c("all", "lower", "sentinel")]))
if(save_plots) dev.off()

## Figure 3. use examples_upper
nodesets <- c("all", "lower", "sentinel", "upper", "random", "lowinput_sentinel")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
tau_imgfile_u <- "./img/examples-upper-taus-complenet22.svg"
pch <- colors <- 1:length(nodesets)
pt.lwd <- 4
pt.cex <- 3
ht <- 6
wd <- 12
if(save_plots) {
    svg(tau_imgfile_u, height = ht, width = wd, bg = "transparent")
} else dev.new(height = ht, width = wd)
par(mar = c(4, 7, 0, 1) + 0.3)#, xpd = TRUE)
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
ylim <- c(0.5, length(signals) + 0.5)
xlim <- c(-1, 1)
plot(NULL, xlab = "Rank Correlation", ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
     bty = "o", cex.lab = fontscale, cex.axis = fontscale)
axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1,
     cex.lab = fontscale, cex.axis = fontscale)
legend("bottomright", cex = 1.25, #inset = c(-.12, 0),# c(-.32, 0),
       legend = c("All", "Lower State", "High Input", "Upper State", "Random", "Low Input"),
       pch = pch, col = colors+1, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
       bty = "n")
segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .5, col = "gray")
df <- examples_upper[[2]]
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
dat$signal_ <- as.numeric(dat$signal)
dat$set_ <- as.numeric(dat$set)
points(signal_ ~ tau, data = dat, col = dat$set_+1, pch = dat$set_, lwd = pt.lwd, cex = pt.cex)
if(save_plots) dev.off()

### Figure 4(SI). use examples_upper
nodesets <- c("all", "upper", "lowinput_sentinel")
signals <- c("avgac")
columns <- paste(signals, nodesets, sep = "_")
u_imgfile <- "./img/examples-multistage-transition-upper-complenet22.svg"
wd <- 12
ht <- 6
lwd <- 3
colors <- c(2, 5, 7)#1:length(nodesets)
if(save_plots) {
    svg(u_imgfile, width = wd, height = ht, bg = "transparent")
} else dev.new(width = wd, height = ht)
par(mar = c(3.5, 4, 1, 4)+.4, xpd = TRUE)
df <- examples_upper[[2]]
df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
df$p_upperstate <- round( (1 - (df$n_lowerstate/max(df$n_lowerstate)))*100 )
plot(df$Ds, df$p_upperstate,#df$p_lowerstate,
     type = "l", lwd = lwd, lty = 1,
     col = 1,#palette()[length(palette())],
     xlim = rev(range(df$Ds)), ylim = c(0, 100),
     xlab = expression(italic(D)), ylab = "% Nodes in Upper State", yaxt = "n",
     cex.lab = 1.75, cex.axis = 1.75)
axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100), cex.axis = 1.75)
legend(##"topleft", #inset = c(0, -.15),
    ##ncol = 4,
    x = 1, y = 95,
    bty = "n",
    lty = 1, lwd = 4,
    col = c(1, colors),
    cex = 1.25,
    legend = c("State", "All", "Upper State", "Low Input")
)
par(new = TRUE)
plotsamples <- 1:nrow(df)
matplot(df$Ds[plotsamples], df[plotsamples, columns],
        type = "l", lwd = lwd, lty = 1,
        col = colors,
        xlim = rev(range(df$Ds)), ylim = c(0, 1), axes = FALSE,
        xlab = "", ylab = "")
axis(4, at = pretty(c(0, 1)), cex.axis = 1.75)
mtext("Average Autocorrelation",side = 4, cex = 1.75, font = 1, line = 3)
if(save_plots) dev.off()
