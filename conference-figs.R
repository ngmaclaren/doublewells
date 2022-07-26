library(igraph)
library(doublewells)

save_plots <- TRUE # FALSE
##palette("Okabe-Ito")
palette("Tableau 10")

load("./data/examples-lower.rda")
##load("./data/examples-upper.rda")

imgfile <- "./complex-networks-2022/LatexAbstractTemplate/conference-fig.pdf"

ht <- 5
wd <- 20
if(save_plots) {
    pdf(imgfile, height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mfrow = c(1, 2))

## A
nodesets <- c("all", "lower", "sentinel")
signals <- c("avgac")
columns <- paste(signals, nodesets, sep = "_")
colors <- 1:length(nodesets)
ylims <- list(
    c(0, 110),
    c(0, 70)
)
par(mar = c(3.5, 4, 2, 4) + 0.4, xpd = TRUE)
df <- examples_lower[[1]]
df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
df$p_upperstate <- round( (1 - (df$n_lowerstate/max(df$n_lowerstate)))*100 )
plot(df$Ds, df$p_upperstate,#df$p_lowerstate,
     type = "l", lwd = 4, lty = 1,
     col = palette()[length(palette())],
     xlim = range(df$Ds), ylim = c(0, 100),
     xlab = expression(italic(D)), ylab = "% Nodes in Upper State", yaxt = "n",
     cex.lab = 1.75, cex.axis = 1.75)
axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100), cex.axis = 1.75)
legend("topright", inset = c(0, -.15),
       ncol = 4, bty = "n",
       lty = 1, lwd = 4,
       col = c(length(palette()), colors),
       cex = 1.25,
       legend = c("State", "All", "Lower", "High Input")
       )
mtext("A", adj = -0.08, line = 1, cex = 1.75, font = 2)
par(new = TRUE)
plotsamples <- 1:nrow(df)
matplot(df$Ds[plotsamples], df[plotsamples, columns],
        type = "l", lwd = 4, lty = 1,
        col = colors,
        xlim = range(df$Ds), ylim = c(0, 1), axes = FALSE,
        xlab = "", ylab = "")
axis(4, at = pretty(c(0, 1)), cex.axis = 1.75)
mtext("Average Autocorrelation",side = 4, cex = 1.75, font = 1, line = 3)

## B
nodesets <- c("all", "lower", "sentinel", "upper", "lowrank", "random")#, "largecorr", "largesd")
signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
pch <- colors <- 1:length(nodesets)
pt.lwd <- 4
pt.cex <- 3
par(mar = c(3.5, 7, 0, 1) + 0.4, xpd = TRUE)
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
ylim <- c(0.5, length(signals) + 0.5)
xlim <- c(-1, 1)
plot(NULL, xlab = expression(italic(tau)), ylab = "", yaxt = "n", ylim = ylim, xlim = xlim,
     bty = "n", cex.lab = 1.75, cex.axis = 1.75)
mtext("B", adj = -0.08, line = -1, cex = 1.75, font = 2)
axis(side = 2, at = 5:1, labels = yticklabels, tick = FALSE, las = 1, cex.lab = 1.75, cex.axis = 1.75)
legend("bottomleft", cex = 1.25, #inset = c(-.1, 0),#c(-.32, 0),
       legend = c("All", "Lower State", "High Input",
                  "Upper State", "Lower Half", "Random"),#, "Large Corr.", "Large SD"),
       pch = pch, col = colors, pt.cex = pt.cex*.75, pt.lwd = pt.lwd*.75,
       bty = "n")
segments(y0 = 1:4 + .5, x0 = -1, x1 = 1, lty = 1, lwd = .5, col = "gray")
segments(x0 = c(-.5, 0, .5), y0 = .5, y1 = length(signals) + .5, lty = 2, lwd = .5, col = "gray")
rdf <- df[, -c(grep("lowinput", colnames(df)))]
dat <- as.data.frame(Kendall_correlations(rdf)$means)
dat$fullname <- rownames(dat)
dat <- cbind(dat, do.call(rbind, strsplit(dat$fullname, "_")))
dat <- dat[, -grep("fullname", colnames(dat))]
colnames(dat) <- c("tau", "signal", "set")
rownames(dat) <- NULL
dat$signal <- factor(dat$signal, levels = rev(signals))
dat$set <- factor(dat$set, levels = nodesets)
dat$signal_ <- as.numeric(dat$signal)
dat$set_ <- as.numeric(dat$set)
points(signal_ ~ tau, data = dat, col = dat$set_, pch = dat$set_, lwd = pt.lwd, cex = pt.cex)

if(save_plots) dev.off()
