x <- seq(-.75, 8.75, length.out = 1000)
y <- x^4/4 - 4*x^3 + 39*x^2/2 - 28*x

palette("Okabe-Ito")

svg("./img/doublewell-potential-1-4-7.svg", height = 7, width = 10, bg = "transparent")
par(mar = rep(0, 4))
plot(x, y, type = "l", bty = "n", axes = FALSE, ylab = "", xlab = "", lwd = 3, col = 4)
dev.off()

y <- x^4/4 - 4*x^3 + 39*x^2/2 - 36*x

svg("./img/doublewell-potential-1-4-7-close.svg", height = 7, width = 10, bg = "transparent")
par(mar = rep(0, 4))
plot(x, y, type = "l", bty = "n", axes = FALSE, ylab = "", xlab = "", lwd = 3, col = 4)
dev.off()

y <- x^4/4 - 4*x^3 + 39*x^2/2 - 38*x
svg("./img/doublewell-potential-1-4-7-transition.svg", height = 7, width = 10, bg = "transparent")
par(mar = rep(0, 4))
plot(x, y, type = "l", bty = "n", axes = FALSE, ylab = "", xlab = "", lwd = 3, col = 4)
dev.off()
