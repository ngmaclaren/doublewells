xeq1 <- function(r) sqrt(r)
xeq2 <- function(r) -sqrt(r)

domain <- seq(0, 10, length.out = 1000)

plot(domain, xeq1(domain), col = "blue", type = "l", lwd = 3,
     xlim = c(-10, 10), ylim = c(-5, 5), xlab = "r", ylab = "x_eq")
lines(domain, xeq2(domain), col = "red", lwd = 3, lty = 2)
points(0, 0, col = "darkgreen", pch = 19)

