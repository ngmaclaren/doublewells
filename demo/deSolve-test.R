library(deSolve)

## Lorenz model

parameters <- c(a = -8/3, b = -10, c = 28)
state <- c(X = 1, Y = 1, Z = 1)

Lorenz <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
        dX <- a*X + Y*Z
        dY <- b*(Y - Z)
        dZ <- -X*Y + c*Y - Z

        list(c(dX, dY, dZ))
    })
}

times <- seq(0, 100, by = 0.01)

out <- deSolve::ode(y = state, times = times, func = Lorenz, parms = parameters)

head(out)

par(oma = c(0, 0, 3, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "X"], out[, "Z"], pch = ".")
mtext(outer = TRUE, side = 3, "Lorenz Model", cex = 1.5)

