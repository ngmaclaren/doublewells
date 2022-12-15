## Code by Neil MacLaren 9/26/2022

library(igraph)
library(doublewells)

save_plots <- TRUE # FALSE
palette("R4")

load("./data/examples-lower.rda")
df <- examples_lower[[1]] # the SSF network
highavgac <- subset(df, avgac_sentinel > 0.65, select = c(n_lowerstate, Ds, avgac_sentinel))
topfiveDs <- highavgac[order(highavgac$avgac_sentinel, decreasing = TRUE), ]$Ds[1:5]

                                        #D <- # 1: X0.945X, 2: 0.325, 3: 0.550
Dmax <- 0.550
Dstep <- 0.005
nD <- 10
Ds <- rev(seq(from = Dmax, by = -Dstep, length.out = nD))

data(powerlaw, package = "doublewells")
g <- powerlaw

r <- c(1, 4, 7)
s <- 0.05
sample_spacing <- 0.1
lag <- 1
dt <- 0.01
##TU <- 100 + 100

nnodes <- vcount(g)
A <- as_adj(g, type = "both", sparse = FALSE)

##T <- TU/dt
s <- s*sqrt(dt)
sample_spacing <- sample_spacing/dt

M <- 200 # the number of samples in each grab
betweenM <- 30/dt # the time steps between taking M samples
ngrabs <- 100 # the number of times to grab M samples

grab_sample <- function(t, M) seq(from = t, by = -sample_spacing, length.out = M)
samples <- lapply(seq(from = 100/dt, to = (100 + ngrabs*50)/dt, length.out = 100),
                  function(t) grab_sample(t, M))

T <- max(unlist(samples))

initialx <- rep(min(r), nnodes)
results <- list()

for(i in 1:length(Ds)) {
    D <- Ds[i]
    print(paste(i, "; D =", Ds[i]))
    x <- initialx
    X <- matrix(0, nrow = T, ncol = nnodes)
    for(t in 1:T) {
        X[t, ] <- x
        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s)) # u?
    }

    calc_sd <- function(X, samp) {# indices
        apply(X[samp, ], 2, sd)
    }
    result <- lapply(samples, function(samp) calc_sd(X, samp))

    result <- do.call(rbind, result) # Matrix with nsamples x nnodes
    means <- colMeans(result)
    sds <- apply(result, 2, sd)
    results[[i]] <- list(means = means, sds = sds)
}

### easiest thing may be to use points/segments to plot in groups of four according to D
points_and_segments <- function(D, result, nodes) {
    n <- length(nodes)
    colors <- 2:(n+1)
    points(x = rep(D, n), y = result$means[nodes], pch = 19, col = colors, cex = .75)
    segments(x0 = rep(D, n),
             y0 = result$means[nodes] - result$sds[nodes],
             y1 = result$means[nodes] + result$sds[nodes],
             lwd = 3, lty = 1, col = colors)
}

lastresult <- results[[length(results)]]
sortnodes <- order(lastresult$means)
selectnodes <- rev(c(1, 10, length(sortnodes) - 1, length(sortnodes)))
nodes <- sortnodes[selectnodes]

wd <- 7; ht <- 7;
ymax <- 0.05
if(save_plots) {
    pdf(paste0("./img/longsim-errorbars-", Dmax, ".pdf"), width = wd, height = ht)
} else {
    dev.new(height = ht, width = wd)
}
plot(NULL, xlim = range(Ds), ylim = c(0, ymax), xlab = "D", ylab = "m_i")
for(i in 1:length(Ds)) points_and_segments(Ds[i], results[[i]], nodes)
if(save_plots) dev.off()


## wd <- 7; ht <- 7;
## xmax <- max(means); ymax <- max(sds);
## if(save_plots) {
##     pdf(paste0("./img/longsim-", D, ".pdf"), width = wd, height = ht)
## } else {
##     dev.new(height = ht, width = wd)
## }
## plot(
##     sds ~ means, xlab = "m_i", ylab = "s_i", main = paste("D =", D),
##     pch = 1, col = "gray20", lwd = 1.5,
##     xlim = c(0, xmax*1.1), ylim = c(0, ymax*1.1)#log = "xy"
## )
## if(save_plots) dev.off()
