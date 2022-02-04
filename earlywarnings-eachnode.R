## Code by Neil MacLaren 2/3/2022

library(igraph)
library(doublewells)

whichstate <- function(x, cutoff = 7) ifelse(x > cutoff, 1, 0)
movingup <- function(x, cutoff = 4) ifelse(x > cutoff, 1, 0)

## set.seed(1)

which_net <- "pa" # "me"
which_param <- "D" # "u"
save_plot <- TRUE

## Graph
generate_new <- FALSE 
if(generate_new) {
    nnodes <- 100
    if(which_net == "me") {
        cprob <- 0.06
        g <- sample_gnp(nnodes, cprob)
        write_graph(g, "./data/testgraph.gml", format = "gml")
    } else if(which_net == "pa") {
        outdist <- c(0, .6, .5, .4, .3, .2, .1, .05)
        g <- sample_pa(nnodes, power = 1.5, out.dist = outdist, directed = FALSE)
        write_graph(g, "./data/testgraph-BA.gml", format = "gml")
    }
} else {
    if(which_net == "me") {
        g <- read_graph("./data/testgraph.gml", format = "gml")
    } else if(which_net == "pa") {
        g <- read_graph("./data/testgraph-BA.gml", format = "gml")
    }
    nnodes <- vcount(g)
}
A <- as_adj(g, type = "both", sparse = FALSE)

## Node Systems
r <- c(1, 4, 7) # double well parameters
s <- 0.005 # sd of noise process
### the below may need to be fixed to support varying u
if(which_net == "me") D <- 0.5 else if(which_net == "pa") D <- 0.21 # connection strength
p <- if(s > 0) 3*s else 0.015 # perturbation strength
u <- rep(0, nnodes) # stress vector

dt <- 0.01
T <- 2000

initialx <- rep(1, nnodes) + noise(nnodes, s)
##X <- matrix(0, nrow = T, ncol = nnodes)
x <- initialx
    
## Stress
stressnode <- select_stressnode(g, "highest")
if(which_param == "u") u[stressnode] <- 1

## Simulation
atrisk <- 1:nnodes
stepT <- 500
wl <- .5*stepT
ews <- list()
Ds <- numeric()
Is <- numeric()
i <- 1
cutoff <- .25*nnodes
n_atrisk <- numeric()
while(length(atrisk) > cutoff) {
    X <- matrix(0, nrow = stepT, ncol = nnodes)

    for(t in 1:stepT) {
        X[t, ] <- x
        x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
    }

    ## Analysis
    atrisk <- 1:nnodes
    upwardbound <- movingup(x)
    atrisk <- atrisk[which(upwardbound == 0)]
    n_atrisk[i] <- length(atrisk)
    if(n_atrisk[i] <= cutoff) break

    ## ## Debugging
    ## print(D)
    ## print(x)
    ## print(atrisk)

    ## sd method
    ew <- apply(X[(wl+1):stepT, ], 2, sd)
    ## ac method
    ##window <- (wl+1):stepT
    ##ew <- apply(X, 2, function(x) cor(x[window], x[window-10]))
    ew[-atrisk] <- NA
    I <- ape::Moran.I(X[stepT, ], A)$observed

    ews[[i]] <- ew
    Is[i] <- I
    Ds[i] <- D

    D <- D + 0.001
    i <- i + 1
}

ews <- do.call(rbind, ews)
k <- degree(g)

## First Plot
if(which_net == "me") {
    filename <- "./img/n-earlywarnings-sd.pdf"
} else if(which_net == "pa") {
    filename <- "./img/n-earlywarnings-sd-PA.pdf"
}

if(save_plot) {
    pdf(filename, height = 7, width = 21)
} else dev.new(width = 21, height = 7)

colorfun <- colorRamp(c("navy", "dodgerblue1"), space = "Lab")
nodecolor <- colorfun(k/max(k))
V(g)$color <- apply(nodecolor, 1, function(x) rgb(t(x), maxColorValue = 255))
par(mai = c(1, .8, .8, .8) + .02)
matplot(
    Ds, ews, type = "p", pch = 1, lty = 1, col = V(g)$color,
    xlab = "D", ylab = expression(sd(x[i]))
)
legend("left", pch = c(1, 1, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
       col = c("dodgerblue1", "navy", "darkorchid"), bty = "n",
       legend = c("Highest degree", "Lowest degree", "# Nodes 'At Risk'"))
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], type = "l", lwd = 2, lty = 1, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(n_atrisk[-length(n_atrisk)])))
mtext("# Nodes", side = 4, line = 3)
if(save_plot) dev.off()

## Second plot
if(which_net == "me") {
    filename <- "./img/n-earlywarnings-I.pdf"
} else if(which_net == "pa") {
    filename <- "./img/n-earlywarnings-I-PA.pdf"
}

if(save_plot) {
    pdf(filename, height = 7, width = 21)
} else dev.new(width = 21, height = 7)

par(mai = c(1, .8, .8, .8) + .02)
plot(Ds, Is, type = "p", lwd = 2, col = "red3",
     xlab = "D", ylab = "Moran's I")
legend("left", pch = c(1, NA), lty = c(NA, 1), lwd = c(2, 2),
       col = c("red3", "darkorchid"), bty = "n",
       legend = c("Moran's I", "# Nodes 'At Risk'"))
par(new = TRUE)
plot(Ds, n_atrisk[-length(n_atrisk)], type = "l", lwd = 2, lty = 1, col = "darkorchid",
     axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(n_atrisk[-length(n_atrisk)])))
mtext("# Nodes", side = 4, line = 3)
if(save_plot) dev.off()
##print(n_atrisk)
