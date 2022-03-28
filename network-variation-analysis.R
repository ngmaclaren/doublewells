## Code by Neil MacLaren, 3/24/2022

library(parallel)
library(nlme)
library(igraph)
library(doublewells)

save_results <- FALSE # TRUE
save_plots <- FALSE # TRUE

networks <- c(
    "max_entropy", "me_islands", "pref_attach", "scale_free_2", "LFR",
    empiricals
)

dataloc <- "./data/sims/"
datafiles <- paste0(dataloc, list.files(dataloc))

for(datafile in datafiles) load(datafile)

dls <- gsub("-", "_", list.files(dataloc))
dls <- gsub(".rda", "", dls)

                                        # Each sim output object is a list of 22
                                        # lists, each of which has a "results" df
                                        # inside it.
## dfs <- lapply(dls, function(dl) lapply(get(dl), function(sr) sr$results))
## for(i in 1:length(dfs)) {
##     for(j in 1:length(networks)) {
##         dfs[[i]][[j]]$network <- networks[j]
##     }
## }

##sim_results <- do.call(rbind, dfs)
sim_results <- lapply(dls, get)

## corr_results <- lapply(sim_results, function(results) {
##     lapply(results, function(x) Kendall_correlations(x[[1]])$means)
## })

nthreads <- detectCores() - 1
threads <- makeCluster(nthreads)
corr_results <- clusterApply(
    cl = threads, x = sim_results,
    fun =  function(results) {
        lapply(results, function(x) {
            cr <- doublewells::Kendall_correlations(x[[1]])$means
        })
    }
)
stopCluster(threads)

df_list <- lapply(corr_results, function(x) as.data.frame(t(do.call(cbind, x))))

varcols <- colnames(df_list[[1]])
for(i in 1:length(df_list)) {
    df_list[[i]]$network <- networks
    df_list[[i]]$run <- i
}

rdf_list <- lapply(df_list, reshape, varying = varcols, v.names = "tau",
                   timevar = "ewi_", idvar = "network", times = varcols,
                   new.row.names = 1:10000, direction = "long")

for(i in 1:length(rdf_list)) {
    rdf <- rdf_list[[i]]
    newcols <-  strsplit(rdf$ewi_, "_")
    newcols <- as.data.frame(do.call(rbind, newcols))
    colnames(newcols) <- c("ewi", "nodeset")
    rdf <- cbind(rdf, newcols)
    rdf$nodeset <- factor(rdf$nodeset, levels = c("all", "lower", "sentinel"),
                          ordered = FALSE)
    rdf_list[[i]] <- rdf
}

rdf <- do.call(rbind, rdf_list)

avgsd <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
             data = rdf, subset = ewi == "avgsd")
avgac <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
             data = rdf, subset = ewi == "avgac")

if(save_results) {
    if(!dir.exists("./r-out")) dir.create("./r-out")
    sink("./r-out/network-variation-results.txt", append = FALSE, type = "output",
         split = TRUE)
}

print(summary(avgsd))
print(summary(avgac))

if(save_results) sink()

get_values <- function(model) {
    ints <- intervals(model, which = "fixed")
    ests <- ints$fixed[, 2]
    ests[2:3] <- ests[1] + ests[2:3]
    lowers <- ints$fixed[, 1]
    lowers[2:3] <- ests[1] + lowers[2:3]
    uppers <- ints$fixed[, 3]
    uppers[2:3] <- ests[1] + uppers[2:3]

    df <- as.data.frame(t(rbind(lowers, ests, uppers)))
    ## df <- df[c(2, 1, 3), ] # needed if lower is the ref in the model
    rownames(df) <- c("all", "lower", "sentinel")
    df
}

plotvals <- list(
    avgsd = get_values(avgsd),
    avgac = get_values(avgac)
)

colors <- c(
    "sienna", # all
    "navy", # lower
    "olivedrab" # sentinel
)

yticklabels <- c("Average Standard Deviation", "Average Autocorrelation")
legendlabels <- c("All Nodes", "Lower State Nodes", "Sentinel Nodes")
ht <- 7
wd <- 8
if(save_plots) {
    pdf("./img/kendall-corr-figure.pdf", width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mar = c(5, 12, 1, 1) + 0.1)
ylims <- c(.5, 2.5)
plot(NULL, ylim = ylims, xlim = c(0.45, .95), #xlim = c(min(lowers)*.9, max(uppers)*1.1),
     ylab = "", xlab = expression(tau), yaxt = "n")
axis(side = 2, tick = FALSE, labels = yticklabels, at = 2:1, las = 1)
ypos <- 2 + c((1/3), 0, -(1/3))
points(plotvals[[1]]$ests, ypos, pch = 16:18, col = colors, cex = 1)
segments(x0 = plotvals[[1]]$lowers, x1 = plotvals[[1]]$uppers, y0 = ypos,
         lty = 1, lwd = 2, col = colors)
abline(h = 1.5, col = "gray", lwd = .5, lty = 1)
ypos <- 1 + c((1/3), 0, -(1/3))
points(plotvals[[2]]$ests, ypos, pch = 16:18, col = colors, cex = 1)
segments(x0 = plotvals[[2]]$lowers, x1 = plotvals[[2]]$uppers, y0 = ypos,
         lty = 1, lwd = 2, col = colors)
legend("bottomleft", legend = legendlabels,
       col = colors, bty = "n", pch = 16:18, lty = 1, lwd = 1)#, pt.cex = 2)
if(save_plots) dev.off()
