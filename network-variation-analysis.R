## Code by Neil MacLaren, 3/24/2022

library(parallel)
library(nlme)
library(igraph)
library(doublewells)

                                        # Should the output be saved?
save_results <- TRUE # FALSE
save_plots <- TRUE # FALSE

                                        # Which networks should be analyzed?
                                        # This is the same list as that in network-variation-sims.R
                                        # `empiricals` is a vector of empirical networks stored as
                                        # part of this package, retrieved from the "networkdata" pkg.
networks <- c(
    "max_entropy", "me_islands", "pref_attach", "scale_free_2", "LFR",
    empiricals
)

                                        # Collect a vector of file names to load...
dataloc <- "./data/sims/"
datafiles <- paste0(dataloc, list.files(dataloc))
                                        # make those files available in this session...
for(datafile in datafiles) load(datafile)
                                        # make a vector of corresponding variable names...
dls <- gsub("-", "_", list.files(dataloc))
dls <- gsub(".rda", "", dls)
                                        # and store the sim results in a list.
sim_results <- lapply(dls, get)

                                        # Using parallel processing, calculate Kendall correlations
                                        # for all simulations. The Kendall_correlations() function is
                                        # documented as part of the doublewells package. Variables to
                                        # be correlated are the bifurcation parameter, D, and several
                                        # early warning indicators.
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
                                        # Store correlation results as a list of data frames
df_list <- lapply(corr_results, function(x) as.data.frame(t(do.call(cbind, x))))
                                        # Store the names of the early warning indicator variables
varcols <- colnames(df_list[[1]])
                                        # And update the data frames with metadata.
for(i in 1:length(df_list)) {
    df_list[[i]]$network <- networks
    df_list[[i]]$run <- i
}
                                        # Reshape for statistical analysis.
                                        # Previous output was a matrix/data frame of correlation values
                                        # for the different indicators. For analysis, a row should
                                        # contain the indicator value and the indicator label, as well
                                        # as the network label.
rdf_list <- lapply(df_list, reshape, varying = varcols, v.names = "tau",
                   timevar = "ewi_", idvar = "network", times = varcols,
                   new.row.names = 1:10000, direction = "long")
                                        # This loop splits the `varcols` variable names for analysis
                                        # convenience. Specifically, it separates the variables into
                                        # a label for the indicator (avg SD, avg AC, etc.) and a label
                                        # for the node set (all, lower, sentinel). 
for(i in 1:length(rdf_list)) {
    rdf <- rdf_list[[i]]
    newcols <-  strsplit(rdf$ewi_, "_")
    newcols <- as.data.frame(do.call(rbind, newcols))
    colnames(newcols) <- c("ewi", "nodeset")
    rdf <- cbind(rdf, newcols)
    rdf$nodeset <- factor(rdf$nodeset, levels = c("all", "lower", "sentinel"),
                          ordered = FALSE)
                                        # Metadata is also added. 
    rdf$run <- i
                                        # And the list of reshaped data frames updated.
    rdf_list[[i]] <- rdf
}
                                        # All the data frames are combined into one for analysis. 
rdf <- do.call(rbind, rdf_list)

                                        # The models themselves.
                                        # These are linear mixed effects models, or repeated measures
                                        # models in ANOVA terms. The fixed effect is the repeated
                                        # measures: the node set (all nodes, only nodes in the lower
                                        # state, or only sentinel nodes). The repeated measures are
                                        # nested within networks, which are in turn nested within
                                        # simulation runs. There are 22 networks and 50 simulation
                                        # runs.

                                        # The model for the average standard deviation indicator
avgsd <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
             data = rdf, subset = ewi == "avgsd")
                                        # and for the average autocorrelation coefficient.
avgac <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
             data = rdf, subset = ewi == "avgac")

                                        # If saving results, open the connection
if(save_results) {
    if(!dir.exists("./r-out")) dir.create("./r-out")
    sink("./r-out/network-variation-results.txt", append = FALSE, type = "output",
         split = TRUE)
}
                                        # Print the results to the console (and connection)
print(summary(avgsd))
print(summary(avgac))
                                        # and turn off the connection.
if(save_results) sink()

                                        # Below here the code makes a ladder/coefficient plot to
                                        # visualize the above results.

                                        # This function calculates the average effect for each node
                                        # set from zero, rather than with respect to the reference
                                        # category for the coefficents.
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
                                        # This block makes the ladder plot based on the calculated
                                        # values from get_values().
plot_errorbars <- FALSE
yticklabels <- c("Average Standard Deviation", "Average Autocorrelation")
legendlabels <- c("All Nodes", "Lower State Nodes", "Sentinel Nodes")
ht <- 7
wd <- 8
if(plot_errorbars) pchs <- 4 else pchs <- 15:17
pchcex <- 2
if(save_plots) {
    pdf("./img/kendall-corr-figure.pdf", width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mar = c(5, 12, 1, 1) + 0.1)
ylims <- c(.5, 2.5)
plot(NULL, ylim = ylims, xlim = c(0.45, .95), #xlim = c(min(lowers)*.9, max(uppers)*1.1),
     ylab = "", xlab = expression(tau), yaxt = "n")
axis(side = 2, tick = FALSE, labels = yticklabels, at = 2:1, las = 1)
ypos <- 2 + c((1/3), 0, -(1/3))
points(plotvals[[1]]$ests, ypos, pch = pchs, col = colors, cex = pchcex)
if(plot_errorbars) {
    segments(x0 = plotvals[[1]]$lowers, x1 = plotvals[[1]]$uppers, y0 = ypos,
             lty = 1, lwd = 2, col = colors)
}
abline(h = 1.5, col = "gray", lwd = .5, lty = 1)
ypos <- 1 + c((1/3), 0, -(1/3))
points(plotvals[[2]]$ests, ypos, pch = pchs, col = colors, cex = pchcex)
if(plot_errorbars) {
    segments(x0 = plotvals[[2]]$lowers, x1 = plotvals[[2]]$uppers, y0 = ypos,
             lty = 1, lwd = 2, col = colors)
}
legend("bottomleft", legend = legendlabels,
       col = colors, bty = "n", pch = pchs)#, lty = 1, lwd = 1)#, pt.cex = 2)
if(save_plots) dev.off()
