## Code by Neil MacLaren, 3/24/2022

library(parallel)
library(nlme)
library(igraph)
library(doublewells)

                                        # Should the output be saved?
save_results <- FALSE # TRUE
save_plots <- FALSE # TRUE

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

                                        # The model for the maximum eigenvalue indicator...
maxeig <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
              data = rdf, subset = ewi == "maxeig")
                                        # the maximum standard deviation...
maxsd <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
             data = rdf, subset = ewi == "maxsd")
                                        # average standard deviation...
avgsd <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
             data = rdf, subset = ewi == "avgsd")
                                        # maximum autocorrelation
maxac <- lme(tau ~ nodeset, random = ~ 1 | run/network/nodeset,
             data = rdf, subset = ewi == "maxac")
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
print(summary(maxeig))
print(summary(maxsd))
print(summary(avgsd))
print(summary(maxac))
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
    colnames(df) <- c("lower", "est", "upper")
    rownames(df) <- c("all", "lower", "sentinel")
    df
}

plotvals <- list(
    maxeig = get_values(maxeig),
    maxsd = get_values(maxsd),
    avgsd = get_values(avgsd),
    maxac = get_values(maxac),
    avgac = get_values(avgac)
)

for(i in 1:length(plotvals)) {
    plotvals[[i]]$set <- rownames(plotvals[[i]])
    plotvals[[i]]$signal <- names(plotvals)[i]
}
plotvals <- do.call(rbind, plotvals)
rownames(plotvals) <- NULL

signals <- c("maxeig", "maxsd", "avgsd", "maxac", "avgac")
plotvals$signal <- factor(plotvals$signal, levels = signals)
nodesets <- c("all", "lower", "sentinel")
plotvals$set <- factor(plotvals$set, levels = nodesets)
plotvals$signal_ <- rev(as.numeric(plotvals$signal))
plotvals$set_ <- as.numeric(plotvals$set)

## colors <- c(
##     "sienna", # all
##     "navy", # lower
##     "olivedrab" # sentinel
## )
five_ten <- c(
    "#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE",
    "#882255", "#44AA99", "#999933", "#AA4499", "#AAAAAA"
)
pal <- five_ten
palette(pal)

                                        # This block makes the ladder plot based on the calculated
                                        # values from get_values().

plot_errorbars <- FALSE # TRUE
## yticklabels <- c("Maximum Autocorrelation", #"Average Standard Deviation",
##                  "Average Autocorrelation")
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
##legendlabels <- c("All Nodes", "Lower State Nodes", "Sentinel Nodes")
legendlabels <- c("All", "Lower State", "Sentinel")
ht <- 4
wd <- 8
##if(plot_errorbars) pchs <- 3 else pchs <- 1:3
pchcex <- 2; pchs <- 1:3
lwd <- 3
if(save_plots) {
    pdf("./img/kendall-corr-figure.pdf", width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mar = c(5, 8, 1, 1) + 0.1)
ylims <- range(plotvals$signal_) + c(-.5, .5)
xlims <- c(.7, 1)
plot(NULL, ylim = ylims, xlim = xlims, #xlim = c(min(lowers)*.9, max(uppers)*1.1),
     ylab = "", xlab = expression(tau), yaxt = "n", bty = "n",
     cex.lab = 1.75, cex.axis = 1.75)
axis(side = 2, tick = FALSE, labels = yticklabels, at = max(plotvals$signal_):1, las = 1,
     cex.axis = 1.75)
##ypos <- 2 + c((1/3), 0, -(1/3))
points(signal_ ~ est, data = plotvals, pch = pchs, lwd = lwd, col = plotvals$set_, cex = pchcex)
if(plot_errorbars) {
    segments(x0 = plotvals$lower, x1 = plotvals$upper, y0 = plotvals$signal_,
             lty = 1, lwd = lwd, col = plotvals$set_)
}
legend("bottomright", legend = legendlabels, col = 1:3, bty = "n", pch = 1:3, pt.lwd = lwd,
       cex = 1.25)
if(save_plots) dev.off()

##abline(h = 1.5, col = "gray", lwd = .5, lty = 1)
## ypos <- 1 + c((1/3), 0, -(1/3))
## points(plotvals[[2]]$ests, ypos, pch = pchs, col = colors, cex = pchcex)
## if(plot_errorbars) {
##     segments(x0 = plotvals[[2]]$lowers, x1 = plotvals[[2]]$uppers, y0 = ypos,
##              lty = 1, lwd = 2, col = colors)
## }
## legend("bottomleft", legend = legendlabels,
##        col = 1:length(legendlabels), bty = "n", pch = pchs)#, lty = 1, lwd = 1)#, pt.cex = 2)



## nodesets <- c("all", "lower", "sentinel")
