## Monday, July 18, 2022: DO NOT RUN. Needs to be updated based on final simulation runs.

## Code by Neil MacLaren, 3/24/2022

library(parallel)
library(nlme)
library(igraph)
library(doublewells)

                                        # Should the output be saved?
save_results <- FALSE # TRUE
save_plots <- TRUE # FALSE

                                        # Which networks should be analyzed?
                                        # This is the same list as that in network-variation-sims.R
                                        # `empiricals` is a vector of empirical networks stored as
                                        # part of this package, retrieved from the "networkdata" pkg.
networks <- c(
    "erdos_renyi", "er_islands", "barabasi_albert", "LFR", "powerlaw", "fitness",
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
for(i in 1:length(sim_results)) {
    for(j in 1:length(networks)) {
        sim_results[[i]][[j]]$results <- as.numeric(sim_results[[i]][[j]]$results)
    }
}

                                        # Monday, July 18, 2022
                                        # Choosing only those networks with more than 1 transition
## check_transitions <- function(result, result_names = networks) {
##     ## old version
##     transitions <- lapply(result, function(df) table(df$n_lowerstate))
##     names(transitions) <- result_names
##     tr <- transitions[order(names(transitions))]
##     ntr <- lapply(tr, function(x) x[which(x >= 15)])
##     lengths(ntr)
## }

for(i in 1:length(sim_results)) names(sim_results[[i]]) <- networks
sim_dfs <- unlist(sim_results, recursive = FALSE)
for(df in sim_dfs) df$results <- as.numeric(df$results)

transitions <- lapply(sim_dfs, function(df) table(df$n_lowerstate))
tr <- transitions[order(names(transitions))]
ntr <- lapply(tr, function(x) x[which(x >= 15)])
lntr <- lengths(ntr)
lntrdf <- data.frame(network = names(lntr), tr = lntr)
count_transitions <- aggregate(tr ~ network, data = lntrdf, FUN = min)
selected_networks <- count_transitions$network[which(count_transitions$tr > 1)]

## level 1 is run, and is a list of lists
## level 2 is network, and is a list of data frames
#### OLD: level 3 is [1]: df [2]: sentinel history
## collect_dfs <- function(run_results, ...) {
##     lapply(run_results, function(dl) dl[[1]])
## }

## ##sim_dfs <- lapply(sim_results, collect_dfs)
## count_transitions <- lapply(sim_dfs, check_transitions)
## transitions <- as.data.frame(do.call(rbind, count_transitions))
## colnames(transitions) <- names(count_transitions[[1]])
## mins <- apply(transitions, 2, min)
## reduced_transitions <- transitions[, which(mins >= 2)]
## selected_networks <- colnames(reduced_transitions)


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
            cr <- doublewells::Kendall_correlations(x)$means # [[1]]
        })
    }
)
stopCluster(threads)
## for(i in 1:length(corr_results)) {
##     for(j in 1:length(networks)) {#(corr_results[[1]])) {
##         ##corr_results[[i]][[j]]["run", 1] <- i
##         ##corr_results[[i]][[j]] <- rbind(corr_results[[i]][[j]], c(i))
##     }
## }

                                        # Store correlation results as a list of data frames
df_list <- lapply(corr_results, function(x) as.data.frame(t(do.call(cbind, x))))
                                        # Store the names of the early warning indicator variables
##varcols <- colnames(df_list[[1]])[-which(colnames(df_list[[1]]) %in% c("run", "network"))]
varcols <- colnames(df_list[[1]])[-which(colnames(df_list[[1]]) %in% c("results"))]
                                        # And update the data frames with metadata.
for(i in 1:length(df_list)) {
    df_list[[i]]$network <- networks
    ##df_list[[i]]$run <- i
}
                                        # Reshape for statistical analysis.
                                        # Previous output was a matrix/data frame of correlation
                                        # values for the different indicators. For analysis, a row
                                        # should contain the indicator value and the indicator label,
                                        # as well as the network label.
rdf_list <- lapply(df_list, reshape, varying = varcols, v.names = "tau",
                   timevar = "ewi_", idvar = "network", times = varcols,
                   #new.row.names = 1:10000,
                   direction = "long")
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
    ##rdf$run <- i
    rdf$rownames <- rownames(rdf)
    newcols <- strsplit(rdf$rownames, split = ".", fixed = TRUE)
    newcols <- do.call(rbind, newcols)
    rdf$network <- newcols[, 1]
    rdf <- rdf[, -grep("rownames", colnames(rdf))]
                                        # And the list of reshaped data frames updated.
    rdf_list[[i]] <- rdf
}
                                        # All the data frames are combined into one for analysis. 
rdf <- do.call(rbind, rdf_list)
rownames(rdf) <- 1:length(rownames(rdf))
                                        # Monday, July 18, 2022
                                        # Filter on `selected_networks
rdf <- rdf[rdf$network %in% selected_networks, ]

                                        # The models themselves.
                                        # These are linear mixed effects models, or repeated measures
                                        # models in ANOVA terms. The fixed effect is the repeated
                                        # measures: the node set (all nodes, only nodes in the lower
                                        # state, or only sentinel nodes). The repeated measures are
                                        # nested within networks, which are in turn nested within
                                        # simulation runs. There are 22 networks and 50 simulation
                                        # runs.

                                        # The model for the maximum eigenvalue indicator...
maxeig <- lme(tau ~ nodeset, random = ~ 1 | network,#/run,#/nodeset,
              data = rdf, subset = ewi == "maxeig")
                                        # the maximum standard deviation...
maxsd <- lme(tau ~ nodeset, random = ~ 1 | network,#/run,#/nodeset,
             data = rdf, subset = ewi == "maxsd")
                                        # average standard deviation...
avgsd <- lme(tau ~ nodeset, random = ~ 1 | network,#/run,#/nodeset,
             data = rdf, subset = ewi == "avgsd")
                                        # maximum autocorrelation
maxac <- lme(tau ~ nodeset, random = ~ 1 | network,#/run,#/nodeset,
             data = rdf, subset = ewi == "maxac")
                                        # and for the average autocorrelation coefficient.
avgac <- lme(tau ~ nodeset, random = ~ 1 | network,#/run,#/nodeset,
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

palette("R4")#conference val: "Tableau 10")

                                        # This block makes the ladder plot based on the calculated
                                        # values from get_values().

plot_errorbars <- TRUE # FALSE
yticklabels <- c("Dom. Eig.", "Max. SD", "Avg. SD", "Max. AC", "Avg. AC")
legendlabels <- c("All", "Lower State", "High Input")
ht <- 6
wd <- 8
if(plot_errorbars) {
    pchcex <- 1; pchs <- 1:3; lty <- 1; lwd <- 5
    plotvals$adjust <- plotvals$signal_ - ((plotvals$set_ - 2)/5)
} else {
    pchcex <- 2; pchs <- 1:3; lwd <- 3#4
}
if(save_plots) {
    pdf("./img/kendall-corr-figure.pdf", width = wd, height = ht)
    ##pdf("./complex-networks-2022/LatexAbstractTemplate/conference-kendall-fig.pdf",
      ##  width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mar = c(5, 8, 0, 1) + 0.1)
ylims <- range(plotvals$signal_) + c(-.33, .33)
xlims <- c(.6, .9)
plot(NULL, ylim = ylims, xlim = xlims, #xlim = c(min(lowers)*.9, max(uppers)*1.1),
     ylab = "", xlab = expression(tau), yaxt = "n", bty = "n",
     cex.lab = 1.75, cex.axis = 1.75)
axis(side = 2, tick = FALSE, labels = yticklabels, at = max(plotvals$signal_):1, las = 1,
     cex.axis = 1.75)
##ypos <- 2 + c((1/3), 0, -(1/3))
if(plot_errorbars) {
    points(adjust ~ est, data = plotvals, pch = pchs, lwd = lwd, col = plotvals$set_, cex = pchcex)
    ## segments(x0 = plotvals$est, y0 = plotvals$signal_ - adjust, y1 = plotvals$signal_ + adjust,
    ##          lwd = lwd, lty = lty, col = plotvals$set_)
    ## rect(xleft = plotvals$lower, xright = plotvals$upper,
    ##      ybottom = plotvals$signal_ - adjust, ytop = plotvals$signal_ + adjust,
    ##      density = NULL, col = NA, border = plotvals$set_, lty = lty, lwd = lwd)
    segments(x0 = plotvals$lower, x1 = plotvals$upper, y0 = plotvals$adjust,
             lty = 1, lwd = lwd, col = plotvals$set_)
    legend("topleft", legend = legendlabels, col = 1:3, bty = "n", lwd = lwd, lty = lty, pch = pchs)
    abline(h = (1:4) + 0.5, col = "gray", lwd = .5, lty = 2)
} else {
    points(signal_ ~ est, data = plotvals, pch = pchs, lwd = lwd, col = plotvals$set_, cex = pchcex)
    legend("bottomright", legend = legendlabels, col = 1:3, bty = "n", pch = 1:3, pt.lwd = lwd,
       cex = 1.25)
}
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




## #### Do not run
## ### Code to observe variation due to network
## hist(rdf$tau)
## summary(rdf$tau)

## dev.new(height = 7, width = 14)
## boxplot(tau ~ as.factor(network), data = rdf)

## plotdf <- subset(rdf, ewi == "avgac")
## plotdf$netno <- as.integer(factor(plotdf$network))
## plotdf$adjust <- as.numeric(factor(plotdf$nodeset))/3 - (1/3)

## svg("./img/network-var-vis-test.svg", height = 7, width = 14)
## palette("Polychrome 36")
## plot(NULL,
##      xlim = range(plotdf$netno), ylim = range(plotdf$tau),
##      xlab = "Network x NodeSet", ylab = "tau")
## points(x = plotdf$netno + plotdf$adjust,
##        y = plotdf$tau,
##        pch = plotdf$netno,
##        col = plotdf$netno)
## dev.off()



## ## Checking intuition on number of transitions
## dl <- sim_results[[1]]
## dfs <- vector("list", 22)
## for(i in 1:length(dl)) dfs[[i]] <- dl[[i]][[1]]
## transitions <- lapply(dfs, function(df) table(df$n_lowerstate))
## names(transitions) <- networks
## tr <- transitions[order(names(transitions))]
## lapply(tr, function(x) x[which(x >= 15)]) # yes, matches.
## ## Try a different sim run
## dl <- sim_results[[7]]
## dfs <- vector("list", 22)
## for(i in 1:length(dl)) dfs[[i]] <- dl[[i]][[1]]
## transitions <- lapply(dfs, function(df) table(df$n_lowerstate))
## names(transitions) <- networks
## tr <- transitions[order(names(transitions))]
## lapply(tr, function(x) x[which(x >= 15)])

## check_transitions <- function(result, result_names = networks) {
##     transitions <- lapply(result, function(df) table(df$n_lowerstate))
##     names(transitions) <- result_names
##     tr <- transitions[order(names(transitions))]
##     ntr <- lapply(tr, function(x) x[which(x >= 15)])
##     lengths(ntr)
## }

## ## level 1 is run
## ## level 2 is network
## ## level 3 is [1]: df [2]: sentinel history
## collect_dfs <- function(run_results, ...) {
##     lapply(run_results, function(dl) dl[[1]])
## }

## test1 <- lapply(sim_results, collect_dfs)
## test2 <- lapply(test1, check_transitions)
## test3 <- as.data.frame(do.call(rbind, test2))
## colnames(test3) <- names(test2[[1]])
## lapply(test3, summary)
## mins <- apply(test3, 2, min)
## test4 <- test3[, which(mins >= 2)]

## test4 <- subset(t(test3), . >= 2)
