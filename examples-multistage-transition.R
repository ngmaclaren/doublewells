## Code by Neil MacLaren 2/25/2022

library(igraph)
library(doublewells)
    
## To save the figure outside of this session
save_plots <- FALSE # TRUE
imgfile <- "./img/examples-multistage-transition.pdf"

## Run new sims with this seed or load a previous version from disk?
run_sims <- FALSE # TRUE
set.seed(123)
outfile <- "./data/examples-results.rda"

## Networks to use for examples
choices <- c("pref_attach", "dolphins")

if(run_sims) {
    data(list = choices) ## load the appropriate networks

    df_list <- list() ## storage list
    ##df_list_alt <- list()

    for(i in 1:length(choices)) {## do the same procedure with both networks
        g <- get(choices[i]) ## convenient renaming for the network object

        df_list[[i]] <- simulation(g) ## run the simulation and store the results
       ## df_list_alt[[i]] <- simulation(g, state_check = 75)
    }

    examples_results <- df_list ## rename results for saving purposes (because objects that have been saved with save() open as the name they were saved under
    save(examples_results, file = outfile)
} else {## if not running the simulations new
    load(outfile)
    df_list <- examples_results ## for convenience, give the results the same name as if sims were run
}

## Plot results
wd <- 10
ht <- 10
ewi <- "avgac"
colors <- c("sienna", "navy", "olivedrab")
ylims <- list(
    c(0, 110),
    c(0, 70)
)
if(save_plots) {
    pdf(imgfile, width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mfrow = c(2, 1), mar = c(5, 4, 4, 4)+.1, xpd = TRUE) #13
for(i in 1:length(choices)) {
    df <- df_list[[i]]$results
    columns <- colnames(df)[grep(ewi, colnames(df))]
    plot(df$Ds, df$n_lowerstate,
         ##type = "p", pch = 1,
         type = "l", lwd = 3, lty = 1,
         col = "darkorchid",
         xlim = range(df$Ds), ylim = ylims[[i]],
         xlab = "D", ylab = "# Nodes in Lower State")
    legend("topright", ncol = 4, bty = "n",
           ##pch = 1:4,
           lty = 1, lwd = 4,
           col = c("darkorchid", colors),
           legend = c("State", "All", "Lower", "Sentinel")
           )
    mtext(LETTERS[i], adj = 0.01, line = -1.5, cex = 1.5, font = 2)
    par(new = TRUE)
    plotsamples <- 1:nrow(df)
    matplot(df$Ds[plotsamples], df[plotsamples, columns],
            ##type = "p", pch = 2:4,
            type = "l", lwd = 2, lty = 1,
            col = colors,
            xlim = range(df$Ds), ylim = c(0, 1), axes = FALSE,
            xlab = "", ylab = "")
    axis(4, at = pretty(c(0, 1)))
    mtext("Average Autocorrelation",side = 4, cex = 1, font = 1, line = 3)
}
if(save_plots) dev.off()

## old code for testing
## cor(df_list[[1]]$Ds[df_list[[1]]$n_lowerstate == 100], df_list[[1]]$avgac_all[df_list[[1]]$n_lowerstate == 100], method = "kendall")
## subset(df_list[[1]], Ds <.3, select = c(n_lowerstate, Ds, avgac_all))
