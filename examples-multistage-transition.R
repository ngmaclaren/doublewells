## Code by Neil MacLaren 2/25/2022

library(igraph)
library(doublewells)
    
                                        # To save the figure outside of this session
save_plots <- FALSE # TRUE
imgfile <- "./img/examples-multistage-transition.pdf"

                                        # Run new sims with this seed or load a previous version from
                                        # disk?
run_sims <- FALSE # TRUE
set.seed(123)
outfile <- "./data/examples-results.rda"
                                        # Save results for tables?
save_results <- FALSE # TRUE
r_outfile <- "./r-out/examples-results.txt"

                                        # Networks to use for examples
choices <- c("pref_attach", "dolphins")

if(run_sims) {
                                        # load the appropriate networks, which are part of the
                                        # doublewells package
    data(list = choices) 
                                        # Store results for both networks
    df_list <- list() 
    for(i in 1:length(choices)) {
                                        # Relabel for convenience
        g <- get(choices[i]) 
                                        # Run the simulation with default settings and store the
                                        # results.
        df_list[[i]] <- simulation(g)#, check_alts = TRUE)
    }
                                        # rename results for saving purposes (because objects that
                                        # have been saved with save() open as the name they were saved
                                        # under.
    examples_results <- df_list 
    save(examples_results, file = outfile)
} else {
                                        # if not running the simulations new
    load(outfile)
                                        # for convenience, give the results the same name as if sims
                                        # were run
    df_list <- examples_results 
}

                                        # Plot results
wd <- 10
ht <- 10
ewi <- "avgac"
colors <- rev( c("sienna", "navy", "olivedrab") )#, "cyan", "chartreuse") )
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
    df$p_lowerstate <- round((df$n_lowerstate/max(df$n_lowerstate))*100)
    columns <- rev( colnames(df)[grep(ewi, colnames(df))] )
    plot(df$Ds, df$p_lowerstate,#df$n_lowerstate,
         type = "l", lwd = 3, lty = 1,
         col = "darkorchid",
         xlim = range(df$Ds), ylim = c(0, 100),#ylims[[i]],
         xlab = expression(italic(D)), ylab = "% Nodes in Lower State")#"D"
    legend("topright", ncol = 4, bty = "n",
           lty = 1, lwd = 4,
           col = c("darkorchid", rev(colors)),
           legend = c("State", "All", "Lower", "Sentinel")#, "Upper", "Anti-Sentinel")
           )
    mtext(LETTERS[i], adj = 0.01, line = -1.5, cex = 1.5, font = 2)
    par(new = TRUE)
    plotsamples <- 1:nrow(df)
    matplot(df$Ds[plotsamples], df[plotsamples, columns],
            type = "l", lwd = 2, lty = 1,
            col = colors,
            xlim = range(df$Ds), ylim = c(0, 1), axes = FALSE,
            xlab = "", ylab = "")
    axis(4, at = pretty(c(0, 1)))
    mtext("Average Autocorrelation",side = 4, cex = 1, font = 1, line = 3)
}
if(save_plots) dev.off()

if(save_results) {
    sink(r_outfile, append = FALSE, type = "output", split = TRUE)
    for(i in 1:length(choices)) {
        df <- df_list[[i]]$results
        corr_results <- Kendall_correlations(df)
        print(choices[i])
        print(corr_results$means)
    }
}



## old code for testing
## cor(df_list[[1]]$Ds[df_list[[1]]$n_lowerstate == 100], df_list[[1]]$avgac_all[df_list[[1]]$n_lowerstate == 100], method = "kendall")
## subset(df_list[[1]], Ds <.3, select = c(n_lowerstate, Ds, avgac_all))
