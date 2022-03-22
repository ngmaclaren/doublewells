## Code by Neil MacLaren 2/25/2022

library(igraph)
library(doublewells)
    
set.seed(123)

save_plots <- FALSE # TRUE
run_sims <- FALSE # TRUE

choices <- c("pref_attach", "dolphins") #"scale_free_3"

if(run_sims) {
    data(list = choices)

    graphlist <- list()
    df_list <- list()
    kendalls_list <- list()

    for(k in 1:length(choices)) {
        g <- get(choices[k])
        graphlist[[k]] <- g

        df_list[[k]] <- simulation(g)#, TU = 75, check_equil = TRUE)
        ## kendalls_list[[k]] <- Kendall_correlations()
    }

    examples_results <- df_list
    save(examples_results, file = "./data/examples-results.rda")
} else {
    load("./data/examples-results.rda")
    df_list <- examples_results
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
    pdf("./img/examples-multistage-transition.pdf", width = wd, height = ht)
} else dev.new(width = wd, height = ht)
par(mfrow = c(2, 1), mar = c(5, 4, 4, 4)+.1, xpd = TRUE) #13
for(k in 1:length(choices)) {
    df <- df_list[[k]]$results
    columns <- colnames(df)[grep(ewi, colnames(df))]
    plot(df$Ds, df$n_lowerstate,
         ##type = "p", pch = 1,
         type = "l", lwd = 3, lty = 1,
         col = "darkorchid",
         xlim = range(df$Ds), ylim = ylims[[k]],
         xlab = "D", ylab = "# Nodes in Lower State")
    legend("topright", ncol = 4, bty = "n",
           ##pch = 1:4,
           lty = 1, lwd = 4,
           col = c("darkorchid", colors),
           legend = c("State", "All", "Lower", "Sentinel")
           )
    mtext(LETTERS[k], adj = 0.01, line = -1.5, cex = 1.5, font = 2)
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

## cor(df_list[[1]]$Ds[df_list[[1]]$n_lowerstate == 100], df_list[[1]]$avgac_all[df_list[[1]]$n_lowerstate == 100], method = "kendall")
## subset(df_list[[1]], Ds <.3, select = c(n_lowerstate, Ds, avgac_all))
