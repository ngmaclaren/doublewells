## Code by Neil MacLaren 3/1/2022

library(igraph)
library(doublewells)

set.seed(123)

save_plots <- FALSE

networks <- c("max_entropy", "me_islands", "scale_free_2", "LFR",
              "weaverbirds", "dolphins", "jpr", "sn_auth", "jazz", "netsci")

data(list = networks)

results <- list()
for(i in 1:length(networks)) {
    net <- networks[i]
    print(net)
    g <- get(net)
    results[[i]] <- simulation(g, D.init = 0.01)
}

Kendall_correlations <- function(df) {
    columns <- colnames(df)[3:length(colnames(df))]
    dfsplit <- split(df, factor(df$n_lowerstate))

    kendalls <- lapply(dfsplit, function(x) {
        cor(x[, c("Ds", columns)], method = "kendall", use = "pairwise.complete.obs")[1, -1]
    })
    kendalls <- as.data.frame(do.call(rbind, kendalls))
    kendalls$n_lowerstate <- as.integer(rownames(kendalls))
    kendalls$n_steps <- sapply(dfsplit, nrow)

    kendalls <- kendalls[kendalls$n_steps > 15, ]

    ##kendalls <- kendalls[order(kendalls$n_steps), ]
    results <- list(
        means = as.matrix(round(colMeans(kendalls[, columns], na.rm = TRUE), 3)),
        sds = as.matrix(round(apply(kendalls[, columns], 2, sd, na.rm = TRUE), 3))
    )
    return(results)
}

if(!dir.exists("./r-out")) dir.create("./r-out")
sink("./r-out/network-variation.txt", append = FALSE, type = "output", split = TRUE)

for(i in 1:length(results)) {## this prints the MEAN Kendall correlation. I will need sd as well.
    print(networks[i])
    corr_results <- Kendall_correlations(results[[i]]$results)
    print("Means")
    print(corr_results$means)
    print("SDs")
    print(corr_results$sds)
}

sink()
