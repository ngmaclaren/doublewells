## Code by Neil MacLaren 3/1/2022

library(igraph)
library(doublewells)

set.seed(123)

save_plots <- FALSE

data(pref_attach)
g <- pref_attach

results <- list()
intensities <- c(.001, .005, .01)
sample_sizes <- c(250, 150, 50)

for(i in 1:length(intensities)) {
    print(intensities[i])
    results[[i]] <- simulation(g, s = intensities[i])
}

for(i in (1:length(sample_sizes) + length(intensities))) {
    print(sample_sizes[i - length(intensities)])
    results[[i]] <- simulation(g, nsamples = sample_sizes[i - length(intensities)])
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
sink("./r-out/parameter-variation.txt", append = FALSE, type = "output", split = TRUE)

for(i in 1:length(results)) {## this prints the MEAN Kendall correlation. I will need sd as well.
    if(i <= 3) {
        print(paste0("Intensity: ", intensities[i]))
        print(paste0("N Samples: ", 250))
    } else if(i > 3) {
        j <- i - 3
        print(paste0("Intensity: ", 0.005))
        print(paste0("N Samples: ", sample_sizes[j]))
    }
    corr_results <- Kendall_correlations(results[[i]]$results)
    print("Means")
    print(corr_results$means)
    print("SDs")
    print(corr_results$sds)
}

sink()
