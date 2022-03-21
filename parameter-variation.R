## Code by Neil MacLaren 3/1/2022

library(igraph)
library(doublewells)

set.seed(123)

save_plots <- FALSE # TRUE
run_sims <- FALSE # TRUE

intensities <- c(.01, .05, .1, .5)
sample_sizes <- c(250, 150, 50, 25)
rs <- list(c(1, 3, 5), c(1, 2.5, 7), c(1, 4, 7), c(1, 5.5, 7))
stopifnot(length(intensities) == length(sample_sizes) & length(sample_sizes) == length(rs))
length.each <- length(intensities)

if(run_sims) {
    data(pref_attach)
    g <- pref_attach

    results <- list()

    for(i in 1:length.each) {
        print(intensities[i])
        results[[i]] <- simulation(g, s = intensities[i])
    }

    for(i in 1:length.each) {
        print(sample_sizes[i])
        results[[i + length.each]] <- simulation(g, nsamples = sample_sizes[i])
    }

    for(i in 1:length.each) {
        print(rs[[i]])
        results[[i + 2*length.each]] <- simulation(g, r = rs[[i]])
    }

    parameter_variation_results <- results
    save(parameter_variation_results, file = "./data/parameter-variation-results.rda")
} else {
    load("./data/parameter-variation-results.rda")
    results <- parameter_variation_results
}

if(!dir.exists("./r-out")) dir.create("./r-out")
sink("./r-out/parameter-variation.txt", append = FALSE, type = "output", split = TRUE)

## modify this to give a data frame with mean, sd, ll, and ul ?
for(i in 1:length(results)) {## this prints the MEAN Kendall correlation. I will need sd as well.
    if(i <= length.each) {
        print(paste0("Intensity: ", intensities[i]))
    } else if(i > length.each & i <= 2*length.each) {
        print(paste0("N Samples: ", sample_sizes[i - length.each]))
    } else print(paste0("r = ", rs[[i - 2*length.each]]))
    corr_results <- Kendall_correlations(results[[i]]$results)
    print("Means")
    print(corr_results$means)
    ##print("SDs")
    ##print(corr_results$sds)
}

sink()
