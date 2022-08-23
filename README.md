This repository contains an R package  (`doublewells`) and simulation files that support ["Early warnings for multistage transitions," by N. G. MacLaren, P. Kundu, and N. Masuda](https://arxiv.org/abs/2208.08961). These files have been tested on Arch (Manjaro) and Ubuntu Linux.

Data and code needed to reproduce our analyses are contained in this repository. Although most of the functions should be ready for use with related analysis, the purpose of the `doublewells` package is primarily to support the specific analyses needed for our project. If you have any difficulty reproducing our results or have questions or input, please don't hesitate to contact the authors or open a new issue in this repository. 

The networks included in the `doublewells` package and can be called by name, _e.g._:

```
library(igraph)
library(doublewells)

choices <- c("powerlaw", "dolphins")
data(list = choices)

dev.new()
par(mfrow = c(1, 2))
plot(powerlaw, main = "Power-law Network")
plot(dolphins, main = "Dolphins Network")
```

To reproduce 
- Figure 1, 2, and 4: Run `examples-sims.R` then `examples-analysis.R`.
- Figure 3 and Table 1: Run `network-variation-sims.R` then `network-variation-analysis.R`. Running time for `network-variation-sims.R` is quite long: approximately 1.1 hr per round of simulations, currently hard coded to 50 rounds. Running time was assessed on a desktop with four Intel CORE i3-3220 CPUs at 3.30 GHz and 8 GB of RAM. 
- Supplemental Figures: Run `parameter-variation.R`.

## Installation and Dependencies

To install `doublewells`, download `doublewells_*.tar.gz` and install from the command line (in the directory where the tarball is located):

```
$ R CMD INSTALL doublewells_0.1.tar.gz # or the current version
```

In addition to the `doublewells` package, this analysis depends on `igraph` (https://igraph.org/r/), `networkdata` (http://networkdata.schochastics.net/), and `parallel`.

We used the NetworkX [implementation](https://networkx.org/documentation/stable/reference/generated/networkx.generators.community.LFR_benchmark_graph.html) of the Lancichinetti-Fortunato-Radicchi model. The file `save-lfr-model.py` depends on NumPy and NetworkX. 
