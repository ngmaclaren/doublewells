An R package and simulation files in support of the double-wells portion of the network resiliency project.

The file `dakos-etal-replication.R` reproduces the basic early warning indicators for the harvest model system in Dakos et al. (2012) and the file `single-system.R` demonstrates their use in a single double-well system.

## Package Install

To install the "doublewells" package, download `doublewells_*.tar.gz` and install from the command line (in the directory where the tarball is located):

```
$ R CMD INSTALL doublewells_0.1.tar.gz # or the current version
```

In addition to the doublewells package, this project depends on the "earlywarnings" package by Dakos et al. (http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0041010, http://cran.r-project.org/web/packages/earlywarnings/index.html) and may in the future depend on the "deSolve" package by Soetaert et al. (doi:10.18637/jss.v033.i09, http://desolve.r-forge.r-project.org/). 
