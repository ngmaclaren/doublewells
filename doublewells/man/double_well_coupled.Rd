\name{double_well_coupled}
\alias{double_well_coupled}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Calculate the next state of a coupled double-well system
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Calculates the next step in a double well simulation with
  network-based coupling (specified in an adjacency matrix) and optional
  stress and noise.
}
\usage{
double_well_coupled(x, r1, r2, r3, D, A, dt, noise = NULL, stress = rep(0, length(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    A numeric vector (row vector) of current states for each node in A.
}
  \item{r1}{
    %%     ~~Describe \code{r1} here~~
    Sets the vicinity of the lower (stable) equilibrium.
  }
  \item{r2}{
    %%     ~~Describe \code{r2} here~~
    Sets the vicinity of the middle (unstable) equilibrium.
  }
  \item{r3}{
    %%     ~~Describe \code{r3} here~~
    Sets the vicinity of the upper (stable) equilibrium.
  }
  \item{D}{
    %%     ~~Describe \code{D} here~~
    The coupling strength.
  }
  \item{A}{
    %%     ~~Describe \code{A} here~~
    The adjacency matrix.
  }
  \item{dt}{
    %%     ~~Describe \code{dt} here~~
    The size of the time step for integration.
  }
  \item{noise}{
    %%     ~~Describe \code{noise} here~~
    Should be a random number generator as a function call. For example
    `rnorm(length(x), 0, 10)`.
  }
  \item{stress}{
    %%     ~~Describe \code{stress} here~~
    A constant representing stress added to the system. Should be a
    numeric vector (row vector) with a value for each `x_i` in the
    system. For example, if only adding stress to `x_1` in a three node
    system, `c(1, 0, 0)`.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  This function implements the differential equation for the double-well
  system using the stochastic Euler method:
  
  \[Delta x_i = (-(x_i - r_1)(x_i - r_2)(x_i - r_3) + D \sum_{j=1}^n
  Ax_j + u_i + e_i) Delta t\]

  The state variable is `x`; `r1`, `r2`, and `r3` are parameters; `D` is
  the coupling strength, assumed to be the same across the system; `A`
  is an adjacency matrix; `u` is additive stress; and `\epsilon` is a
  noise term. For convenience, use `r1` < `r2` < `r3`. The `noise`
  argument if not NULL should be a function that generates a single
  random value.
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
  A numeric vector of length `n` representing the next state `x` of the
  system for all `x_i`.
}
\references{
  %% ~put references to the literature/web site here ~
  Kundu et al. (2021). https://arxiv.org/abs/2110.07140
}
\author{
  %%  ~~who you are~~
  Neil G. MacLaren
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x <- rep(1, 3); r1 <- 1; r2 <- 2; r3 <- 5; dt <- 0.001; D <- 0.1
A <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), byrow = TRUE, nrow = 3)
T <- 100
state <- matrix(0, ncol = length(x), nrow = T)
for(t in 1:T) {
  state[t, ] <- x
  x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise = rnorm(length(x), 0, 10))
}
matplot(1:T, state, type = "o")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
