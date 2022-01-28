\name{double_well}
\alias{double_well}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Calculate the next state for a single double-well system
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Calculates the next step in a double well simulation with no coupling
  and optional stress and noise.
}
\usage{
double_well(x, r1, r2, r3, dt, noise = NULL, stress = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    The current system state.
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
  \item{dt}{
    %%     ~~Describe \code{dt} here~~
    The size of the time step for integration.
  }
  \item{noise}{
    %%     ~~Describe \code{noise} here~~
    Should be a random number generator as a function call. For example
    `rnorm(1, 0, 10)`.
  }
  \item{stress}{
    %%     ~~Describe \code{stress} here~~
    A constant representing stress added to the system. 
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  This function implements the differential equation for the double-well
  system using the stochastic Euler method:

  \[Delta x = (-(x - r_1)(x - r_2)(x - r_3) + u + e) Delta t\]

  The state variable is `x`; `r1`, `r2`, and `r3` are parameters; `u` is
  additive stress; and epsilon is a noise term. For convenience, use
  `r1` < `r2` < `r3`. The `noise` argument if not NULL should be a
  function that generates a single random value.
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
  A numeric vector of length 1 representing the next state `x` of the system.
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
x <- 1; r1 <- 1; r2 <- 2; r3 <- 5; dt <- 0.001
T <- 100
state <- numeric(T)
for(t in 1:T) {
  state[t] <- x
  x <- double_well(x, r1, r2, r3, dt, noise = rnorm(1, 0, 10))
}
plot(1:T, state, type = "o")
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
