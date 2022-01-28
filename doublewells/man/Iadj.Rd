\name{Iadj}
\alias{Iadj}

\title{
  Calculate a variant of Moran's I. 
}
\description{
  Calculate a variant of Moran's I that removes the expected value of
  x_i. A is a 2D array of connection weights (all w_ij ∈ {0, 1} for
  now), t is the chosen time, X is a 2D array where the rows are t_t and
  the columns are x_i.
}
\usage{
Iadj(X, A, t = NULL, nsteps = NULL, times = NULL)
}
\arguments{
  \item{X}{
    X is a 2D array where the rows are t_t and the columns are x_i.
}
  \item{A}{
    A is a 2D array of connection weights (all w_ij ∈ {0, 1} for now).
  }
  \item{t}{
    A time at which to calculate the index. Sets the analysis point
    furthest forward in time.
  }
  \item{nsteps}{
    From t, how many steps to look backwards in time. 
  }
  \item{times}{
    Alternately, provide a vector of times at which to sample X.
  }
}
\details{
  The function calculates a variant of Moran's I that removes the effect
  of variance in E(x_i). 
}
\value{
  A numeric vector of length one giving the value of the index at time t.
}
\references{
  %% ~put references to the literature/web site here ~
  
}
\author{
  %%  ~~who you are~~
  Neil G. MacLaren and Naoki Masuda
}
\note{
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%% x <- rep(1, 3); r1 <- 1; r2 <- 2; r3 <- 5; dt <- 0.001; D <- 0.1
%% A <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), byrow = TRUE, nrow = 3)
%% T <- 100
%% state <- matrix(0, ncol = length(x), nrow = T)
%% for(t in 1:T) {
%%   state[t, ] <- x
%%   x <- double_well_coupled(x, r1, r2, r3, D, A, dt, noise = rnorm(length(x), 0, 10))
%% }
%% matplot(1:T, state, type = "o")
%% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
}
