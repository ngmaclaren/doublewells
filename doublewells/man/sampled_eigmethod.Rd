\name{sampled_eigmethod}
\alias{sampled_eigmethod}

\title{
  Calculate an early warning indicator from a covariance matrix
}
\description{
  Calculate an early warning indicator in a multivariate system that can
  be sampled over time. A covariance matrix of the x_i is built from the
  values x_{i,t} and the dominant eigenvalue of that covariance matrix
  is returned.
}
\usage{
sampled_eigmethod(X, A, samples)
}
\arguments{
  \item{X}{
    X is a 2D array where the rows are t_t and the columns are x_i.
}
  \item{A}{
    A is a 2D array of connection weights (all w_ij ∈ {0, 1} for now).
  }
  \item{times}{
    A vector of times at which to sample X.
  }
}
\details{
  Calculate the dominant eigenvalue of a covariance matrix in a
  multivariate system evolving over time. Actually, X could be any
  matrix and `samples` is just a vector of row indices. The function
  takes X at `samples` rows, calculates a covariance matrix, and returns
  the dominant eigenvalue of the covariance matrix.
}
\value{
  A numeric vector of length one giving the value of the index at time t.
}
\references{
  %% ~put references to the literature/web site here ~
  
}
\author{
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
