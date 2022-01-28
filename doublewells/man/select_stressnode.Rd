\name{select_stressnode}
\alias{select_stressnode}

\title{
  Choose a node in a graph.
}
\description{
  Given a graph, g, choose a node according to one of a set of
  conditions. In the future, will add the capability to add stress to
  a set of nodes following the same decision rules.
}
\usage{
select_stressnode(g, add_stress_to = NULL)
}
\arguments{
  \item{g}{
    A graph.
}
  \item{add_stress_to}{
    One of "high", "low", "highest", "lowest", or NULL. See below.
  }
}
\details{
  Select a node in a graph according to a condition. The conditions are:
  (1) "high", choose a node with degree in the top quintile of degrees
  in the graph, (2) "low", same, but for the lowest quintile; (3)
  "highest" and (4) "lowest" choose the highest or lowest degree,
  respectively; and (5) NULL choose a random node. In all cases, if more
  than one node matches the decision rule, a random node from the set of
  nodes meeting the criteria is returned.
}
\value{
  An igraph node.
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
