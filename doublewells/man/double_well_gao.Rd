\name{double_well_gao}
\alias{double_well_gao}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Calculate the next state of a coupled double-well system
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Calculate the next state of a coupled double-well system using the Gao
  et al. approximation.
}
\usage{
double_well_gao(x, r1, r2, r3, D, A, beta_eff, dt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{r1}{
%%     ~~Describe \code{r1} here~~
}
  \item{r2}{
%%     ~~Describe \code{r2} here~~
}
  \item{r3}{
%%     ~~Describe \code{r3} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{beta_eff}{
%%     ~~Describe \code{beta_eff} here~~
}
  \item{dt}{
%%     ~~Describe \code{dt} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
