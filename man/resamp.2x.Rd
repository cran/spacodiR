\name{resamp.2x}

\alias{resamp.2x}

\title{randomizing a community phylogenetics matrix: '2x' of Hardy (2008)}
\description{\code{resamp.2x} is used for resampling data within a community dataset}

\usage{resamp.2x(obj, level=0.1)}

\arguments{
  \item{obj}{a community dataset in \pkg{spacodiR} format (see \code{\link{as.spacodi}})}
  \item{level}{a proportion specifying the extent of data shuffling}
}

\details{
  A resampling procedure for a species-by-plots matrix, based on Gotelli swapping. 
  Shuffles abundances within a pair of plots and for a pair of species.  The \code{level} defines 
  the degree of sampling, with larger values dictating a higher level of reshuffling.  For instance, 
  if \code{level = 0.4} and the dataset involves 5 species and 10 plots, a total of 20 (\code{0.4x5x10}) Gotelli-like swaps are performed.  The number of swaps is always rounded up such that at least one swap is performed (i.e., when \code{level=0}).
 
}
\value{A shuffled dataset}

\references{
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
  
  GOTELLI NJ. 2000. Null model analysis of species co-occurrence patterns. Ecology 81:2606-2621.
}
\author{Jonathan Eastman}

\seealso{\code{\link{spacodi.by.nodes}} for permutation tests of community diversity and additional randomization methods; 
see \code{\link{spacodi.calc}} for the main underlying function; \code{\link[picante]{randomizeMatrix}}}
\examples{
data(sp.example)
attach(sp.example)
spl

# shuffle dataset
resamp.2x(obj=spl, level=0.2) 
}
