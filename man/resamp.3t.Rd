\name{resamp.3t}

\alias{resamp.3t}

\title{randomizing a community phylogenetics matrix: '3t' of Hardy (2008)}
\description{\code{resamp.3t} is used for resampling data within a community dataset}


\usage{resamp.3t(obj, dmat=NULL)}

\arguments{
  \item{obj}{a community dataset in \pkg{spacodiR} format (see \code{\link{as.spacodi}})}
  \item{dmat}{an optional dataframe of distances between plots; row names and column names should be identical}

}

\details{A resampling procedure for a species-by-plots matrix, where observed abundances within species are shuffled to adjacent plots. 
This procedure thus assumes meaningful arrangement of plots in space. If a distance matrix is supplied, the likelihood of shuffling to 
a particular plot is proportional to the distance between the plots.     
}

\value{A shuffled dataset}

\references{
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
}

\author{Jonathan Eastman}

\seealso{\code{\link{spacodi.by.nodes}} for permutation tests of community diversity and additional randomization methods; 
see \code{\link{spacodi.calc}} for the main underlying function; \code{\link[picante]{randomizeMatrix}}}

\examples{
# define a community dataset
data(sp.example)
attach(sp.example)
spl

# define a distance matrix
foo <- matrix(runif((ncol(spl)->ss)^2,0,100),ss,ss)
foo[upper.tri(foo)] <- foo[lower.tri(foo)]
diag(foo) <- 0
dmat <- as.data.frame(foo)
row.names(dmat) <- names(spl)
names(dmat) <- row.names(dmat)


# shuffle dataset
resamp.3t(obj=spl, dmat=dmat) 
spl ## comparison with original
}
