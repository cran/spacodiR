\name{resamp.1a}

\alias{resamp.1a}

\title{randomizing a community phylogenetics matrix: '1a' of Hardy (2008)}
\description{\code{resamp.1a} is used for resampling data within a community dataset}

\usage{resamp.1a(obj, abund.class.ratio = 4)}

\arguments{
  \item{obj}{a community dataset in \pkg{spacodiR} format (see \code{\link{as.spacodi}})}
  \item{abund.class.ratio}{a ratio defining abundance classes}
}

\details{A resampling procedure for a species-by-plots matrix, where species are shuffled within abundance classes. 
Species are grouped into distinct abundance classes characterized by a fixed ratio: 
\code{abund.class.ratio} = maximal abundance / minimal abundance.  For instance, if \code{abund.class.ratio = 4},
the limits between abundance classes could be 1, 4, 16, ... .  Species are randomly permuted within each 
class, which maintains most of the abundance phylogenetic structure originally present in a dataset.  
}

\value{A shuffled dataset}

\references{
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
}

\author{Timothy Paine and Jonathan Eastman}

\seealso{\code{\link{spacodi.by.nodes}} for permutation tests of community diversity and additional randomization methods; 
see \code{\link{spacodi.calc}} for the main underlying function; \code{\link[picante]{randomizeMatrix}}}

\examples{
data(sp.example)
attach(sp.example)
spl

# shuffle dataset
resamp.1a(obj=spl, abund.class.ratio=3) 
}
