\name{match.spacodi.data}

\alias{match.spacodi.data}

\title{pruning a tree to match datasets used in SPACoDi}

\description{
 \code{match.spacodi.data} is a utility for automating the process of pruning datasets for use in \pkg{spacodiR}. Note that nearly all \pkg{spacodiR} functions require that
 trees have present (at least) all species sampled in the community dataset.}

\usage{match.spacodi.data(sp.plot, phy=NULL, sp.traits=NULL, prune=TRUE, verbose=FALSE)}

\arguments{
  \item{sp.plot}{a community dataset in \pkg{spacodiR} format (see \code{\link{as.spacodi}})}
  \item{phy}{a phylogenetic tree of class \code{phylo}; see \code{\link[ape]{read.tree}}}
  \item{sp.traits}{a species-by-trait(s) dataframe}
  \item{prune}{Boolean; whether plots with fewer than two species sampled are pruned}
  \item{verbose}{Boolean; whether dropped plots are written to the console}
}

\value{
A list of pruned dataset(s): 
\itemize{
    \item \code{sp.plot}
	
    \item \code{sp.tree}
	
	\item \code{sp.traits}
	}
	Note: if \code{prune=TRUE}, the list of plots pruned will be printed to the console.
}

\author{Jonathan Eastman}

\examples{
# load a species-by-plots matrix, along with a tree
data(sp.example)
attributes(sp.example)
attach(sp.example)
spl
phy
trt

# prune out undersampled plots
spl[,2]=0
match.spacodi.data(spl) -> sp.plot.new
sp.plot.new


# match datasets where sp.traits is smaller than the remainder
match.spacodi.data(sp.plot=spl, phy=phy, sp.traits=trt[1:6,])
							
}


