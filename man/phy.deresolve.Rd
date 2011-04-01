\name{phy.deresolve}
\alias{phy.deresolve}
\title{collapse phylogenetic nodes into multichotomies}
\description{collapses phylogenetic dichotomies within a temporal range into multichotomies}
\usage{
phy.deresolve(phy, time.range = c(0, 0), relative = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class \code{phylo}; see \code{\link[ape]{read.tree}}}
  \item{time.range}{a supplied temporal range, within which nodes will be collapsed}
  \item{relative}{whether relative branch lengths are used (e.g., root is 1; present is 0)}

}
\details{
If \code{time.range} is supplied as a single value, \code{phy.deresolve} will assume that the user intended a range whose minimum is 
zero (i.e., present; see \strong{Examples} below).  If no nodes are present within the supplied \code{time.range}, the tree will be returned unadulterated.}
\value{A phylogenetic tree of class \code{phylo} }
\author{Jonathan Eastman, based on \code{\link[ape]{di2multi}} by Emmanuel Paradis}
\examples{
data(sp.example)
attach(sp.example)

## An example with an internal temporal bin
# arbitrarily transform tree to a root age of 100
phy$edge.length=phy$edge.length*(100/max(branching.times(phy)))

# collapse nodes within the middle half of the tree's history
mulphy=phy.deresolve(phy=phy, time.range=c(0.05, 0.95))

# compare result with original tree
plot(mulphy, show.tip.label=FALSE)
mtext("multifurcations in the middle nine-tenths")
dev.new()
plot(phy, show.tip.label=FALSE)
mtext("original phylogeny")

## An example with time.range given by a single value
mulphy=phy.deresolve(phy=phy, time.range=0.25) # affects the most recent quarter of the tree
dev.new()
plot(mulphy, show.tip.label=FALSE)
mtext("multifurcations in the most recent quarter")

## An example showing effect of 'relative'
phy$edge.length=phy$edge.length*(100/max(branching.times(phy)))
dev.new()
plot(dd<-phy.deresolve(phy=phy, time.range=c(0,1), relative=FALSE))
mtext("relative = FALSE")
dev.new()
plot(rr<-phy.deresolve(phy=phy, time.range=c(0,1), relative=TRUE))
mtext("relative = TRUE")



}
