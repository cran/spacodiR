\name{phy.nodetimes}

\alias{phy.nodetimes}

\title{finding splitting times for nodes}

\description{\code{phy.nodetimes} is a function to collect all nodes (and times) that are within a particular time fraction of the total temporal extent of a phylogeny}

\usage{phy.nodetimes(phy, time.range = c(0, 0), proportion = TRUE)}

\arguments{
  \item{phy}{a phylogenetic tree of class \code{phylo}; see \code{\link[ape]{read.tree}}}
  \item{time.range}{the temporal slice of the tree defining the bounds for extracted nodes}
  \item{proportion}{Boolean; whether \code{time.range} is regarded as proportions of tree height}
}

\details{Given a temporal fraction of a phylogeny, this function will return all contained nodes (and 
branching times if desired). If used with \code{\link{spacodi.by.nodes}}, one could conduct significance testing of \code{Bst} through 
particular time slices of the phylogeny. Note: when using \code{proportion=TRUE}, '0' is regarded as the present, '1' as the root of the tree (e.g., \code{start.time=1.0} and 
\code{stop.time=0.8} gathers nodes for the basalmost or most rootward 20 percent of the temporal extent of the tree).}

\value{
A named vector of times (whose names are node labels; see \code{\link[ape]{nodelabels}} for the function used to label nodes).
}

\author{Jonathan Eastman}


\examples{
## an example of comparing observed and expected Bst at nodes ##
data(sp.example)
attach(sp.example)
s=spl
p=phy

# get nodes for a time-slice
n.t=phy.nodetimes(phy=p,time.range=c(0.25,0.75))

# get Bsts for nodes
o=spacodi.by.nodes(sp.plot=s, phy=p, obs.only=TRUE)[[1]]
obs=o[!is.na(match(o[,"node.ID"], names(n.t))),1]

# compile expected Bst-values from data reshuffling
spacodi.by.nodes(sp.plot=s, phy=p, n.rep=5, method="1s")[[2]]->exps
e=exps[!is.na(match(row.names(exps), names(n.t))),]
exp=unname(unlist(e))

# randomization test of observed and expected Bst
resamp.test(obs=obs, exp=exp, iter=100, two.tailed=TRUE)
}