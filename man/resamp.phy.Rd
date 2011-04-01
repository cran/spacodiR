\name{resamp.phy}
\alias{resamp.phy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{partial phylogeny randomization for tips}
\description{\code{resamp.phy} is used to shuffle tips subtended by a set of internal nodes determining by divergence time or specified by the user}

\usage{resamp.phy(phy, node = NULL, time.threshold = 1, proportion = TRUE)}


\arguments{
  \item{phy}{a phylogenetic tree of class \code{phylo}}
  \item{node}{numeric value(s), specifying the internal node(s) whose tips to shuffle; see \code{\link[ape]{nodelabels}}}
  \item{time.threshold}{either absolute or relative divergence times, used to define clades within which tips will be shuffled}
  \item{proportion}{if \code{TRUE}, divergence times are treated as relative (\code{0: present}; \code{1: root})}
}

\details{Either a numeric vector is supplied for \code{node} or a \code{time.threshold}.  If given a set of \code{node}s, this 
function will naively shuffle tips descended from the nodes in the order supplied (without regard to whether any internal node 
in the vector is a descendant of any other node in the \code{node} vector).  If given a \code{time.threshold}, tips will be reshuffled 
within non-nested clades that have a rootmost node that occurs within the range \code{[0, time.threshold]}.  Note that regard to 
absolute divergence times can be enforced with \code{proportion=FALSE}. Note further that \code{resamp.phy(phy=phy, node=NULL, node.threshold=1, proportion=TRUE)} 
achieves the same effect as \code{\link{resamp.1s}}.}


\value{a phylogenetic tree whose tips have been shuffled (without any modification of the underlying topology)}

\references{
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
}

\author{Jonathan Eastman}

\seealso{see \code{\link[ape]{nodelabels}} for values given for internal nodes and how these are retrieved}

\examples{
# load a tree
data(sp.example)
attach(sp.example)
plot(phy)
mtext("original phylogeny")

# reshuffle within a time range
time=1/3
bb=branching.times(phy)
bb=bb/max(bb)
nodes=(Ntip(phy)+1):max(phy$edge)
nodes[bb<=time]=1

dev.new()
plot(resamp.phy(phy, time.threshold=time, proportion=TRUE))
mtext("reshuffled phylogeny showing affected nodes")
nodelabels(cex=ifelse(nodes==1, 2, NA), col=ifelse(nodes==1, 1, NA), pch=19)

}
