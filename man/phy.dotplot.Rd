\name{phy.dotplot}
\alias{phy.dotplot}
\title{EXPERIMENTAL: plotting members of communities (or plots) on a phylogeny}
\description{\code{phy.dotplot} is used to visually display species membership in communities, generating one plot per community}
\usage{phy.dotplot(sp.plot, phy, edge.width = 0.2, lab.adj = c(0, 0), tips.adj = c(0.5, 3), tips.cex = 1.0, pch = 21, print.labs = FALSE, outfile=NULL,...)}
\arguments{
  \item{sp.plot}{a community dataset in \pkg{spacodiR} format (see \code{\link{as.spacodi}})}
  \item{phy}{a phylogenetic tree of class \code{phylo}; see \code{\link[ape]{read.tree}}}
  \item{edge.width}{a value specifying the thickness of branches on the plotted phylogen(ies)}
  \item{lab.adj}{a vector of two values, adjusting plotted group names, with respect to the x- and y-axes}
  \item{tips.adj}{a vector of two values, adjusting the position of plotted symbols (referencing the x- and y-axes)}
  \item{tips.cex}{a value for character expansion of plotted symbols for present species}
  \item{pch}{a value for the plotting character for species-presence symbols: see \code{\link[graphics]{par}}}
  \item{print.labs}{Boolean; returns array of tip labels, group membership of species, and tree to the console}
  \item{outfile}{an optional .pdf file to which to write output}
  \item{\dots}{additional plotting options to be specified}
}
\details{For several groups (fewer than sixteen, for optimal visibility), this function will generate a plot of phylogenies, group labels for each 
tree, and symbols demarcating a particular species that is present within the group. This may be useful in visual representations of phylogenetic overdispersion 
or underdispersion within plots (or other groups). Particularly the \code{tips.adj} option may require fiddling in order for symbols to appear \emph{precisely} where they should.}

\value{A collection of trees (or a single tree), where species present within a predefined group are represented with a dark circle at the tip (if \code{pch=21}).}

\author{Jonathan Eastman}

\examples{
data(sp.example)
attach(sp.example)

phy.dotplot(sp.plot=spl[,1:9], phy=phy, tips.adj=c(0.50,0.55), lab.adj=c(0,1))
}