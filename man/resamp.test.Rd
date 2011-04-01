\name{resamp.test}
\alias{resamp.test}
\title{statistical comparison of sets of values by randomization}
\description{compares means by bootstrap resampling of differences between empirical distributions}
\usage{
resamp.test(obs = obs, exp = exp, mu = 0, iter = 10000, two.tailed = FALSE, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{a vector of numeric values}
  \item{exp}{a vector of numeric values}
  \item{mu}{the true difference in means}
  \item{iter}{number of randomization comparisons to perform}
  \item{two.tailed}{Boolean; as default, the test is performed under a one-tailed assumption}
  \item{na.rm}{whether to remove NA data}

}
\details{
If a single value is supplied for \code{obs}, this test equates to finding the quantile in \code{exp} in which \code{obs} would be found (under a one-tailed test); 
see \strong{Examples} and also \code{\link[stats]{ecdf}}}
\value{
A list, whose contents are determined by the above argument: 
  \item{unnamed value}{if \code{two.tailed=TRUE}, this is the two-tailed p-value}
  \item{diffs}{the full resampling distribution of differences between \code{obs} and \code{exp}, given \code{mu} }
  \item{greater}{if \code{two.tailed=FALSE}, this is the p-value associated with the righthand tail}
  \item{lesser}{if \code{two.tailed=FALSE}, this is the p-value associated with the lefthand tail}
}
\author{Jonathan Eastman}
\examples{

## A GENERAL EXAMPLE
# a comparison between two distributions
a=rnorm(n=1000, mean=1, sd=0.5)
b=rnorm(n=1000, mean=0, sd=1)
resamp.test(obs=a, exp=b, two.tailed=FALSE)

# a comparison of a single value to a normal distribution
a=3
b=rnorm(n=1000, mean=0, sd=1)
resamp.test(obs=a, exp=b, two.tailed=FALSE)

# compare above result with ecdf(), in which we compute an empirical 
f=ecdf(b)
print(1-f(a))		# analogous to a one-tailed test as above

## A COMMUNITY PHYLOGENETIC EXAMPLE
## an example of comparing observed and expected Bst at nodes ##
data(sp.example)
attach(sp.example)
s=spl
p=phy

# get nodes for a time-slice
n.t=phy.nodetimes(phy=p,time.range=c(0.25, 0.75))

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