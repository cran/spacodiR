\name{spacodi.calc}
\alias{spacodi.calc}


\title{measuring spatial and phylogenetic structuring of diversity in communities }

\description{Considering species-, phylogenetic-, or trait-diversities, this function measures diversity structuring of community samples.}

\usage{spacodi.calc(sp.plot, phy = NULL, sp.traits = NULL, all.together = TRUE, prune = TRUE, pairwise = FALSE, ...)}

\arguments{
  \item{sp.plot}{a community dataset in \pkg{spacodiR} format (see \code{\link{as.spacodi}})}
  \item{phy}{a phylogenetic tree of class \code{phylo} or evolutionary distance matrix between species (see \code{\link[ape]{cophenetic.phylo}}) }
  \item{sp.traits}{a species-by-trait(s) dataframe or a species traits distance matrix (see \code{\link[stats]{dist}})}
  \item{all.together}{Boolean; whether to treat all traits together or separately}
  \item{prune}{Boolean; whether to dynamically prune datasets if mismatches occur} 
  \item{pairwise}{Boolean; whether to return pairwise diversity measures amongst all plots} 
  \item{\dots}{additional arguments to be passed to \code{\link{match.spacodi.data}}}
}

\details{
  \code{spacodi.calc} \strong{requires} a community dataset (species-by-plots matrix; \code{sp.plot}) of numerical abundance, relative abundance, or presence | absence for plots. 
  \code{spacodi.calc} returns statistics of diversity partitioning of plots, considering species diversity and, if additional information is 
  provided, \strong{either} trait or phylogenetic diversities among plots.  If \code{phy=NULL} and \code{sp.traits=NULL}, a measure of partitioning for species 
  diversity will be returned.
  
  In treating each pair of plots as a community unto its own, \code{pairwise=TRUE} will return estimates for diversity structuring for all pairwise combinations of plots. 
  
  If a phylogeny or trait dataset is supplied with species that are not present in the community dataset (i.e., \code{sp.plot}) or vice versa, 
  the user has the option to dynamically prune these datasets to match (\code{prune=TRUE}). If \code{prune=FALSE} and dataset mismatches occur, 
  the function will inevitably return \code{NaN} where plots have fewer than two distinct species sampled. 

  For proper display, please view the package manual online (\url{http://cran.r-project.org/web/packages/spacodiR/spacodiR.pdf})
  
  GLOBAL MEASURES
  \itemize{
    \item{} \eqn{N}: number of local communities sampled
	\item{} \eqn{n_k}: number of individuals sampled in local community \eqn{k}
	\item{} \eqn{f_{ik}}: observed relative abundance of species \eqn{i} in the local community \eqn{k} (\eqn{\sum_i{f_{ik}=1}})
	\item{} \eqn{p_{ik}}: presence (1) or absence (0) of species \eqn{i} in the local community \eqn{k}
	\item{} \eqn{\delta_{ij}}: phyletic or functional (trait) distance between species \eqn{i} and \eqn{j}
  } 
  
  INDICES
  \itemize{
	\item{indices using species abundances:}
	Mean phyletic or functional distance between two individuals in the same local community is expressed as follows:
	\deqn{D_k = \frac{n_k}{n_k-1} \sum_i \sum_j{\delta_{ij}f_{ik}f_{jk} }}  
	This measure is sample-size corrected by the term (\eqn{\frac{n_k}{n_k-1}}) and is only applied where data are counts of individuals (not where data represent a relative measure of species abundances). 
	Mean functional or phyletic distance between individuals of different species in the same local community is computed as follows:
	\deqn{D_k* = \sum_i\sum_{j \neq i}{\delta_{ij}f_{ik}f_{jk} } / \sum_i\sum_{j \neq i}{f_{ik}f_{jk} } }
	Rao's quadratic entropy within samples (average inter-individual distance among samples) is computed either with (\eqn{D_S}) or without (\eqn{D_S*}) sample-size correction:
	\deqn{D_S = \frac{1}{N} \sum_{k \leq N} D_k }
	\deqn{D_S^* = \frac{1}{N} \sum_{k \leq N} D_k^* }
	Rao's quadratic entropy among samples (average distance between individuals from different samples) is as follows, including or excluding intraspecific comparisons (respectively):
	\deqn{D_T = \frac{1}{N(N-1)} \sum_k \sum_{l \neq k} \sum_i \sum_j \delta_{ij}f_{ik}f_{jl} }
	\deqn{D_T^* = \frac{1}{N(N-1)} \sum_k \sum_{l \neq k} [ \sum_i \sum_{j \neq i} \delta_{ij}f_{ik}f_{jl} / \sum_i \sum_{j \neq i} f_{ik}f_{jl} ] }
	
	\itemize{
		\item{}\eqn{P_{ST}} integrates both species turnover and phylogenetic turnover, weighting species according to local abundances. 
		\deqn{ P_{ST}=T_{ST}=1-\frac{D_S}{D_T} } for phylogenetic (\eqn{P_{ST}}) or functional (\eqn{T_{ST}}) distances between individuals.
	
		\item{}\eqn{I_{ST}} is a measure of species turnover accounting for local abundances, but with different properties than measures of species turnover based on species sharing. It is the analogue of \eqn{G_{ST}} in population genetics.
		\deqn{ I_{ST}=1-\frac{D_S}{D_T} } in the special case of spatial partitioning of species diversity, where \eqn{\delta_{ij}=0} if \eqn{i=j} and where \eqn{\delta_{ij}=1} if \eqn{i \neq j} 
	
		\item{}\eqn{B_{ST}} and \eqn{U_{ST}} integrate only phylogenetic or functional (phenotypic) turnover, respectively, weighting species according to local abundances.
		\deqn{ B_{ST}=U_{ST}=1-\frac{D_S^*}{D_T^*} } for phylogenetic (\eqn{B_{ST}}) or functional (\eqn{U_{ST}}) distances between individuals, excluding intraspecific comparisons.
	}

	\item{indices using species occurrences:}
	If we let \eqn{p_{ik}=1} where species \eqn{i} occurs in sample \eqn{k} and \eqn{p_{ik}=0} otherwise, we compute the following local measure of phyletic or functional distinctiveness: 
	\deqn{\Delta_k=\sum_i \sum_{j \neq i} {\delta_{ij}p_{ik}p_{jk} } / \sum_i \sum_{j \neq i} {p_{ik}p_{jk} }} 
	\eqn{\Delta_k} is the mean interspecific distance in local community \eqn{k} and mean distance between distinct species within local communities is then \deqn{\Delta_S = \frac{1}{N} \sum_{k \leq N} \delta_k}
	Mean distance between distinct species among local communities is: \deqn{\Delta_T = \frac{1}{N}  \sum_k \sum_{l \neq k} [ \sum_i \sum_{j \neq i} \delta_{ij}p_{ik}p_{jl} / \sum_i \sum_{j \neq i} p_{ik}p_{jl} ] }
	This formula corrects a typographic error in eq. 16 of Hardy and Senterre (2007) where \eqn{j \neq i} were missing.
	
	\eqn{\Pi_{ST}} expresses phylogenetic turnover without accounting for species local abundances (it is a presence-absence version of \eqn{B_{ST}}):
	\deqn{\Pi_{ST} = 1 - \frac{\Delta_S}{\Delta_T}} 
	
	\item{assumptions:}
	\itemize{
		\item{}\eqn{N} is much is much smaller than the total number of local communities constituting the whole community of which the sampled local communities are representative (this justifies that \eqn{\Delta_T} and \eqn{\Delta_T} are estimated excluding pairs of species from the same local community)	
		\item{}\eqn{n_k} is much smaller than the total number of individuals constituting local community \eqn{k} (this justifies the \eqn{\frac{n_k}{n_k-1}} correction factors)
	}
  }
}

\value{ 
  A named list of at least one element (\code{Ist}) is returned. The size of the returned list is wholly dependent upon given arguments.
  
  SPECIES DIVERSITY STRUCTURING
  \itemize{
  
    \item \code{Ist}: a measure of local species identity excess between 
	individuals, expressing species turnover. It is a form of spatial 
	partition of Gini-Simpson diversity (equivalent to \code{Fst} in population 
	genetics). \code{Ist} considers \strong{only} abundances (or presences) in the 
	species-by-plots matrix.
  } 
  
  PHYLOGENETIC DIVERSITY STRUCTURING
  \itemize{
   
    \item \code{Pst}: a measure of local phyletic proximity excess between 
	individuals, expressing species + phylogenetic turnover. It is a form 
	of spatial partition of Rao's quadratic entropy (equivalent to \code{Nst} in 
	population genetics). \code{Tst} is the analogue for trait data, estimating 
	the spatial partitioning of mean trait-divergence between individuals.
	
	\item \code{Bst}: a measure of local phyletic proximity excess between 
	individuals of distinct species, expressing phylogenetic turnover 
	(independently of species turnover). \code{Ust} is the analogue for 
	trait data, estimating the spatial partitioning of mean trait-divergence 
	between individuals that belong to distinct species.
	
	\item \code{PIst}: Bst analogue for presence/absence data, expressing
	phylogenetic turnover (independently of species turnover). \code{TAUst}
	is the analogue for trait data, estimating mean trait-divergence between 
	distinct species.
		
  }
    
  TRAIT DIVERSITY STRUCTURING
  \itemize{
	\item Measures analogous to those under \code{PHYLOGENETIC DIVERSITY STRUCTURING} can be 
	computed from trait data.  For trait data, these analogues are \code{Tst} (see \code{Pst}), \code{Ust} 
	(see \code{Bst}), and \code{TAUst} (see \code{PIst}). Note: elsewhere, \code{Ust} will be referred to 
	as \code{T*st} but here has been renamed to avoid issues of indexing in R. Trait values 
	are not assumed to follow any particular model of evolution; rather, distances between observed 
	species traits are expected to be uniform in distribution. 
	
    \item If \code{all.together=TRUE}, all traits will be used to generate distance a distance 
	matrix for sampled species.  Where \code{all.together=FALSE} is used, output is generated for each 
	trait independently.  
  }
  
  INTERPRETATION
  \itemize{
  	\item \code{spatial clustering}: species within plots are more 
	phylogenetically related on average than species from distinct plots where \code{ Pst > Ist},
	\code{ Bst > 0}, or  \code{ PIst > 0}. Species are functionally more similar locally than those from
	distinct plots where \code{ Tst > Ist}, \code{ Ust > 0}, or  \code{ TAUst > 0}
  
    \item \code{spatial overdispersion}: species within plots are less 
	phylogenetically related on average than species from distinct plots where \code{ Pst < Ist}, 
	\code{ Bst < 0}, or  \code{ PIst < 0}. Species are functionally less similar locally than are species from
	distinct plots where \code{ Tst < Ist}, \code{ Ust < 0}, or  \code{ TAUst < 0} 
  }	
}

\references{
  HARDY OJ and B SENTERRE. 2007. Characterizing the 
  phylogenetic structure of communities by an additive partitioning of 
  phylogenetic diversity. Journal of Ecology 95:493-506.
  
  HARDY OJ. 2008. Testing the spatial phylogenetic 
  structure of local communities: statistical performances of 
  different null models and test statistics on a locally neutral 
  community. Journal of Ecology 96:914-926.
  
  HARDY OJ and L JOST. 2008. Interpreting and estimating measures of community 
  phylogenetic structuring. Journal of Ecology 96:849-852. 
}

\author{Olivier Hardy, Timothy Paine, and Jonathan Eastman}

\seealso{\code{\link{match.spacodi.data}}; \code{\link{as.spacodi}}}

\examples{
# load a species-by-plots matrix, along with a tree
data(sp.example)
attributes(sp.example)
attach(sp.example)
spl
phy

# community diversity statistics of Hardy and Senterre (2007): tree-based
spacodi.calc(sp.plot = spl, phy = phy)

# community diversity statistics: trait-based with pairwise comparisons
spacodi.calc(sp.plot = spl, phy = phy, pairwise=TRUE)

# community diversity for a pair of traits
spacodi.calc(sp.plot = spl, sp.traits = trt, all.together=TRUE)

# community diversity for a pair of traits, each singly
spacodi.calc(sp.plot = spl, sp.traits = trt, all.together=FALSE)

# Ist: using abundance data only				
spacodi.calc(sp.plot = spl)	

# calculations with missing taxa between tree and sp.plot
# excluding the last five species in sp.plot, 
spacodi.calc(sp.plot = spl[1:15,], phy = phy, prune=TRUE)

# as before but with 'manual' pruning of the datasets
match.spacodi.data(sp.plot=spl[1:15,],phy=phy) -> prn.data
spacodi.calc(sp.plot=prn.data$sp.plot, phy=prn.data$sp.tree)
prn.data$sp.plot
prn.data$sp.tree


							
}
