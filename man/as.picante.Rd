\name{as.picante}

\alias{as.picante}

\title{converting between data formats for community phylogenetics}
\description{\code{as.picante} converts from \code{spacodi} or \code{phylocom} data formats into \code{picante} format.}


\usage{as.picante(data, outfile = NULL)}

\arguments{
  \item{data}{a community dataset in \code{phylocom} or \code{spacodi} format}
  \item{outfile}{an optional text file to which to write output}
}

\details{This utility converts a community dataset (either from \code{phylocom} or \pkg{spacodiR}) 
into a format interpretable by \pkg{picante} (see \code{\link[picante]{picante-package}}). \code{phylocom} format is also referred 
to as \code{triplet}-formatting, where plots are within the first column, abundances in the second, and species names in the
third column of the dataframe. The user has the option to save an output file, defined by \code{outfile}.  \code{SPACoDi} format is 
similar to that for \pkg{picante}, where dataframes between these packages are transposed.  \code{SPACoDi} format should have species as row names.
}

\value{An array, formatted for use in \code{picante}}

\author{Jonathan Eastman}

\seealso{\code{\link{as.spacodi}} and \code{\link{as.phylocom}} for converting between \code{phylocom} 
and \code{SPACoDi} formats; see \code{\link[picante]{picante-package}} 
for an R-port of \code{phylocom}}

\examples{
# call example data from SPACoDi
data(sp.example)
attach(sp.example)
spl->d.spacodi  
d.spacodi ## SPACoDi format

# convert to phylocom
as.phylocom(data=spl, picante=FALSE)->d.phylocom
d.phylocom ## phylocom format

# convert dataset to picante
as.picante(data=d.phylocom)->d.picante
d.picante ## picante format

# convert back to SPACoDi 
as.spacodi(data=d.picante)

# run standardized effect size mean nearest taxon distances in picante
ses.mntd(as.picante(spl), cophenetic(phy))

}