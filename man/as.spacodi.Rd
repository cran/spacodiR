\name{as.spacodi}

\alias{as.spacodi}

\title{converting between data formats for community phylogenetics}
\description{\code{as.spacodi} converts from \code{picante} or \code{phylocom} data formats into \code{spacodi} format.}
\usage{as.spacodi(data, outfile = NULL)}

\arguments{
  \item{data}{a community dataset in \code{phylocom} or \code{picante} format}
  \item{outfile}{an optional text file to which to write output}
}

\details{
This utility converts a community dataset (either from \pkg{phylocom} or \pkg{picante} (see \code{\link[picante]{picante}}) 
into a format interpretable by either this R-package \pkg{spacodiR} or the external program \code{SPACoDi}, 
a Windows executable (available at \url{http://ebe.ulb.ac.be/ebe/Software.html}). Note also that the 
community-dataset format used here is also that called for by the package \pkg{vegan}; see \code{\link[vegan]{vegandocs}}. 

\code{phylocom} format is also referred to as \code{triplet}-formatting, where plots are within the first column, 
abundances in the second, and species names in the third column of the dataframe. \code{picante} format is simply 
the transpose of \code{spacodiR}-formatting of the community dataset: in \code{spacodiR}, species are expected as 
the row names of the dataframe, where plots are represented as the column names. The user has the option to save 
an output file, defined by \code{outfile}. 
}

\value{An array, formatted for use in \code{spacodiR}}

\author{Jonathan Eastman}

\seealso{\code{\link{as.phylocom}}; \code{\link{as.picante}}}

\examples{
# call example data from spacodiR
data(sp.example)
attach(sp.example)
spl->d.spacodi  
d.spacodi ## SPACoDi format

# convert to phylocom
as.phylocom(data=d.spacodi, picante=FALSE)->d.phylocom
d.phylocom ## phylocom format

# convert dataset to picante
as.picante(data=d.phylocom)->d.picante
d.picante ## picante format

# convert back to SPACoDi 
as.spacodi(data=d.picante)

}