\name{as.phylocom}

\alias{as.phylocom}

\title{converting between data formats for community phylogenetics}
\description{\code{as.phylocom} converts from \code{spacodi} or \code{picante} data formats into \code{phylocom} format.}
\usage{as.phylocom(data, picante=FALSE, outfile = NULL)}

\arguments{
  \item{data}{a community dataset in either \pkg{spacodiR} or \pkg{picante} format (with species as row names and plots as column names or vice versa)}
  \item{picante}{Boolean; if current data format is \code{spacodi}, argument should be \code{picante=FALSE}; otherwise, format is assumed to be for \code{picante}, with plots as row names}
  \item{outfile}{an optional text file to which to write output}
}

\details{This utility converts a species-by-plots matrix into \code{triplet} format, which is readable by the external program \code{phylocom}. If
\code{picante=TRUE}, the data are expected to be in the form used for \pkg{picante} (i.e., a plots-by-species matrix; \code{\link[picante]{picante-package}}).  If the user selects 
\code{picante=FALSE}, the data are expected to be in the form used for \pkg{spacodiR} (i.e., a species-by-plots matrix).
The user has the option to save an output file, defined by \code{outfile}. 
}

\value{A named array, formatted for use in \code{phylocom}; note that while the R-object returned by this function has column names, if output is written to a file, 
the header is dropped (as appropriate for use in the external \code{phylocom} executable: \url{http://www.phylodiversity.net/phylocom/}).}

\references{WEBB CO, DD ACKERLY and SW KEMBEL. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Bioinformatics 24:2098-2100.}

\author{Jonathan Eastman}

\seealso{\code{\link{as.spacodi}} and \code{\link{as.picante}} for converting between \code{phylocom} 
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

}