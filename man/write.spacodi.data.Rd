\name{write.spacodi.data}

\alias{write.spacodi.data}

\title{saving community phylogenetics datasets to a file}
\description{\code{write.spacodi.data} generates community phylogenetic datasets to be used in the external program \code{SPACoDi}}
\usage{write.spacodi.data(sp.plot, outfile)}
\arguments{
  \item{sp.plot}{a community dataset formatted for the R-package \pkg{spacodiR} (see \code{\link{as.spacodi}})}
  \item{outfile}{a formatted file for the Windows executable \code{SPACoDi}}
}

\details{This utility writes a species-by-plots matrix into a format readable 
by the external program \code{SPACoDi}, a Windows executable (available at 
  \url{http://ebe.ulb.ac.be/ebe/Software.html}).  
}

\author{Jonathan Eastman}

\seealso{\code{\link{as.spacodi}} for converting from \code{phylocom} to \pkg{spacodiR} formats;
see \code{\link{spacodi.calc}} for community diversity measures of empirical data
}

\examples{
# generate a community-phylogenetics dataset
data(sp.example)
attach(sp.example)

# save the dataset to working directory
write.spacodi.data(sp.plot=spl, outfile="spacodi.formatted.txt")

}