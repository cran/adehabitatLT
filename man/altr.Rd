\name{altr}
\alias{altr}
\alias{dltr}
\title{Attach or Detach Bursts in Objects of Class 'ltraj' to Search
  Path}
\description{
  The burst in the object of class \code{ltraj} is attached to the R
  search path.  This means that the parameters describing the movement
  are searched by R when evaluating a variable, so 
  parameters can be accessed by simply giving their
  names.

}
\usage{
altr(object, bu = burst(object)[1])
dltr(object, bu = burst(object)[1])
}
\arguments{
  \item{object}{object of class \code{ltraj}}
  \item{bu}{a character string indicating the name of the burst to be
    attached or detached}
}
\value{
  The 'environment' is returned invisibly with a '"name"' attribute.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{attach}}}
\examples{
\dontrun{
data(ibex)

altr(ibex)
plot(dt, dist)
dltr(ibex)

altr(ibex, bu = "A160")
hist(rel.angle)
dltr(ibex, bu = "A160")

## same as
altr(ibex, 2)
hist(rel.angle)
dltr(ibex, 2)
}
}
\keyword{data}
