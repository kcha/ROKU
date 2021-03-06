% Generated by roxygen2 (4.0.2): do not edit by hand
\name{ROKU_as_df}
\alias{ROKU_as_df}
\title{Convert ROKU results as data frame}
\usage{
ROKU_as_df(r, meta = NULL)
}
\arguments{
\item{r}{a list returned by \code{\link{ROKU}}}

\item{meta}{a optional data frame of additional metadata (column-wise) that will
be appended in front of the resulting data frame}
}
\value{
a data frame
}
\description{
Take ROKU list and convert to data frame
}
\examples{
rk <- ROKU(psidata)
df <- ROKU_as_df
head(df)
}
\seealso{
\code{\link{ROKU}}
}

