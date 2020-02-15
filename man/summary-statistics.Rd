\name{A9-summary-statistics} % \name{Time-Series and Panel Computations}
\alias{A9-summary-statistics}
% \alias{tscomp}
\title{\code{collapse} Summary Statistics} %

\description{
\code{collapse} provides the following functions to efficiently summarize data:
\itemize{
\item \code{\link{qsu}}, shorthand for quick-summary, is an extremely fast summary command inspired by the (xt)summarize command in the STATA statistical software. It computes a set of 7 statistics (nobs, mean, sd, min, max, skewness and kurtosis) using a numerically stable one-pass method. Statistics can be computed weighted, by groups, and also within-and between entities (for multilevel / panel-data).

\item \code{\link{descr}} computes a concise and detailed description of a data.frame, including frequency tables for categorical variables and various statistics and quantiles for numeric variables. It is inspired by \code{Hmisc::describe}, but about 10x faster.

\item \code{\link{pwcor}}, \code{\link{pwcov}} and \code{\link{pwNobs}} compute pairwise correlations, covariances and observation counts, respectively. Pairwise correlations and covariances can be computed together with observation counts and p-values, and output as 3D array (default) or list of matrices. A major feature of \code{pwcor} and \code{pwcov} is the print method displaying all of these statistics in a single correlation table.
}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ Function / S3 Generic }   \tab\tab \emph{ Methods }   \tab\tab \emph{ Description }  \cr
                 \code{\link{qsu}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame}  \tab\tab Fast (grouped, weighted, panel-decomposed) summary statistics \cr
                 \code{\link{descr}} \tab\tab No methods, for data.frame's or lists of vectors  \tab\tab detailed statistical description of data.frame \cr
                 \code{\link{pwcor}} \tab\tab No methods, for matrices or data.frame's \tab\tab pairwise correlations \cr
                 \code{\link{pwcov}} \tab\tab No methods, for matrices or data.frame's \tab\tab pairwise covariances \cr
                 \code{\link{pwNobs}} \tab\tab No methods, for matrices or data.frame's \tab\tab pairwise observation counts \cr
}
}
\seealso{
\link[=A1-fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\keyword{ts}
\keyword{manip}
\keyword{documentation}


