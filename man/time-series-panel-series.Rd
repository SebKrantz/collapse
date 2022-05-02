\name{time-series-panel-series} % \name{Time Series and Panel Computations}
\alias{A7-time-series-panel-series}
\alias{time-series-panel-series}
% \alias{tscomp}
\title{Time Series and Panel Series} % \emph{collapse}

\description{
\emph{collapse} provides a flexible and powerful set of functions and classes to work with time-dependent data:
\itemize{
\item \code{\link[=findex_by]{findex_by/iby}} creates an 'indexed_frame': a flexible structure that can be imposed upon any data-frame like object and facilitates \bold{indexed (time-aware) computations on time series and panel data}. Indexed frames are composed of 'indexed_series', which can also be created from vector and matrix-based objects using the \code{reindex} function. Further functions \code{findex/ix}, \code{unindex}, \code{is_irregular} and \code{to_plm} help operate these classes, check for irregularity, and ensure \emph{plm} compatibility. Methods are defined for various time series, data transformation and data manipulation functions in \emph{collapse}.

\item \code{\link{timeid}} efficiently converts numeric time sequences, such as 'Date' or 'POSIXct' vectors, to a \bold{time-factor / integer id}, where a unit-step represents the greatest common divisor of the underlying sequence.

\item \code{\link{flag}}, and the lag- and lead- operators \code{\link{L}} and \code{\link{F}} are S3 generics to efficiently compute sequences of \bold{lags and leads} on regular or irregular / unbalanced time series and panel data.
\item Similarly, \code{\link{fdiff}}, \code{\link{fgrowth}}, and the operators \code{\link{D}}, \code{\link{Dlog}} and \code{\link{G}} are S3 generics to efficiently compute sequences of suitably lagged / leaded and iterated \bold{differences, log-differences and growth rates}. \code{\link[=fdiff]{fdiff/D/Dlog}} can also compute \bold{quasi-differences} of the form \eqn{x_t - \rho x_{t-1}}.
\item \code{\link{fcumsum}} is an S3 generic to efficiently compute \bold{cumulative sums} on time series and panel data. In contrast to \code{\link{cumsum}}, it can handle missing values and supports both grouped and indexed / ordered computations.
\item \code{\link{psmat}} is an S3 generic to efficiently convert panel-vectors / 'indexed_series' and data frames / 'indexed_frame's to \bold{panel series matrices and 3D arrays}, respectively (where time, individuals and variables receive different dimensions, allowing for fast indexation, visualization, and computations).
\item \code{\link{psacf}}, \code{\link{pspacf}} and \code{\link{psccf}} are S3 generics to compute estimates of the \bold{auto-, partial auto- and cross- correlation or covariance functions} for panel-vectors / 'indexed_series', and multivariate versions for data frames / 'indexed_frame's.
}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ S3 Generic }   \tab\tab \emph{ Methods }   \tab\tab \emph{ Description }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"}
                 \code{\link[=findex_by]{findex_by/iby}}, \code{findex/ix}, \code{reindex}, \code{unindex}, \code{is_irregular}, \code{to_plm} \tab\tab For vectors, matrices and data frames / lists.   \tab\tab Fast and flexible time series and panel data classes 'indexed_series' and 'indexed_frame'. \cr
                 \code{\link{timeid}} \tab\tab For time sequences represented by integer or double vectors / objects.  \tab\tab Generate integer time-id/factor \cr

                 \code{\link[=flag]{flag/L/F}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab Compute (sequences of) lags and leads \cr
                 \code{\link[=fdiff]{fdiff/D/Dlog}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab Compute (sequences of lagged / leaded and iterated) (quasi-)differences or log-differences \cr
                 \code{\link[=fgrowth]{fgrowth/G}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab Compute (sequences of lagged / leaded and iterated) growth rates (exact, via log-differencing, or compounded) \cr
                 \code{\link{fcumsum}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab Compute cumulative sums \cr
                 \code{\link{psmat}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab Convert panel data to matrix / array \cr
                 \code{\link{psacf}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab Compute ACF on panel data \cr
                 \code{\link{pspacf}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab Compute PACF on panel data \cr
                 \code{\link{psccf}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab Compute CCF on panel data
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=data-transformations]{Data Transformations}
}
\keyword{ts}
\keyword{manip}
\keyword{documentation}


