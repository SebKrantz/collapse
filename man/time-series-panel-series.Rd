\name{time-series-panel-series} % \name{Time-Series and Panel Computations}
\alias{time-series-panel-series}
% \alias{tscomp}
\title{\code{collapse} Time-Series and Panel-Series} %

\description{
\code{collapse} provides the following functions to work with time-dependent data:
\itemize{
\item \code{\link{flag}}, and the lag- and lead- operators \code{\link{L}} and \code{\link{F}} are S3 generics to efficiently compute sequences of \bold{lags and leads} on ordered or unordered time-series and panel data.
\item \code{\link{fdiff}}, \code{\link{fgrowth}}, and the operators \code{\link{D}} and \code{\link{G}} are S3 generics to efficiently compute sequences of suitably lagged and iterated \bold{differences and growth rates or log-differences} on ordered or unordered tim-series and panel data. They can also be used to compute forward (leaded) differences or growth rates.
\item \code{\link{psmat}} is an S3 generic to efficiently expand panel-vectors or \code{plm::pseries} and data.frame's or \code{plm::pdata.frame}'s to \bold{panel-series matrices and 3D arrays}, respectively.
\item \code{\link{psacf}}, \code{\link{pspacf}} and \code{\link{psccf}} are S3 generics to compute estimates of the \bold{auto-, partial auto- and cross- correlation or covariance functions} for panel-vectors or \code{plm::pseries}, and multivariate versions for data.frame's or \code{plm::pdata.frame}'s.
}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ S3 Generic }   \tab\tab \emph{ Methods }   \tab\tab \emph{ Description }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"}
                 \code{\link[=flag]{flag/L/F}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab compute (sequences of) lags and leads \cr
                 \code{\link[=fdiff]{fdiff/D}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab compute (sequences of lagged and iterated) differences \cr
                 \code{\link[=fgrowth]{fgrowth/G}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab compute (sequences of lagged and iterated) growth rates or log-differences \cr
                 \code{\link{psmat}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab convert panel-data to matrix/array \cr
                 \code{\link{psacf}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab compute ACF on panel-data \cr
                 \code{\link{pspacf}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab compute PACF on panel-data \cr
                 \code{\link{psccf}} \tab\tab \code{default, pseries, data.frame, pdata.frame} \tab\tab compute CCF on panel-data
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=data-transformations]{Data Transformations}
}



