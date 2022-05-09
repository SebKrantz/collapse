\name{summary-statistics} % \name{Time Series and Panel Computations}
\alias{A9-summary-statistics}
\alias{summary-statistics}
% \alias{tscomp}
\title{Summary Statistics} % \emph{collapse}

\description{
\emph{collapse} provides the following functions to efficiently summarize and examine data:
\itemize{
\item \code{\link{qsu}}, shorthand for quick-summary, is an extremely fast summary command inspired by the (xt)summarize command in the STATA statistical software. It computes a set of 7 statistics (nobs, mean, sd, min, max, skewness and kurtosis) using a numerically stable one-pass method. Statistics can be computed weighted, by groups, and also within-and between entities (for multilevel / panel data).

\item \code{\link{qtab}}, shorthand for quick-table, is a faster and more versatile alternative to \code{\link{table}}. Notably, it also supports tabulations with frequency weights, as well as computing a statistic over combinations of variables. 'qtab's inherit the 'table' class, allowing for seamless application of 'table' methods.

\item \code{\link{descr}} computes a concise and detailed description of a data frame, including (sorted) frequency tables for categorical variables and various statistics and quantiles for numeric variables. It is inspired by \code{Hmisc::describe}, but about 10x faster.

\item \code{\link{pwcor}}, \code{\link{pwcov}} and \code{\link{pwnobs}} compute (weighted) pairwise correlations, covariances and observation counts on matrices and data frames. Pairwise correlations and covariances can be computed together with observation counts and p-values. The elaborate print method displays all of these statistics in a single correlation table.

\item \code{\link{varying}} very efficiently checks for the presence of any variation in data (optionally) within groups (such as panel-identifiers). A variable is variant if it has at least 2 distinct non-missing data points.

% \item \code{\link{fFtest}} is a fast implementation of the R-Squared based F-test, to test \bold{exclusion restrictions} in linear models potentially involving multiple large factors (fixed effects). It internally utilizes \code{\link{fhdwithin}} to project out factors while counting the degrees of freedom.

}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ Function / S3 Generic }   \tab\tab \emph{ Methods }   \tab\tab \emph{ Description }  \cr
                 \code{\link{qsu}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame}  \tab\tab Fast (grouped, weighted, panel-decomposed) summary statistics \cr
                 \code{\link{qtab}} \tab\tab No methods, for data frames or vectors  \tab\tab Fast (weighted) cross tabulation \cr
                 \code{\link{descr}} \tab\tab No methods, for data frames or lists of vectors  \tab\tab Detailed statistical description of data frame \cr
                 \code{\link{pwcor}} \tab\tab No methods, for matrices or data frames \tab\tab Pairwise (weighted) correlations \cr
                 \code{\link{pwcov}} \tab\tab No methods, for matrices or data frames \tab\tab Pairwise (weighted) covariances \cr
                 \code{\link{pwnobs}} \tab\tab No methods, for matrices or data frames \tab\tab Pairwise observation counts \cr
                \code{\link{varying}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df} \tab\tab Fast variation check \cr
%                \code{\link{fFtest}} \tab\tab No methods, it's a standalone test to which data needs to be supplied.  \tab\tab Fast F-test of exclusion restrictions in linear models (with factors variables) \cr
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=fast-statistical-functions]{Fast Statistical Functions}
}
\keyword{manip}
\keyword{documentation}


