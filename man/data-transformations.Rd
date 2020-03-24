\name{A6-data-transformations}
\alias{A6-data-transformations}
\alias{.OPERATOR_FUN}
\title{\code{collapse} Data Transformations} %

\description{
\code{collapse} provides an ensemble of functions to perform common data transformations efficiently and user friendly:
\itemize{
\item \code{\link{dapply}} applies functions to \bold{rows or columns} of matrices and data.frame's.
\item \code{\link{BY}} is an S3 generic for \bold{Split-Apply-Combine computing} and can perform aggregation as well as grouped transformations. (for aggregation please also see \code{\link{collap}} and \link[=A1-fast-statistical-functions]{Fast Statistical Functions}).

\item \code{\link{TRA}} is an S3 generic to efficiently perform (groupwise) \bold{replacement and sweeping out of statistics}. %The basic syntax is \code{TRA(x, xag, g)} where \code{x} is data to be transformed, \code{xag} is some set of aggregate statistics to tranform \code{x} and \code{g} is an optional grouping vector for grouped transformations.
Supported operations are:
  \tabular{lllll}{\emph{ Integer-id }  \tab\tab \emph{ String-id }  \tab\tab \emph{ Description }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"} % \Sexpr{"\u200B"} \Sexpr{"\u200B"}
                 1 \tab\tab "replace_fill"   \tab\tab replace and overwrite missing values \cr
                 2 \tab\tab "replace" \tab\tab replace but preserve missing values \cr
                 3 \tab\tab "-"   \tab\tab subtract \cr
                 4 \tab\tab "-+"  \tab\tab subtract group-statistics but add group-frequency weighted average of group statistics \cr
                 5 \tab\tab "/"   \tab\tab divide \cr
                 6 \tab\tab "\%"     \tab\tab compute percentages \cr
                 7 \tab\tab "+" \tab\tab add \cr
                 8 \tab\tab "*" \tab\tab multiply \cr
                 9 \tab\tab "\%\%" \tab\tab modulus \cr
                 10 \tab\tab "-\%\%" \tab\tab subtract modulus

}
All of \code{collapse}'s \link[=A1-fast-statistical-functions]{Fast Statistical Functions} have a built-in \code{TRA} argument for faster access (i.e. you can compute (groupwise) statistics and use them to transform your data with a single function call).

\item \code{\link[=fscale]{fscale/STD}} is an S3 generic to perform (groupwise and / or weighted) \bold{scaling / standardizing} of data and is orders of magnitude faster than \code{base::scale}.

\item \code{\link[=fwithin]{fwithin/W}} is an S3 generic to efficiently perform (groupwise and / or weighted) \bold{within-transformations / demeaning / centering} of data. Similarly \code{\link[=fbetween]{fbetween/B}} computes (groupwise and / or weighted) \bold{between-transformations / averages}.

\item \code{\link[=HDW]{fHDwithin/HDW}}, shorthand for 'higher-dimensional within transform', is an S3 generic to efficiently \bold{center data on multiple groups} and partial-out linear models (possibly involving many levels of fixed effects and interactions). In other words, \code{\link[=HDW]{fHDwithin/HDW}} efficiently computes \bold{residuals} from (potentially complex) linear models. Similarly \code{\link[=HDB]{fHDbetween/HDB}}, shorthand for 'higher-dimensional between transformation', computes the corresponding means or \bold{fitted values}.

\item \code{\link{fFtest}} is a fast implementation of the R-Squared based F-test, to test \bold{exclusion restrictions} on linear models potentially involving multiple large factors (fixed effects). It internally utilizes \code{\link{fHDwithin}} to project out factors while counting the degrees of freedom.

\item \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}} and \code{\link[=fgrowth]{fgrowth/G}} are S3 generics to compute sequences of \bold{lags / leads} and suitably lagged and iterated \bold{differences} and \bold{growth rates} on time-series and panel data. More in \link[=A7-time-series-panel-series]{Time-Series and Panel-Series}.

\item \code{STD, W, B, HDW, HDB, L, D} and \code{G} are parsimonious wrappers around the \code{f-} functions above representing the corresponding transformation 'operators'. They have additional capabilities when applied to data-frames (i.e. variable selection, formula input, auto-renaming and id-variable preservation), and are easier to employ in regression formulas, but are otherwise identical in functionality.
}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ Function / S3 Generic }   \tab\tab \emph{ Methods }  \tab\tab \emph{ Description }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"} % \Sexpr{"\u200B"} \Sexpr{"\u200B"}
  \code{\link{dapply}} \tab\tab No methods, works with matrices and data frames \tab\tab apply functions to rows or columns \cr
  \code{\link{BY}} \tab\tab \code{default, matrix, data.frame, grouped_df} \tab\tab Split-Apply-Combine computing \cr
  \code{\link{TRA}} \tab\tab \code{default, matrix, data.frame, grouped_df} \tab\tab replace and sweep out statistics \cr
  \code{\link[=fscale]{fscale/STD}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df} \tab\tab scale / standardize data \cr
  \code{\link[=fwithin]{fwithin/W}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab demean / center data \cr
  \code{\link[=fbetween]{fbetween/B}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab compute means / average data \cr
  \code{\link[=HDW]{fHDwithin/HDW}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab\tab high-dimensional centering and lm residuals \cr
  \code{\link[=HDB]{fHDbetween/HDB}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab\tab high-dimensional averages and lm fitted values \cr
  \code{\link{fFtest}} \tab\tab No methods, it's a standalone test to which data needs to be supplied.  \tab\tab fast F-test of exclusion restrictions on linear models (involving factors) \cr
  \code{\link[=flag]{flag/L/F}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab (sequences of) lags / leads \cr
  \code{\link[=fdiff]{fdiff/D}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab (sequences of lagged/leaded and iterated) differences \cr
  \code{\link[=fdiff]{fgrowth/G}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab (sequences of lagged/leaded and iterated) growth rates or log-differences
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=A1-fast-statistical-functions]{Fast Statistical Functions}, \code{\link{collap}}, \link[=A7-time-series-panel-series]{Time-Series and Panel-Series}
}
\keyword{manip}
\keyword{documentation}


