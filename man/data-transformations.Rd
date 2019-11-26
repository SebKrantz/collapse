\name{data-transformations}
\alias{data-transformations}

\title{Data Transformations} % \code{collapse} Data Transformations

\description{
\code{collapse} provides an ensemble of functions to perform common data transformations efficiently and user friendly:
\itemize{
\item \code{\link{dapply}} applies functions to \bold{rows or columns} of matrices and data.frame's.
\item \code{\link{BY}} is an S3 generic for \bold{Split-Apply-Combine computing} and can perform aggregation as well as grouped transformations. (for aggregation please also see \code{\link{collap}} and \link[=fast-statistical-functions]{Fast Statistical Functions}).

\item \code{\link{TRA}} is an S3 generic to efficiently perform (groupwise) \bold{replacement and sweeping out of statistics}. %The basic syntax is \code{TRA(x, xag, g)} where \code{x} is data to be transformed, \code{xag} is some set of aggregate statistics to tranform \code{x} and \code{g} is an optional grouping vector for grouped transformations.
Supported operations are:
  \tabular{lll}{\emph{ Integer-id } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ String-id } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
                 1 \tab "replace_fill" \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab replace and overwrite missing values \cr
                 2 \tab "replace" \tab replace but preserve missing values \cr
                 3 \tab "-"   \tab subtract \cr
                 4 \tab "-+"  \tab subtract group-statistics but add average of group statistics \cr
                 5 \tab "/"   \tab divide \cr
                 6 \tab "\%"     \tab compute percentages \cr
                 7 \tab "+" \tab add \cr
                 8 \tab "*" \tab multiply
}
All of \code{collapse}'s \link[=fast-statistical-functions]{Fast Statistical Functions} have a built-in \code{TRA} argument for faster access (i.e. you can compute (groupwise) statistics and use them to transform your data with a single function call).

\item \code{\link[=fscale]{fscale/STD}} is an S3 generic to perform (groupwise and / or weighted) \bold{scaling / standardizing} of data and is orders of magnitude faster than \code{base::scale}.

\item \code{\link{W}}, shorthand for 'within transform', is an S3 generic to efficiently perform (groupwise and / or weighted) \bold{demeaning / centering} of data. Similarly \code{\link{B}}, shorthand for 'between transform', computes (groupwise and / or weighted) \bold{averages} of data.

\item \code{\link{HDW}}, shorthand for 'higher-dimensional within transform', is an S3 generic to efficiently \bold{center data on multiple groups} and partial-out linear models (possibly involving many levels of fixed effects and interactions). In other words, \code{\link{HDW}} efficiently computes \bold{residuals} from (potentially complex) linear models. Similarly \code{\link{HDB}}, shorthand for 'higher-dimensional between transformation', computes the corresponding means or \bold{fitted values}.

\item \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}} and \code{\link[=fgrowth]{fgrowth/G}} are S3 generics to compute sequences of \bold{lags / leads} and suitably lagged and iterated \bold{differences} and \bold{growth rates} on time-series and panel data. More in \link[=time-series-panel-series]{Time-Series and Panel-Series}.
}
}
\section{Table of Functions}{
  \tabular{lll}{\emph{ Function / S3 Generic } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Methods } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
  \code{\link{dapply}} \tab No methods, works with matrices and data frames \tab apply functions to rows or columns \cr
  \code{\link{BY}} \tab \code{default, matrix, data.frame, grouped_df} \tab Split-Apply-Combine computing \cr
  \code{\link{TRA}} \tab \code{default, matrix, data.frame, grouped_df} \tab replace and sweep out statistics \cr
  \code{\link[=fscale]{fscale/STD}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df} \tab scale / standardize data \cr
  \code{\link{W}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab demean / center data \cr
  \code{\link{B}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab compute means / average data \cr
  \code{\link{HDW}} \tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab high-dimensional centering and lm residuals \cr
  \code{\link{HDB}} \tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab high-dimensional averages and lm fitted values \cr
  \code{\link[=flag]{flag/L/F}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab (sequences of) lags / leads \cr
  \code{\link[=fdiff]{fdiff/D}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab (sequences of lagged/leaded and iterated) differences \cr
  \code{\link[=fdiff]{fgrowth/G}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab (sequences of lagged/leaded and iterated) growth rates or log-differences
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \code{\link{collap}}, \link[=time-series-panel-series]{Time-Series and Panel-Series}
}



