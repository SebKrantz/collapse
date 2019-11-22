\name{Data Transformation}
\alias{collapse-data-transformations}
\alias{datatrans}

\title{\code{collapse} Data Transformations}

\description{
\code{collapse} provides an ensemble of functions to perform common data transformations efficiently and user friendly:
\itemize{
\item \code{dapply} applies functions to \bold{rows or columns} of matrices and data.frame's.
\item \code{BY} is an S3 generic for \bold{Split-Apply-Combine computing} and can perform aggregation as well as grouped transformations. (for aggregation please also see \code{collap} and \code{Fast Statistical Functions}).
\item \code{TRA} is an S3 generic to efficiently perform (groupwise) \bold{replacement and sweeping out of statistics}. %The basic syntax is \code{TRA(x, xag, g)} where \code{x} is data to be transformed, \code{xag} is some set of aggregate statistics to tranform \code{x} and \code{g} is an optional grouping vector for grouped transformations.
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
All of \code{collapse}'s \code{\link{statistical functions}} have a built-in \code{TRA} argument for faster access (i.e. you can compute (groupwise) statistics and use them to transform your data with a single function call). See also \code{\link{fast-statistical-functions}}
\item \code{fscale / STD} is an S3 generic to perform (groupwise and / or weighted) \bold{scaling / standardizing} of data and is orders of magnitude faster than \code{base::scale}.
\item \code{W}, shorthand for 'within transformation', is an S3 generic to efficiently perform (groupwise and / or weighted) \bold{demeaning / centering} of data. Similarly \code{B}, shorthand for 'between transformation', computes (groupwise and / or weighted) \bold{averages} of data.
\item \code{HDW}, shorthand for 'higher-dimensional within transformation', is an S3 generic to efficiently \bold{center data on multiple groups} and partial-out linear models (possibly involving many levels of fixed effects and interactions). In other words, \code{HDW} efficiently computes \bold{residuals} from (potentially complex) linear models. Similarly \code{HDB}, shorthand for 'higher-dimensional between transformation', computes the corresponding means or \bold{fitted values}.
\item \code{flag / L / F}, \code{fdiff / D} and \code{fgrowth / G} are S3 generics to compute sequences of \bold{lags / leads} and suitably lagged and iterated \bold{differences} and \bold{growth rates} on time-series and panel data. See also \code{\link{time-series-computations}}.
}
}
\section{Table of Functions}{
  \tabular{lll}{\emph{ Function / S3 Generic } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Methods } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
  dapply \tab No methods, works for matrices and data frames \tab apply functions to rows or columns \cr
                 BY \tab default, matrix, data.frame, grouped_df \tab Split-Apply-Combine computing \cr
                 TRA \tab default, matrix, data.frame, grouped_df \tab replace and sweep out statistics \cr
                 fscale / STD \tab default, matrix, data.frame, pseries, pdata.frame, grouped_df \tab scale / standardize data \cr
                 W \tab default, matrix, data.frame, pseries, pdata.frame, grouped_df  \tab demean / center data \cr
                 B \tab default, matrix, data.frame, pseries, pdata.frame, grouped_df  \tab compute means / average data \cr
                 HDW \tab default, matrix, data.frame, pseries, pdata.frame \tab High-dimensional centering and lm residuals \cr
                 HDB \tab default, matrix, data.frame, pseries, pdata.frame \tab High-dimensional averages and lm fitted values \cr
                 flag / L, flead / F \tab default, matrix, data.frame, pseries, pdata.frame, grouped_df  \tab (sequences of) lags / leads \cr
                 fdiff / D \tab default, matrix, data.frame, pseries, pdata.frame, grouped_df  \tab (sequences of lagged/leaded and iterated) differences \cr
                 fgrowth / G \tab default, matrix, data.frame, pseries, pdata.frame, grouped_df  \tab (sequences of lagged/leaded and iterated) growth rates or log-differences
}
}



