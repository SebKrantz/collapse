\name{data-transformations}
\alias{A6-data-transformations}
\alias{data-transformations}
\alias{.OPERATOR_FUN}
\title{Data Transformations} % \emph{collapse}

\description{
\emph{collapse} provides an ensemble of functions to perform common data transformations efficiently and user friendly:
\itemize{
\item \code{\link{dapply}} \bold{applies functions to rows or columns} of matrices and data frames, preserving the data format.
\item \code{\link{BY}} is an S3 generic for efficient \bold{Split-Apply-Combine computing}, similar to \code{\link{dapply}}.

\item A set of arithmetic operators facilitates \bold{row-wise} \code{\link{\%rr\%}}, \code{\link{\%r+\%}}, \code{\link{\%r-\%}}, \code{\link{\%r*\%}}, \code{\link{\%r/\%}} and
\bold{column-wise} \code{\link{\%cr\%}}, \code{\link{\%c+\%}}, \code{\link{\%c-\%}}, \code{\link{\%c*\%}}, \code{\link{\%c/\%}} \bold{replacing and sweeping operations} involving a vector and a matrix or data frame / list. Since v1.7, the operators \code{\link{\%+=\%}}, \code{\link{\%-=\%}}, \code{\link{\%*=\%}} and \code{\link{\%/=\%}} do column- and element- wise math by reference, and the function \code{\link{setop}} can also perform sweeping out rows by reference.

\item \code{\link[=TRA]{(set)TRA}} is a more advanced S3 generic to efficiently perform \bold{(groupwise) replacing and sweeping out of statistics}, either by creating a copy of the data or by reference. %The basic syntax is \code{TRA(x, xag, g)} where \code{x} is data to be transformed, \code{xag} is some set of aggregate statistics to tranform \code{x} and \code{g} is an optional grouping vector for grouped transformations.
Supported operations are:
  \tabular{lllll}{\emph{ Integer-id }  \tab\tab \emph{ String-id }  \tab\tab \emph{ Description }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"} % \Sexpr{"\u200B"} \Sexpr{"\u200B"}
                 0 \tab\tab "replace_NA"   \tab\tab replace missing values \cr
                 1 \tab\tab "replace_fill"   \tab\tab replace everything \cr
                 2 \tab\tab "replace" \tab\tab replace data but preserve missing values \cr
                 3 \tab\tab "-"   \tab\tab subtract \cr
                 4 \tab\tab "-+"  \tab\tab subtract group-statistics but add group-frequency weighted average of group statistics \cr
                 5 \tab\tab "/"   \tab\tab divide \cr
                 6 \tab\tab "\%"     \tab\tab compute percentages \cr
                 7 \tab\tab "+" \tab\tab add \cr
                 8 \tab\tab "*" \tab\tab multiply \cr
                 9 \tab\tab "\%\%" \tab\tab modulus \cr
                 10 \tab\tab "-\%\%" \tab\tab subtract modulus

}
All of \emph{collapse}'s \link[=fast-statistical-functions]{Fast Statistical Functions} have a built-in \code{TRA} argument for faster access (i.e. you can compute (groupwise) statistics and use them to transform your data with a single function call).

\item \code{\link[=fscale]{fscale/STD}} is an S3 generic to perform (groupwise and / or weighted) \bold{scaling / standardizing} of data and is orders of magnitude faster than \code{\link{scale}}.

\item \code{\link[=fwithin]{fwithin/W}} is an S3 generic to efficiently perform (groupwise and / or weighted) \bold{within-transformations / demeaning / centering} of data. Similarly \code{\link[=fbetween]{fbetween/B}} computes (groupwise and / or weighted) \bold{between-transformations / averages} (also a lot faster than \code{\link{ave}}).

\item \code{\link[=HDW]{fhdwithin/HDW}}, shorthand for 'higher-dimensional within transform', is an S3 generic to efficiently \bold{center data on multiple groups and partial-out linear models} (possibly involving many levels of fixed effects and interactions). In other words, \code{\link[=HDW]{fhdwithin/HDW}} efficiently computes \bold{residuals} from linear models. Similarly \code{\link[=HDB]{fhdbetween/HDB}}, shorthand for 'higher-dimensional between transformation', computes the corresponding means or \bold{fitted values}.

%\item \code{flm} is an efficient function for bare-bones (weighted) \bold{linear model fitting}. It supports 6 different fitting methods, 4 from base R, and 2 utilizing the \emph{RcppArmadillo} or \emph{RcppEigen} packages.

\item \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D/Dlog}} and \code{\link[=fgrowth]{fgrowth/G}} are S3 generics to compute sequences of \bold{lags / leads} and suitably lagged and iterated (quasi-, log-) \bold{differences} and \bold{growth rates} on time series and panel data. \code{\link{fcumsum}} flexibly computes (grouped, ordered) cumulative sums. More in \link[=time-series-panel-series]{Time Series and Panel Series}.

\item \code{STD, W, B, HDW, HDB, L, D, Dlog} and \code{G} are parsimonious wrappers around the \code{f-} functions above representing the corresponding transformation 'operators'. They have additional capabilities when applied to data-frames (i.e. variable selection, formula input, auto-renaming and id-variable preservation), and are easier to employ in regression formulas, but are otherwise identical in functionality.
}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ Function / S3 Generic }   \tab\tab \emph{ Methods }  \tab\tab \emph{ Description }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"} % \Sexpr{"\u200B"} \Sexpr{"\u200B"}
  \code{\link{dapply}} \tab\tab No methods, works with matrices and data frames \tab\tab Apply functions to rows or columns \cr
  \code{\link{BY}} \tab\tab \code{default, matrix, data.frame, grouped_df} \tab\tab Split-Apply-Combine computing \cr
  \code{\link[=arithmetic]{\%(r/c)(r/+/-/*//)\%}} \tab\tab No methods, works with matrices and data frames / lists \tab\tab Row- and column-arithmetic \cr
  \code{\link[=TRA]{(set)TRA}} \tab\tab \code{default, matrix, data.frame, grouped_df} \tab\tab Replace and sweep out statistics (by reference) \cr
  \code{\link[=fscale]{fscale/STD}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df} \tab\tab Scale / standardize data \cr
  \code{\link[=fwithin]{fwithin/W}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab Demean / center data \cr
  \code{\link[=fbetween]{fbetween/B}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab Compute means / average data \cr
  \code{\link[=HDW]{fhdwithin/HDW}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab\tab High-dimensional centering and lm residuals \cr
  \code{\link[=HDB]{fhdbetween/HDB}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab\tab High-dimensional averages and lm fitted values \cr
  % \code{\link{flm}} \tab\tab No methods, for matrices \tab\tab Linear model fitting \cr
  \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D/Dlog}}, \code{\link[=fdiff]{fgrowth/G}}, \code{\link{fcumsum}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab\tab (Sequences of) lags / leads, differences, growth rates and cumulative sums
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=time-series-panel-series]{Time Series and Panel Series}
}
\keyword{manip}
\keyword{documentation}


