\name{collapse-documentation}
\alias{collapse-documentation}
\docType{package}
\title{Collapse Documentation & Overview}
\description{
The following table fully summarizes the content of \code{\link{collapse}}. The documentation follows a hierarchical structure: This is the main overview page, linking to topical overview pages and associated functions (unless fully documented on the topic page). Typing \code{?FUN} always brings up the most relevant documentation page. % Functions with separate documentation entries (apart from the topic page) are linked. % Each topic further has it's own overview page in the documentation. % , linking to functions % , i.e. only functions with separate pages are linked here
}

\section{Topics and Functions}{

\tabular{lllll}{ \emph{ Topic } \Sexpr{"\u200B"} \Sexpr{"\u200B"} \tab\tab \emph{ Main Features / Keywords} \tab\tab \emph{ Functions }  \cr

  \link[=fast-statistical-functions]{Fast Statistical Functions} \tab\tab Fast (grouped and weighted) statistical functions for vector, matrix, data.frame and grouped_df (dplyr compatible).  \tab\tab \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fNobs}}, \code{\link{fNdistinct}} \cr \cr \cr

\link[=quick-grouping]{Quick Grouping} \tab\tab Fast (ordered or unordered) and reusable groupings from vectors, data.frames, lists... \tab\tab \code{GRP}, \code{as.factor.GRP}, \code{group_names.GRP} \cr \cr \cr

\link[=select-replace-vars]{Quick Select and Replace Variables} \tab\tab Quick, flexible and secure select and replace variables from data.frames / data.tables / tibbles / grouped tibbles / pdata.frame etc... %, get data, variables names, variable indices
\tab\tab \code{get_vars}, \code{num_vars}, \code{cat_vars}, \code{char_vars}, \code{fact_vars}, \code{logi_vars}, \code{Date_vars} \cr

\link[=quick-conversion]{Quick Data Converion} \tab\tab Fast conversions, data.frame <> data.table, matrix <> list / data.frame / data.table, array > matrix / data.frame / data.table, list > data.frame / data.table, vector > factor / matrix / data.frame / data.table. \tab\tab \code{qDF}, \code{qDT}, \code{qM}, \code{qF}, \code{qG}, \code{mrtl}, \code{mctl} \cr \cr \cr

\link[=advanced-data-aggregation]{Advanced Data Aggregation} \tab\tab Fast and easy aggregation of multi-type data, using multiple-functions, and fully customizeable. Can perform weighted and parallelized aggregations. \tab\tab \code{collap} \cr \cr \cr

\link[=data-transformations]{Data Transformations} \tab\tab Efficient row- column- data-apply, Split-Apply-Combine computing, fast (grouped and weighted) replace, sweep, scale / standardize, (higher-dimensional) within- and between-transformations (e.g. centering and averaging), complex lin. prediction and partialling out. + methods for \code{dplyr} and \code{plm}. \tab\tab \code{\link{dapply}}, \code{\link{BY}}, \code{\link{TRA}}, \code{\link[=fscale]{fscale/STD}}, \code{\link[=fbetween]{fbetween/B}}, \code{\link[=fwithin]{fwithin/W}}, \code{\link[=HDB]{fHDbetween/HDB}}, \code{\link[=HDW]{fHDwithin/HDW}}. \cr \cr \cr

\link[=time-series-panel-series]{Time-Series and Panel-Series} \tab\tab Efficient (sequences of), lags / leads and (iterated) differences and growth rates / log-differences on (unordered) time-series and panel-data. Panel-data > (ts-)array conversions, multivariate panel- auto, partial- and cross - correlation functions. + methods for \code{dplyr} and \code{plm}. \tab\tab \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}}, \code{\link[=fgrowth]{fgrowth/G}}, \code{\link{psmat}}, \code{\link{psacf}}, \code{\link{pspacf}}, \code{\link{psccf}}  \cr

\link[=list-processing]{List Processing} \tab\tab (Recursive) list-identification, search and extract list-elements / list-subsetting, data-apply, and row-binding / unlisting in 2-dimensions.  \tab\tab \code{\link{is.regular}}, \code{\link{is.unlistable}}, \code{\link{ldepth}}, \code{\link{has_elem}}, \code{\link{get_elem}}, \code{\link{atomic_elem}}, \code{\link{list_elem}}, \code{\link{reg_elem}}, \code{\link{irreg_elem}}, \code{\link{rapply2d}}, \code{\link{unlist2d}} \cr

\link[=quick-summary]{Quick Summary Statistics} \tab\tab Extremely fast (one-pass, grouped and weighted), summary statistics for cross-sectional and panel-data (within- and between-groups). + methods for \code{plm}. \tab\tab \code{qsu} \cr \cr \cr \cr

%Visual Data Exploration \tab\tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab\tab histplot, densplot, scplot, corplot \cr

\link[=recode-replace]{Recode and Replace Values} \tab\tab Recode multiple values (exact or regex matching) and replace \code{NaN/Inf/-Inf} and outliers (according to 1- or 2-sided threshold or column standard-deviation) in vectors, matrices or data.frames. \tab\tab \code{Recode}, \code{replace_non_finite}, \code{replace_outliers} \cr \cr \cr

\link[=small-helpers]{(Small) Helper Functions} \tab\tab Set and extract variable labels, display variable names, labels and classes, add prefix or postfix to column names, not-in operator, faster nlevels for factors, fast unique vector elements, remove \code{NA}'s from vector, insert \code{NA}'s at random into matrix-like objects, check exact equality of multiple objects or of all elements in a list, pairwise correlations and covariances, seq_along rows or columns or matrix-like objects, return object with dimnames, row- or colnames set, identify categorical and date/time objects, convert factors (or all factors in a list of columns) to numeric. \tab\tab \code{vlabels}, \code{namlab}, \code{add_stub}, \code{\%!in\%}, \code{fnlevels}, \code{funique}, \code{na_rm}, \code{na_insert}, \code{all.identical}, \code{pwcor}, \code{pwcov}, \code{seq_row}, \code{seq_col}, \code{setDimnames}, \code{setRownames}, \code{setColnames}, \code{is.categorical}, \code{is.Date}, \code{as.numeric_factor} \cr \cr \cr


Data and Global Macros \tab\tab Groningen Growth and Development Centre 10-Sector Database, World Bank World Development dataset, and some global macros containing links to the topical documentation pages (including this page), all exported objects (excluding exported S3 methods and documentation entries), all generic functions, the 2 datasets, all fast functions, all fast statistical (scalar-valued) functions, and all transformation operators / operator-like functions.  \tab\tab \code{\link{GGDC10S}, \link{wlddev}, .COLLAPSE_TOPICS, .COLLAPSE_ALL, .COLLAPSE_GENERIC, .COLLAPSE_DATA, .FAST_FUN, .FAST_STAT_FUN, .OPERATOR_FUN} \cr

}

\bold{Maintainer}: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}
}

\seealso{
  \link{collapse-package}
}


