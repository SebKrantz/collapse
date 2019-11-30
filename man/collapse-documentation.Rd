\name{collapse-documentation}
\alias{collapse-documentation}

\title{Collapse Documentation & Overview}
\description{
The following table fully summarizes the content of \code{collapse}. The documentation follows a hierarchical structure: This is the top page, linking to topical pages, linking to functions (unless fully documented on the topic page, i.e. only functions with separate pages are linked here). Typing \code{?FUN} always brings up the most relevant documentation page. % Functions with separate documentation entries (apart from the topic page) are linked. % Each topic further has it's own overview page in the documentation.
}

\section{Topics and Functions}{

\tabular{lll}{ \emph{ Topic } \Sexpr{"\u200B"} \Sexpr{"\u200B"} \tab \emph{ Main Features / Keywords} \tab \emph{ Functions }  \cr

  \link[=fast-statistical-functions]{Fast Statistical Functions} \tab Fast (grouped and weighted) statistical functions for vector, matrix, data.frame and grouped_df (dplyr compatible).  \tab \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fNobs}}, \code{\link{fNdistinct}} \cr

\link[=quick-grouping]{Quick Grouping} \tab Fast (ordered or unordered) and reusable groupings from vectors, data.frames, lists... \tab \code{GRP} \cr

\link[=select-replace-vars]{Quick Select and Replace Variables} \tab Quick, flexible and secure select and replace variables from data.frames / data.tables / tibbles / grouped tibbles / pdata.frame etc... %, get data, variables names, variable indices
\tab \code{get_vars}, \code{num_vars}, \code{cat_vars}, \code{char_vars}, \code{fact_vars}, \code{logi_vars}, \code{Date_vars} \cr

\link[=quick-conversion]{Quick Data Converion} \tab Fast conversions, vector > factor / matrix / data.frame / data.table, matrix <> list / data.frame / data.table, array > matrix / data.frame / data.table, list > data.frame / data.table, data.frame <> data.table. \tab \code{qDF}, \code{qDT}, \code{qM}, \code{qF}, \code{qG}, \code{mrtl}, \code{mctl} \cr \cr \cr

\link[=advanced-data-aggregation]{Advanced Data Aggregation} \tab Fast and easy aggregation of multi-type data, using multiple-functions, and fully customizeable. Can perform weighted and parallelized aggregations. \tab \code{collap} \cr \cr \cr

\link[=data-transformations]{Data Transformations} \tab Efficient row- column- data-apply, Split-Apply-Combine computing, fast (grouped and weighted) replace, sweep, scale / standardize, (higher-dimensional) within- and between-transformations (e.g. centering and averaging), complex lin. prediction and partialling out \tab \code{\link{dapply}}, \code{\link{BY}}, \code{\link{TRA}}, \code{\link[=fscale]{fscale/STD}}, \code{\link[=fbetween]{fbetween/B}}, \code{\link[=fwithin]{fwithin/W}}, \code{\link[=HDB]{fHDbetween/HDB}}, \code{\link[=HDW]{fHDwithin/HDW}}. \cr \cr \cr

\link[=time-series-panel-series]{Time-Series and Panel-Series} \tab Efficient (sequences of), lags / leads and (iterated) differences and growth rates / log-differences on (unordered) time-series and panel-data. Panel-data > (ts-)array conversions, panel- auto, partial- and cross - correlation functions \tab \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}}, \code{\link[=fgrowth]{fgrowth/G}}, \code{\link{psmat}}, \code{\link{psacf}}, \code{\link{pspacf}}, \code{\link{psccf}}  \cr \cr \cr

\link[=list-processing]{List Processing} \tab (Recursive) row-binding / unlisting in 2-dimensions, data-apply, search and extract list-elements / list-subsetting, list-identification.  \tab \code{\link{unlist2d}}, \code{\link{rapply2d}}, \code{\link{get_elem}}, \code{\link{has_elem}}, \code{\link{atomic_elem}}, \code{\link{list_elem}}, \code{\link{reg_elem}}, \code{\link{irreg_elem}}, \code{\link{is.regular}}, \code{\link{is.unlistable}}, \code{\link{ldepth}} \cr

\link[=quick-summary]{Quick Summary Statistics} \tab Extremely fast (one-pass, grouped and weighted), summary statistics for cross-sectional and panel-data (within- and between-groups) \tab \code{qsu} \cr \cr \cr \cr

%Visual Data Exploration \tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab histplot, densplot, scplot, corplot \cr

\link[=recode-replace]{Recode and Replace Values} \tab Recode multiple values (exact or regex matching) and replace \code{NaN/Inf/-Inf} and outliers (according to 1- or 2-sided threshold or column standard-deviation) in vectors, matrices or data.frames. \tab \code{Recode}, \code{replace_non_finite}, \code{replace_outliers} \cr \cr \cr

\link[=small-helpers]{(Small) Helper Functions} \tab Set and extract variable labels, remove \code{NA} from vector, check multiple objects or all list-elements identical, pairwise cor cov, seq_along row col, return object with dimnames, row- or col- names set. \tab \code{vlabels}, \code{namlab}, \code{na.rm}, \code{all.identical}, \code{pwcor}, \code{pwcov}, \code{seq_row}, \code{seq_col}, \code{setDimnames}, \code{setRownames}, \code{setColnames} \cr \cr \cr


Data and Global Macros \tab Groningen Growth and Development Centre 10-Sector Database, World Bank World Development dataset, and some global macros containing links to the topical documentation pages (including this page), all exported objects (excluding exported S3 methods and documentation entries), all generic functions, the 2 datasets, all fast functions, all fast statistical (scalar-valued) functions, and all transformation operators / operator-like functions.  \tab \code{GGDC10S, wlddev, .COLLAPSE_TOPICS, .COLLAPSE_ALL, .COLLAPSE_GENERIC, .COLLAPSE_DATA, .FAST_FUN, .FAST_STAT_FUN, .OPERATOR_FUN} \cr

}
}
\author{
Maintainer: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}
}
\seealso{
  \link{collapse-package}
}


