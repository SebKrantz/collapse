\name{collapse-documentation}
\alias{collapse-documentation}

\title{Collapse Documentation & Overview}
\description{
The following table fully summarizes the content of \code{collapse}. The documentation follows a hierarchical structure: This is the top page, linking to topical pages, linking to functions (unless fully documented on the topic page, i.e. only functions with separate pages are linked here). Typing \code{?FUN} always brings up the most relevant documentation page. % Functions with separate documentation entries (apart from the topic page) are linked. % Each topic further has it's own overview page in the documentation.
}

\section{Topics and Functions}{

\tabular{lll}{ \emph{ Topic } \Sexpr{"\u200B"} \Sexpr{"\u200B"} \tab \emph{ Main Features / Keywords} \tab \emph{ Functions }  \cr

  \link[=fast-statistical-functions]{Fast Statistical Functions} \tab Fast (grouped and weighted) statistical functions for vector, matrix, data.frame and grouped_df (dplyr compatible).  \tab \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fNobs}}, \code{\link{fNdistinct}} \cr

\link[=GRP]{Quick Grouping} \tab Fast (ordered or unordered) and reusable groupings from vectors, data.frames, lists... \tab \code{GRP} \cr

\link[=get_vars]{Quick Select and Replace Variables} \tab Quick and secure select and replace variables from data.frames / data.tables / tibbles / grouped tibbles / pdata.frame etc... %, get data, variables names, variable indices
\tab \code{get_vars}, \code{num_vars}, \code{cat_vars}, \code{char_vars}, \code{fact_vars}, \code{logi_vars}, \code{Date_vars} \cr

\link[=qDT]{Quick Data Converion} \tab Fast conversions, vector > factor / matrix / data.frame / data.table, matrix <> list / data.frame / data.table, array > matrix / data.frame / data.table, list > data.frame / data.table, data.frame <> data.table. \tab \code{qDF}, \code{qDT}, \code{qM}, \code{qF}, \code{qG}, \code{mrtl}, \code{mctl} \cr \cr \cr

\link[=collap]{Advanced Data Aggregation} \tab Fast and easy aggregation of multi-type data, using multiple-functions, and fully customizeable. Can perform weighted and parallelized aggregations. \tab \code{collap} \cr \cr \cr

\link[=data-transformations]{Data Transformations} \tab Efficient row- column- data-apply, Split-Apply-Combine computing, fast (grouped and weighted) replace, sweep, scale / standardize, (higher-dimensional) within- and between-transformations (e.g. centering and averaging), complex lin. prediction and partialling out \tab \code{\link{dapply}}, \code{\link{BY}}, \code{\link{TRA}}, \code{\link[=fscale]{fscale/STD}}, \code{\link{W}}, \code{\link{B}}, \code{\link{HDW}}, \code{\link{HDB}}. \cr \cr \cr

\link[=time-series-panel-series]{Time-Series and Panel-Series} \tab Efficient (sequences of), lags / leads and (iterated) differences and growth rates / log-differences on (unordered) time-series and panel-data. Panel-data > (ts-)array conversions, panel- auto, partial- and cross - correlation functions \tab \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}}, \code{\link[=fgrowth]{fgrowth/G}}, \code{\link{psmat}}, \code{\link{psacf}}, \code{\link{pspacf}}, \code{\link{psccf}}  \cr \cr \cr

\link[=list-processing]{List Processing} \tab (Recursive) row-binding / unlisting in 2-dimensions, data-apply, search and extract list-elements / list-subsetting, list-identification.  \tab \code{\link{unlist2d}}, \code{\link{rapply2d}}, \code{\link{get_elem}}, \code{\link{has_elem}}, \code{\link{atomic_elem}}, \code{\link{list_elem}}, \code{\link{reg_elem}}, \code{\link{irreg_elem}}, \code{\link{is.regular}}, \code{\link{is.unlistable}}, \code{\link{ldepth}} \cr

\link[=qsu]{Summary Statistics} \tab Extremely fast (one-pass, grouped and weighted), summary statistics for cross-sectional and panel-data (within- and between-groups) \tab \code{qsu} \cr

%Visual Data Exploration \tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab histplot, densplot, scplot, corplot \cr
\link[=vlabels]{Small (Helper) Functions} \tab set and extract variable labels, return object with dimnames, row- or col- names set, seq row col, pairwise cov cor, vector remove NA's, check all elements identical, recode values, replace NA/NAN/Inf,-Inf. \tab \code{\link{vlabels}}, \code{\link{namlab}}, \code{\link{seq_row}}, \code{\link{seq_col}}, \code{\link{setRownames}}, \code{\link{setColnames}}, \code{\link{setDimnames}}, \code{\link{pwcor}}, \code{\link{pwcov}}, \code{\link{na.rm}}, \code{\link{all.identical}}, Recode, replace_non_finite \cr
}
}
\author{
Maintainer: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}
}
\seealso{
  \link{collapse-package}
}


