\name{A0-collapse-documentation}
\alias{A0-collapse-documentation}
\alias{collapse-documentation}
\alias{.COLLAPSE_TOPICS}
\alias{.COLLAPSE_ALL}
\alias{.COLLAPSE_GENERIC}
\alias{.COLLAPSE_DATA}
\docType{package}
\title{Collapse Documentation & Overview}
\description{
The following table fully summarizes the contents of \code{\link{collapse}}. The documentation follows a hierarchical structure: This is the main overview page, linking to topical overview pages and associated function pages (unless functions are documented on the topic page). Calling \code{?FUN} brings up the documentation page for \code{FUN}, with links to associated topic pages and closely related functions. % Calling \code{help(FUN)} still brings up the right / most relevant page documenting the function. %  % Functions with separate documentation entries (apart from the topic page) are linked. % Each topic further has it's own overview page in the documentation. % , linking to functions % , i.e. only functions with separate pages are linked here
}

\section{Topics and Functions}{

\tabular{lllll}{ \emph{ Topic }  \tab\tab \emph{ Main Features / Keywords} \tab\tab \emph{ Functions }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"}

  \link[=A1-fast-statistical-functions]{Fast Statistical Functions} \tab\tab Fast (grouped and weighted) statistical functions for vector, matrix, data.frame and grouped_df (dplyr compatible).  \tab\tab \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fNobs}}, \code{\link{fNdistinct}} \cr \cr \cr

\link[=A2-fast-grouping]{Fast Grouping} \tab\tab Fast (ordered or unordered) groupings from vectors, data.frames, lists. 'GRP' objects are extremely efficient inputs for \code{collapse}'s fast functions (to optimize different / repeated computations over the same groups).  \tab\tab \code{GRP}, \code{as.factor.GRP}, \code{group_names.GRP}, \code{is.GRP} \cr \cr \cr

\link[=A3-select-replace-vars]{Quick Select and Replace Variables} \tab\tab Quick and flexible select and replace (or add) variables from (to) data.frames / data.tables / tibbles etc... (speed about 2x '[' for selecting and 4x '[<-' for replacing). %, get data, variables names, variable indices
\tab\tab \code{get_vars/gv}, \code{add_vars/av}, \code{num_vars/nv}, \code{cat_vars}, \code{char_vars}, \code{fact_vars}, \code{logi_vars}, \code{Date_vars} \cr

\link[=A4-quick-conversion]{Quick Data Conversion} \tab\tab Quick conversions: data.frame <> data.table | matrix <> list, data.frame, data.table | array > matrix, data.frame, data.table | list > data.frame, data.table | vector > factor, matrix, data.frame, data.table. \tab\tab \code{qDF}, \code{qDT}, \code{qM}, \code{qF}, \code{qG}, \code{mrtl}, \code{mctl} \cr \cr \cr

\link[=A5-advanced-aggregation]{Advanced Data Aggregation} \tab\tab Fast and easy (weighted and parallelized) aggregation of multi-type data, with (multiple) functions applied to numeric and categorical columns. Also supports fully customized aggregation tasks mapping functions directly to columns. \tab\tab \code{collap}, \code{collapv}, \code{collapg} \cr \cr \cr

\link[=A6-data-transformations]{Data Transformations} \tab\tab Efficient row- and column- data-apply and Split-Apply-Combine computing. Fast (grouped and weighted) replacing and sweeping of statistics, scaling / standardizing, (higher-dimensional) within- and between-transformations (i.e. centering and averaging), complex linear prediction and partialling out. A fast F-test for linear models with (large) factors. Additional methods for grouped_df (\code{dplyr}) and pseries, pdata.frame (\code{plm}). \tab\tab \code{\link{dapply}}, \code{\link{BY}}, \code{\link{TRA}}, \code{\link[=fscale]{fscale/STD}}, \code{\link[=fbetween]{fbetween/B}}, \code{\link[=fwithin]{fwithin/W}}, \code{\link[=HDB]{fHDbetween/HDB}}, \code{\link[=HDW]{fHDwithin/HDW}}, \code{\link{fFtest}} \cr \cr \cr

\link[=A7-time-series-panel-series]{Time-Series and Panel-Series} \tab\tab Fast (sequences of) lags / leads and (lagged / leaded and iterated) differences and growth rates / log-differences on (unordered) time-series and panel-data. Panel-data to (ts-)array conversions. Multivariate panel- auto-, partial- and cross- correlation functions. Additional methods for grouped_df (\code{dplyr}) and pseries, pdata.frame (\code{plm}). \tab\tab \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}}, \code{\link[=fgrowth]{fgrowth/G}}, \code{\link{psmat}}, \code{\link{psacf}}, \code{\link{pspacf}}, \code{\link{psccf}}  \cr

\link[=A8-list-processing]{List Processing} \tab\tab (Recursive) list search and identification, search and extract list-elements / list-subsetting, rapply to lists of data.frame's / data objects, and (fast) generalized recursive row-binding / unlisting in 2-dimensions / to data.frame.  \tab\tab \code{\link{is.regular}}, \code{\link{is.unlistable}}, \code{\link{ldepth}}, \code{\link{has_elem}}, \code{\link{get_elem}}, \code{\link{atomic_elem}}, \code{\link{list_elem}}, \code{\link{reg_elem}}, \code{\link{irreg_elem}}, \code{\link{rapply2d}}, \code{\link{unlist2d}} \cr

\link[=A9-summary-statistics]{Summary Statistics} \tab\tab Extremely fast (one-pass, grouped and weighted), summary statistics for cross-sectional and complex multilevel / panel data, with additional methods for pseries and pdata.frame (\code{plm}). Efficient detailed description of data.frame. Pairwise correlations and covariances (with observation count, p-value and pretty printing), pairwise observation count. \tab\tab \code{\link{qsu}}, \code{\link{descr}}, \code{\link{pwcor}}, \code{\link{pwcov}}, \code{\link{pwNobs}} \cr \cr \cr \cr

%  (within- and between-groups) ww

%Visual Data Exploration \tab\tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab\tab histplot, densplot, scplot, corplot \cr

\link[=AA1-recode-replace]{Recode and Replace Values} \tab\tab Recode multiple values (exact or regex matching) and replace \code{NaN/Inf/-Inf} and outliers (according to 1- or 2-sided threshold or column standard-deviation) in vectors, matrices or data.frames. \tab\tab \code{Recode}, \code{replace_non_finite}, \code{replace_outliers} \cr \cr \cr

\link[=AA2-small-helpers]{Small (Helper) Functions} \tab\tab Set and extract variable labels, extract variable classes and C storage types, display variable names, labels and classes, add / remove prefix or postfix to / from column names, not-in operator, matching with error message for non-matched, faster nlevels for factors, fast unique vector elements, faster interactions, remove \code{NA}'s from vector, insert \code{NA}'s at random into matrix-like objects, check exact or near / numeric equality of multiple objects or of all elements in a list, \code{seq_along} rows or columns or matrix-like objects, return object with dimnames, row- or colnames set, identify categorical and date(-time) objects, convert factors (or all factors in a list) to numeric or character by converting levels. \tab\tab \code{vlabels}, \code{vclasses}, \code{vtypes}, \code{namlab}, \code{add_stub}, \code{rm_stub}, \code{\%!in\%}, \code{ckmatch}, \code{fnlevels}, \code{funique}, \code{finteraction}, \code{na_rm}, \code{na_insert}, \code{all_identical}, \code{all_obj_equal}, \code{seq_row}, \code{seq_col}, \code{setDimnames}, \code{setRownames}, \code{setColnames}, \code{is.categorical}, \code{is.Date}, \code{as.numeric_factor}, \code{as.character_factor} \cr \cr \cr


Data and Global Macros \tab\tab Groningen Growth and Development Centre 10-Sector Database, World Bank World Development dataset, and some global macros containing links to the topical documentation pages (including this page), all exported objects (excluding exported S3 methods), all generic functions, the 2 datasets, all fast functions, all fast statistical (scalar-valued) functions, and all transformation operators / operator-like functions.  \tab\tab \code{\link{GGDC10S}, \link{wlddev}, .COLLAPSE_TOPICS, .COLLAPSE_ALL, .COLLAPSE_GENERIC, .COLLAPSE_DATA, .FAST_FUN, .FAST_STAT_FUN, .OPERATOR_FUN} \cr

}

}
\section{Details}{
The added top-level documentation infrastructure in \code{collapse} allows you to effectively navigate the package (as in other commercial software documentations like Mathematica). Calling \code{?FUN} brings up the documentation page documenting the function as in other R packages. You can also call topical documentation pages directly from the console. The links to these pages are contained in \code{.COLLAPSE_TOPICS} (i.e. calling \code{help(.COLLAPSE_TOPICS[1])} brings up this page).
}
\author{
\bold{Maintainer}: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}
}
\seealso{
  \link{collapse-package}
}
% \keyword{package}
\keyword{manip}
\keyword{documentation}

