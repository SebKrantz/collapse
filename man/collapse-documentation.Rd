\name{A0-collapse-documentation}
\alias{A0-collapse-documentation}
\alias{collapse-documentation}
\alias{.COLLAPSE_TOPICS}
\alias{.COLLAPSE_ALL}
\alias{.COLLAPSE_GENERIC}
\alias{.COLLAPSE_DATA}
% \docType{package}
\title{Collapse Documentation & Overview}
\description{
The following table fully summarizes the contents of \emph{\link{collapse}}. The documentation is structured hierarchically: This is the main overview page, linking to topical overview pages and associated function pages (unless functions are documented on the topic page). % Calling \code{?FUN} brings up the documentation page for \code{FUN}, with links to associated topic pages and closely related functions. % Calling \code{help(FUN)} still brings up the right / most relevant page documenting the function. %  % Functions with separate documentation entries (apart from the topic page) are linked. % Each topic further has it's own overview page in the documentation. % , linking to functions % , i.e. only functions with separate pages are linked here
}

\section{Topics and Functions}{

\tabular{lllll}{ \emph{ Topic }  \tab\tab \emph{ Main Features / Keywords} \tab\tab \emph{ Functions }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"}

  \link[=A1-fast-statistical-functions]{Fast Statistical Functions} \tab\tab Fast (grouped and weighted) statistical functions for vector, matrix, data frame and grouped data frames (class 'grouped_df', \emph{dplyr} compatible).  \tab\tab \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{fnth}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fNobs}}, \code{\link{fNdistinct}} \cr \cr \cr

\link[=A2-fast-grouping-ordering]{Fast Grouping and Ordering} \tab\tab Fast (ordered) groupings from vectors, data frames, lists. 'GRP' objects are extremely efficient inputs for programming with \emph{collapse}'s fast functions. \code{fgroup_by} can attach them to a data frame, for fast dplyr-style grouped computations. In addition fast radix-sort based ordering, unique values/rows, factor generation, vector grouping, interactions, generalized run-length type grouping and grouping of time-sequences.   % (to optimize different / repeated computations over the same groups).
\tab\tab \code{\link{GRP}}, \code{\link{as.factor_GRP}}, \code{\link{GRPnames}}, \code{\link{is.GRP}}, \code{\link{fgroup_by}}, \code{\link{fgroup_vars}}, \code{\link{fungroup}}, \code{\link[=radixorder]{radixorder(v)}}, \code{\link{funique}}, \code{\link{qF}}, \code{\link{qG}}, \code{\link{is.qG}}, \code{\link{fdroplevels}}, \code{\link{finteraction}}, \code{\link{groupid}}, \code{\link{seqid}} \cr \cr \cr
% \code{GRP} creates 'GRP' objects, and \code{fgroup_by} can be used to attach them to a data frame (analogous to \code{dplyr::group_by})

% (speed about 2x '[' for selecting and 4x '[<-' for replacing). %, get data, variables names, variable indices

\link[=A3-fast-data-manipulation]{Fast Data Manipulation} \tab\tab Fast and flexible select, subset, summarise, transform, sort/reorder and rename data, including modifying/adding columns by reference, automated replacing/adding with lists of transformed columns, and computing columns saved as a new dataset. In addition a set of (standard evaluation) functions for fast selecting, replacing or adding data frame columns, including shortcuts to select and replace variables by data type.
\tab\tab \code{\link[=fselect]{fselect(<-)}}, \code{\link[=fsubset]{fsubset/ss}}, \code{\link{fsummarise}}, \code{\link[=ftransform]{(f/set)transform(v)(<-)}}, \code{\link{fcompute}}, \code{\link[=roworder]{roworder(v)}}, \code{\link[=colorder]{colorder(v)}}, \code{\link[=frename]{(f/set)rename}}, \code{\link[=get_vars]{get_vars(<-)}}, \code{\link[=add_vars]{add_vars(<-)}}, \code{\link[=num_vars]{num_vars(<-)}}, \code{\link[=cat_vars]{cat_vars(<-)}}, \code{\link[=char_vars]{char_vars(<-)}}, \code{\link[=fact_vars]{fact_vars(<-)}}, \code{\link[=logi_vars]{logi_vars(<-)}}, \code{\link[=Date_vars]{Date_vars(<-)}} \cr \cr \cr

\link[=A4-quick-conversion]{Quick Data Conversion} \tab\tab Quick conversions: data.frame <> data.table <> tibble | matrix <> list, data.frame, data.table (row- or column- wise), tibble | array > matrix, data.frame, data.table, tibble | list > data.frame, data.table, tibble | vector > factor, matrix, data.frame, data.table, tibble; and converting factors / all factor columns. \tab\tab \code{qDF}, \code{qDT}, \code{qTBL}, \code{qM}, \code{qF}, \code{mrtl}, \code{mctl}, \code{as.numeric_factor}, \code{as.character_factor} \cr \cr \cr

\link[=A5-advanced-aggregation]{Advanced Data Aggregation} \tab\tab Fast and easy (weighted and parallelized) aggregation of multi-type data, with (multiple) functions applied to numeric and categorical columns. Also supports fully customized aggregation tasks mapping functions to columns. \tab\tab \code{collap(v/g)} \cr \cr \cr

\link[=A6-data-transformations]{Data Transformations} \tab\tab Fast row- and column- arithmetic and (object preserving) apply functionality for matrices and data frames. Fast (grouped) replacing and sweeping of statistics and (grouped and weighted) scaling / standardizing, (higher-dimensional) within- and between-transformations (i.e. centering and averaging), linear prediction and partialling out. Additional methods for grouped_df (\emph{dplyr}) and pseries, pdata.frame (\emph{plm}). \tab\tab \code{\link[=arithmetic]{\%(r/c)r\%}}, \code{\link[=arithmetic]{\%(r/c)(+/-/*//)\%}}, \code{\link{dapply}}, \code{\link{BY}}, \code{\link{TRA}}, \code{\link[=fscale]{fscale/STD}}, \code{\link[=fbetween]{fbetween/B}}, \code{\link[=fwithin]{fwithin/W}}, \code{\link[=HDB]{fHDbetween/HDB}}, \code{\link[=HDW]{fHDwithin/HDW}} \cr \cr \cr

Linear Models \tab\tab Fast (weighted) linear model fitting with 6 different solvers and a fast F-test to test exclusion restrictions on linear models with (large) factors. \tab\tab \code{\link{flm}}, \code{\link{fFtest}} \cr \cr \cr

\link[=A7-time-series-panel-series]{Time Series and Panel Series} \tab\tab Fast (sequences of) lags / leads and (lagged / leaded and iterated) differences, quasi-differences, (quasi-) log-differences and (compounded) growth rates on (unordered, irregular) time series and panel data. Panel data to (ts-)array conversions. Multivariate panel- auto-, partial- and cross-correlation functions. Additional methods for grouped_df (\emph{dplyr}) and pseries, pdata.frame (\emph{plm}). \tab\tab \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D/Dlog}}, \code{\link[=fgrowth]{fgrowth/G}}, \code{\link{psmat}}, \code{\link{psacf}}, \code{\link{pspacf}}, \code{\link{psccf}}  \cr

\link[=A8-list-processing]{List Processing} \tab\tab (Recursive) list search and identification, search and extract list-elements / list-subsetting, splitting, apply functions to lists of data frames / data objects, and (fast) generalized recursive row-binding / unlisting in 2-dimensions / to data frame.  \tab\tab \code{\link{is.regular}}, \code{\link{is.unlistable}}, \code{\link{ldepth}}, \code{\link{has_elem}}, \code{\link{get_elem}}, \code{\link[=atomic_elem]{atomic_elem(<-)}}, \code{\link[=list_elem]{list_elem(<-)}}, \code{\link{reg_elem}}, \code{\link{irreg_elem}}, \code{\link{rsplit}}, \code{\link{rapply2d}}, \code{\link{unlist2d}} \cr

\link[=A9-summary-statistics]{Summary Statistics} \tab\tab Fast (grouped and weighted), summary statistics for cross-sectional and complex multilevel / panel data. Efficient detailed description of data frame. Fast check of variation in data (within groups / dimensions). (Weighted) pairwise correlations and covariances (with observation count, p-value and pretty printing), pairwise observation count. Some additional methods for grouped_df (\emph{dplyr}) pseries and pdata.frame (\emph{plm}). \tab\tab \code{\link{qsu}}, \code{\link{descr}}, \code{\link{varying}}, \code{\link{pwcor}}, \code{\link{pwcov}}, \code{\link{pwNobs}} \cr \cr \cr \cr

%  (within- and between-groups) ww

%Visual Data Exploration \tab\tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab\tab histplot, densplot, scplot, corplot \cr

\link[=AA1-recode-replace]{Recode and Replace Values} \tab\tab Recode multiple values (exact or regex matching) and replace \code{NaN/Inf/-Inf} and outliers (according to 1- or 2-sided threshold or column standard-deviation) in vectors, matrices or data frames. \tab\tab \code{recode_num}, \code{recode_char}, \code{replace_NA}, \code{replace_Inf}, \code{replace_outliers} \cr \cr \cr

\link[=AA2-small-helpers]{Small (Helper) Functions} \tab\tab Fast missing value detection, insertion and removal, faster nlevels for factors, fast nrow, ncol, dim (for data frames) and \code{seq_along} rows or columns, non-standard concatenation, set and extract variable labels, extract variable classes and C storage types, display variable names, labels and classes, add / remove prefix or postfix to / from column names, not-in operator, matching with error message for non-matched, check exact or near / numeric equality of multiple objects or of all elements in a list, return object with dimnames, row- or colnames efficiently set, or with all attributes removed, C-level functions to set and duplicate / copy attributes, identify categorical and date(-time) objects, Choleski (fast) inverse of symmetric PD matrix. \tab\tab \code{allNA}, \code{missing_cases}, \code{na_insert}, \code{na_rm}, \code{na_omit}, \code{fnlevels}, \code{fnrow}, \code{fncol}, \code{fdim}, \code{seq_row}, \code{seq_col}, \code{.c}, \code{vlabels(<-)}, \code{vclasses}, \code{vtypes}, \code{namlab}, \code{add_stub}, \code{rm_stub}, \code{\%!in\%}, \code{ckmatch}, \code{all_identical}, \code{all_obj_equal}, \code{setDimnames}, \code{setRownames}, \code{setColnames}, \code{unattrib}, \code{setAttrib}, \code{copyAttrib}, \code{copyMostAttrib}, \code{is.categorical}, \code{is.Date}, \code{cinv} \cr \cr \cr


Data and Global Macros \tab\tab Groningen Growth and Development Centre 10-Sector Database, World Bank World Development dataset, and some global macros containing links to the topical documentation pages (including this page), all exported objects (excluding exported S3 methods), all generic functions, the 2 datasets, all fast functions, all fast statistical (scalar-valued) functions, and all transformation operators / operator-like functions.  \tab\tab \code{\link{GGDC10S}, \link{wlddev}, .COLLAPSE_TOPICS, .COLLAPSE_ALL, .COLLAPSE_GENERIC, .COLLAPSE_DATA, .FAST_FUN, .FAST_STAT_FUN, .OPERATOR_FUN} \cr

\link[=collapse-options]{Package Options} \tab\tab Set the action taken by generic functions encountering unknown arguments. The default is \code{"warning"}. Other choices are \code{"error"}, \code{"message"} or \code{"none"}, where \code{"none"} enables silent swallowing. \tab\tab \code{options(collapse_} \code{unused_arg_action)} \cr

}

}
\section{Details}{
The added top-level documentation infrastructure in \emph{collapse} allows you to effectively navigate the package. % (as in other commercial software documentations like Mathematica).
Calling \code{?FUN} brings up the documentation page documenting the function, which contains links to associated topic pages and closely related functions. You can also call topical documentation pages directly from the console. The links to these pages are contained in the global macro \code{.COLLAPSE_TOPICS} (e.g. calling \code{help(.COLLAPSE_TOPICS[1])} brings up this page).
}
\author{
\bold{Maintainer}: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}
}
\seealso{
  \link{collapse-package}
}
% \keyword{package}
\keyword{documentation}

