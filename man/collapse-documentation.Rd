\name{collapse-documentation}
\alias{A0-collapse-documentation}
\alias{collapse-documentation}
\alias{.COLLAPSE_TOPICS}
\alias{.COLLAPSE_ALL}
\alias{.COLLAPSE_GENERIC}
\alias{.COLLAPSE_DATA}
% \docType{package}
\title{Collapse Documentation & Overview}
\description{
The following table fully summarizes the contents of \emph{\link{collapse}}. The documentation is structured hierarchically: This is the main overview page, linking to topical overview pages and associated function pages (unless functions are documented on the topic page). % Calling \code{?FUN} brings up the documentation page for \code{FUN}, with links to associated topic pages and closely related functions. % Calling \code{help(FUN)} still brings up the right / most relevant page documenting the function. %  % Functions with separate documentation entries (apart from the topic page) are linked. % Each topic further has its own overview page in the documentation. % , linking to functions % , i.e. only functions with separate pages are linked here
}

\section{Topics and Functions}{

\tabular{lllll}{ \emph{ Topic }  \tab\tab \emph{ Main Features / Keywords} \tab\tab \emph{ Functions }  \cr % \Sexpr{"\u200B"} \Sexpr{"\u200B"}

  \link[=fast-statistical-functions]{Fast Statistical Functions} \tab\tab Fast (grouped and weighted) statistical functions for vector, matrix, data frame and grouped data frames (class 'grouped_df', \emph{dplyr} compatible).  \tab\tab \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{fnth}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fnobs}}, \code{\link{fndistinct}} \cr \cr \cr

\link[=fast-grouping-ordering]{Fast Grouping and Ordering} \tab\tab Fast (ordered) groupings from vectors, data frames, lists. 'GRP' objects are efficient inputs for programming with \emph{collapse}'s fast functions. \code{fgroup_by} can attach them to a data frame, for fast dplyr-style grouped computations. Fast splitting of vectors based on 'GRP' objects. Fast radix-based ordering and hash-based grouping (the workhorses behind \code{GRP}). Fast matching (rows) and unique values/rows, group counts, factor generation, vector grouping, interactions, dropping unused factor levels, generalized run-length type grouping and grouping of integer sequences and time vectors.   % (to optimize different / repeated computations over the same groups).
\tab\tab \code{\link{GRP}}, \code{\link{as_factor_GRP}}, \code{\link{GRPN}}, \code{\link{GRPid}}, \code{\link{GRPnames}}, \code{\link{is_GRP}}, \code{\link{fgroup_by}}, \code{\link{fgroup_vars}}, \code{\link{fungroup}}, \code{\link{gsplit}}, \code{\link{greorder}}, \code{\link[=radixorder]{radixorder(v)}}, \code{\link{group}}, \code{\link{fmatch}}, \code{\link{ckmatch}}, \code{\link[=fmatch]{\%!in\%}}, \code{\link[=fmatch]{\%[!]iin\%}}, \code{\link{funique}},  \code{\link{fnunique}}, \code{\link{fduplicated}}, \code{\link{any_duplicated}}, \code{\link[=fcount]{fcount(v)}}, \code{\link{qF}}, \code{\link{qG}}, \code{\link{is_qG}}, \code{\link{finteraction}}, \code{\link{fdroplevels}}, \code{\link{groupid}}, \code{\link{seqid}}, \code{\link{timeid}} \cr \cr \cr
% \code{GRP} creates 'GRP' objects, and \code{fgroup_by} can be used to attach them to a data frame (analogous to \code{dplyr::group_by})

% (speed about 2x '[' for selecting and 4x '[<-' for replacing). %, get data, variables names, variable indices

\link[=fast-data-manipulation]{Fast Data Manipulation} \tab\tab Fast and flexible select, subset, summarise, mutate/transform, sort/reorder, combine, join, reshape, rename and relabel data. Some functions modify by reference and/or allow assignment. In addition a set of (standard evaluation) functions for fast selecting, replacing or adding data frame columns, including shortcuts to select and replace variables by data type.
\tab\tab \code{\link[=fselect]{fselect(<-)}}, \code{\link[=fsubset]{fsubset/ss}}, \code{\link{fsummarise}}, \code{\link{fmutate}}, \code{\link{across}}, \code{\link[=ftransform]{(f/set)transform(v)(<-)}}, \code{\link[=fcompute]{fcompute(v)}}, \code{\link[=roworder]{roworder(v)}}, \code{\link[=colorder]{colorder(v)}}, \code{\link{rowbind}}, \code{\link{join}}, \code{\link{pivot}}, \code{\link[=frename]{(f/set)rename}}, \code{\link[=relabel]{(set)relabel}}, \code{\link[=get_vars]{get_vars(<-)}}, \code{\link[=add_vars]{add_vars(<-)}}, \code{\link[=num_vars]{num_vars(<-)}}, \code{\link[=cat_vars]{cat_vars(<-)}}, \code{\link[=char_vars]{char_vars(<-)}}, \code{\link[=fact_vars]{fact_vars(<-)}}, \code{\link[=logi_vars]{logi_vars(<-)}}, \code{\link[=date_vars]{date_vars(<-)}} \cr \cr \cr

\link[=quick-conversion]{Quick Data Conversion} \tab\tab Quick conversions: data.frame <> data.table <> tibble <> matrix (row- or column-wise) <> list | array > matrix, data.frame, data.table, tibble | vector > factor, matrix, data.frame, data.table, tibble; and converting factors / all factor columns. \tab\tab \code{qDF}, \code{qDT}, \code{qTBL}, \code{qM}, \code{qF}, \code{mrtl}, \code{mctl}, \code{as_numeric_factor}, \code{as_character_factor} \cr \cr \cr

\link[=advanced-aggregation]{Advanced Data Aggregation} \tab\tab Fast and easy (weighted and parallelized) aggregation of multi-type data, with different functions applied to numeric and categorical variables. Custom specifications allow mappings of functions to variables + renaming. \tab\tab \code{collap(v/g)} \cr \cr \cr

\link[=data-transformations]{Data Transformations} \tab\tab Fast row- and column- arithmetic and (object preserving) apply functionality for vectors, matrices and data frames. Fast (grouped) replacing and sweeping of statistics (by reference) and (grouped and weighted) scaling / standardizing, (higher-dimensional) between- and within-transformations (i.e. averaging and centering), linear prediction and partialling out. %Additional methods for grouped_df (\emph{dplyr}) and pseries, pdata.frame (\emph{plm}).
\tab\tab \code{\link[=arithmetic]{\%(r/c)r\%}}, \code{\link[=arithmetic]{\%(r/c)(+/-/*//)\%}}, \code{\link{dapply}}, \code{\link{BY}}, \code{\link[=TRA]{(set)TRA}}, \code{\link[=fscale]{fscale/STD}}, \code{\link[=fbetween]{fbetween/B}}, \code{\link[=fwithin]{fwithin/W}}, \code{\link[=HDB]{fhdbetween/HDB}}, \code{\link[=HDW]{fhdwithin/HDW}} \cr \cr \cr

Linear Models \tab\tab Fast (weighted) linear model fitting with 6 different solvers and a fast F-test to test exclusion restrictions on linear models with (large) factors. \tab\tab \code{\link{flm}}, \code{\link{fFtest}} \cr \cr \cr

\link[=time-series-panel-series]{Time Series and Panel Series} \tab\tab Fast and class-agnostic indexed time series and panel data objects, check for irregularity in time series and panels, and efficient time-sequence to integer/factor conversion. Fast (sequences of) lags / leads and (lagged / leaded and iterated, quasi-, log-) differences, and (compounded) growth rates on (irregular) time series and panel data. Flexible cumulative sums. Panel data to array conversions. Multivariate panel- auto-, partial- and cross-correlation functions. %Additional methods for grouped_df (\emph{dplyr}) and pseries, pdata.frame (\emph{plm}).
\tab\tab
\code{\link{findex_by}}, \code{\link{findex}}, \code{\link{unindex}}, \code{\link{reindex}}, \code{\link{is_irregular}}, \code{\link{to_plm}}, \code{\link{timeid}},
\code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D/Dlog}}, \code{\link[=fgrowth]{fgrowth/G}}, \code{\link{fcumsum}}, \code{\link{psmat}}, \code{\link{psacf}}, \code{\link{pspacf}}, \code{\link{psccf}}  \cr \cr \cr

\link[=summary-statistics]{Summary Statistics} \tab\tab Fast (grouped and weighted) summary statistics for cross-sectional and panel data. Fast (weighted) cross tabulation. Efficient detailed description of data frame. Fast check of variation in data (within groups / dimensions). (Weighted) pairwise correlations and covariances (with obs. and p-value), pairwise observation count. %Some additional methods for grouped_df (\emph{dplyr}) pseries and pdata.frame (\emph{plm}).
\tab\tab \code{\link{qsu}}, \code{\link{qtab}}, \code{\link{descr}}, \code{\link{varying}}, \code{\link{pwcor}}, \code{\link{pwcov}}, \code{\link{pwnobs}} \cr \cr \cr \cr

Other Statistical \tab\tab Fast euclidean distance computations, (weighted) sample quantiles, and range of vector. \tab\tab \code{\link{fdist}}, \code{\link{fquantile}}, \code{\link{frange}} \cr \cr \cr

\link[=list-processing]{List Processing} \tab\tab (Recursive) list search and checks, extraction of list-elements / list-subsetting, fast (recursive) splitting, list-transpose, apply functions to lists of data frames / data objects, and generalized recursive row-binding / unlisting in 2-dimensions / to data frame.  \tab\tab \code{\link{is_unlistable}}, \code{\link{ldepth}}, \code{\link{has_elem}}, \code{\link{get_elem}}, \code{\link[=atomic_elem]{atomic_elem(<-)}}, \code{\link[=list_elem]{list_elem(<-)}}, \code{\link{reg_elem}}, \code{\link{irreg_elem}}, \code{\link{rsplit}}, \code{\link{t_list}}, \code{\link{rapply2d}}, \code{\link{unlist2d}}, \code{\link{rowbind}} \cr

%  (within- and between-groups) ww

%Visual Data Exploration \tab\tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab\tab histplot, densplot, scplot, corplot \cr

\link[=recode-replace]{Recode and Replace Values} \tab\tab Recode multiple values (exact or regex matching) and replace \code{NaN/Inf/-Inf} and outliers (according to 1- or 2-sided threshold or standard-deviations) in vectors, matrices or data frames. Insert a value at arbitrary positions into vectors, matrices or data frames. \tab\tab \code{recode_num}, \code{recode_char}, \code{replace_na}, \code{replace_inf}, \code{replace_outliers}, \code{\link{pad}} \cr \cr \cr


\link[=efficient-programming]{(Memory) Efficient Programming} \tab\tab Efficient comparisons of a vector/matrix with a value, and replacing values/rows in vector/matrix/DF (avoiding logical vectors or subsets), faster generation of initialized vectors, and fast mathematical operations on vectors/matrices/DF's with no copies at all. Fast missing value detection, (random) insertion and removal, lengths and C storage types, greatest common divisor of vector, \code{nlevels} for factors, \code{nrow}, \code{ncol}, \code{dim} (for data frames) and \code{seq_along} rows or columns. Fast vectorization of matrices and lists, and choleski inverse of symmetric PD matrix. \tab\tab
\code{anyv}, \code{allv}, \code{allNA}, \code{whichv}, \code{whichNA}, \code{\%==\%},
\code{\%!=\%}, \code{copyv}, \code{setv}, \code{alloc}, \code{setop}, \code{\%+=\%}, \code{\%-=\%}, \code{\%*=\%}, \code{\%/=\%}, \code{missing_cases}, \code{na_insert}, \code{na_rm}, \code{na_omit}, \code{vlengths}, \code{vtypes}, \code{vgcd}, \code{fnlevels}, \code{fnrow}, \code{fncol}, \code{fdim}, \code{seq_row}, \code{seq_col}, \code{vec}, \code{cinv} \cr \cr \cr

\link[=small-helpers]{Small (Helper) Functions} \tab\tab Multiple-assignment, non-standard concatenation, set and extract variable labels and classes, display variable names and labels together, add / remove prefix or postfix to / from column names, check exact or near / numeric equality of multiple objects or of all elements in a list, get names of functions called in an expression, return object with dimnames, row- or colnames efficiently set, or with all attributes removed, C-level functions to set and shallow-copy attributes, identify categorical (non-numeric) and date(-time) objects. \tab\tab \code{massign}, \code{\%=\%}, \code{.c}, \code{vlabels(<-)}, \code{setLabels}, \code{vclasses}, \code{namlab}, \code{add_stub}, \code{rm_stub}, \code{all_identical}, \code{all_obj_equal}, \code{all_funs}, \code{setDimnames}, \code{setRownames}, \code{setColnames}, \code{unattrib}, \code{setAttrib}, \code{setattrib}, \code{copyAttrib}, \code{copyMostAttrib}, \code{is_categorical}, \code{is_date} \cr \cr \cr


Data and Global Macros \tab\tab Groningen Growth and Development Centre 10-Sector Database, World Bank World Development dataset, and some global macros containing links to the topical documentation pages (including this page), all exported objects (excluding exported S3 methods and depreciated functions), all generic functions (excluding depreciated), the 2 datasets, depreciated functions, all fast functions, all fast statistical (scalar-valued) functions, and all transformation operators (these are not infix functions but function shortcuts resembling operators in a statistical sense, such as the lag/lead operators \code{L}/\code{F}, both wrapping \code{flag}, see \code{\link{.OPERATOR_FUN}}).  \tab\tab \code{\link{GGDC10S}, \link{wlddev}, .COLLAPSE_TOPICS, .COLLAPSE_ALL, .COLLAPSE_GENERIC, .COLLAPSE_DATA, .COLLAPSE_OLD, .FAST_FUN, .FAST_STAT_FUN, .OPERATOR_FUN} \cr\cr\cr


\link[=collapse-options]{Package Options} \tab\tab \code{set_collapse}/\code{get_collapse} can be used to globally set/get the defaults for \code{na.rm}, \code{nthreads} and \code{sort}, etc., arguments found in many functions, and to globally control the namespace with options 'mask' and 'remove': 'mask' can be used to mask base R/dplyr functions by export copies of equivalent \emph{collapse} functions starting with \code{"f"}, removing the leading \code{"f"} (e.g. exporting \code{subset <- fsubset}). 'remove' allows removing arbitrary functions from the exported namespace. \code{options("collapse_unused_arg_action")} sets the action taken by generic statistical functions when unknown arguments are passed to a method. The default is \code{"warning"}.  \tab\tab \code{set_collapse}, \code{get_collapse} \cr\cr\cr

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

