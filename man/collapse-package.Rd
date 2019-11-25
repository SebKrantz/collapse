\name{collapse-package}
\alias{collapse-package}
\alias{collapse}
\docType{package}
\title{ %   \code{collapse:}
Advanced and Fast Data Transformation and Exploration % in R % and Fast  % Advanced Data Aggregation, Transformation and Exploration in R % Fast Advanced Data Transformation and Exploration % Advanced and Fast % in R % Advanced Data Aggregation, Transformation and Exploration in R % or: Advanced and Fast Data Manipulation and Exploration. Or: Advanced and Fast Data Transformation and Exploration.. Or: Advanced Data Manipulation and Exploration. Or: Advanced and Fast Data Aggregation, Transformation and Exploration. Or: Fast Advanced Data Manipulation and Exploration. (Nah, 2 adjectives is bad, and you have to earn the title Fast !!) % Advanced and Fast Data Manipulation and Exploration (Manipulation sounds so meticulous, but your package moves in broad strokes)
  % Fazit: The title is Great !!
}
\description{
\code{collapse} is a C/C++ based package for data manipulation in R. It's aims are 1. to facilitate complex data transformation and exploration tasks and 2. to help make R code fast, flexible, parsimonious and programmer friendly. % \code{collapse} is a fast  %to facilitate (advanced) data manipulation in R   % To achieve the latter,
% collapse provides a broad set.. -> Nah, it's not a misc package
}
\details{
% Put this below bullet points ?? % \code{collapse} provides a carefully conceived
\code{collapse} provides a compact set of functions organized into several topics (see below). Most are S3 generic with methods for vectors, matrices and data.frames. Inputs are quickly passed to compiled C/C++ code, enabling flexible and parsimonious coding at extreme speeds.

The package avoids lazy-evaluation and exports it's core methods, for maximum programmability.  % Most are S3 generic with methods for common \code{R} objects (vectors, matrices, data.frames, ...) % high computation  %(aggregation and transformations ~10x \code{data.table} on data <1 Mio obs.).
Smart attribute handling and additional methods ensure compatibility and support for \code{dplyr}, \code{data.table} and the \code{plm} approach to panel-data. % and \code{base/stats}. % extra methods warrant, provide ?

%Beyond speed, flexibility and parsimony in coding, a central objective of \code{collapse} is to facilitate advanced / complex operations on data.
\bold{Key functionality:}     (More details in table below) % Key areas Key topics addressed by \code{collapse} are: % where \code{collapse} offers innovative solutions are:

\tabular{ll}{ % \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
                 (1) \tab \emph{Advanced data programming}: A set of Fast statistical functions supporting groupwise and weighted computations on vectors, matrices and data.frames. Quick (ordered) and reusable grouping, quick data-converions, and quick select and replace variables. \cr
                 (2) \tab \emph{Advanced aggregation}: Fast and easy multi-data-type, multi-function, weighted, parallelized and fully customized data aggregation. \cr
                 (3) \tab \emph{Advanced transformations}: Efficient (grouped, weighted) sweeping out of statistics, scaling, centering, higher-dimensional centering, complex linear prediction and partialling-out. \cr
                 (4) \tab \emph{Advanced time-computations}: Efficient (sequences of) lags/leads, and (iterated) differences and growth rates on (unordered) time-series and panel data. Auto, partial and cross-correlation functions for panel data. Panel data to (ts-)array conversions. \cr
                 (5) \tab \emph{List Processing}: Recursive list-extraction/subsetting, apply, and row-binding/unlisting in 2D. \cr
                 (6) \tab \emph{Advanced data exploration}: Fast (grouped, weighted) summary statistics for cross-sectional and panel data. \cr % Elegant (multivariate, grouped, weighted and panel-decomposed) histograms, density, scatter and correlation plots.  \cr
}
\code{collapse} is built using \code{Rcpp} and imports \code{C} functions from \code{data.table} and \code{lfe}. %, and uses \code{ggplot2} visualizations.

}
\section{Functions and Topics}{
The following table fully summarizes the content of \code{collapse}. The documentation follows a hierarchical structure - this is the top page, linking to topical pages, linking to functions (unless already fully documented on the topic page). % Functions with separate documentation entries (apart from the topic page) are linked. % Each topic further has it's own overview page in the documentation.

\tabular{lll}{ \emph{ Topic } \Sexpr{"\u200B"} \Sexpr{"\u200B"} \tab \emph{ Main Features / Keywords} \tab \emph{ Functions }  \cr

  \link{Fast Statistical Functions} \tab Fast (grouped and weighted) statistical functions for vector, matrix, data.frame and grouped_df (dplyr compatible).  \tab \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fNobs}}, \code{\link{fNdistinct}} \cr

\link{Quick Grouping} \tab Fast (ordered or unordered) and reusable groupings from vectors, data.frames, lists. \tab \code{GRP} \cr

\link{Quick Select and Replace Variables} \tab Quickly and securely select and replace variables from data.frames / data.tables / tibbles / grouped tibbles / pdata.frame etc... %, get data, variables names, variable indices
\tab \code{get_vars}, \code{num_vars}, \code{cat_vars}, \code{char_vars}, \code{fact_vars}, \code{logi_vars}, \code{Date_vars} \cr

\link{Quick Data Converion} \tab Fast conversions, vector > factor / matrix / data.frame / data.table, matrix > list / data.frame / data.table, array > matrix / data.frame / data.table, list > data.frame or data.table. \tab \code{qDF}, \code{qDT}, \code{qM}, \code{qF}, \code{qG}, \code{mrtl}, \code{mctl} \cr \cr \cr

\link{Advanced Data Aggregation} \tab Fast and easy aggregation of multi-type data, using multiple-functions, and fully customizeable. Can perform weighted and parallelized aggregations. \tab \code{collap} \cr \cr \cr

\link{Data Transformations} \tab Efficient row- column- data-apply, Split-Apply-Combine computing, fast (grouped and weighted) replace, sweep, scale / standardize, (higher-dimensional) within- and between-transformations (e.g. centering and averaging), complex lin. prediction and partialling out \tab \code{\link{dapply}}, \code{\link{BY}}, \code{\link{TRA}}, \code{\link[=fscale]{fscale/STD}}, \code{\link{W}}, \code{\link{B}}, \code{\link{HDW}}, \code{\link{HDB}}. \cr \cr \cr

\link{Time-Series and Panel-Series} \tab Efficiently compute (sequences of), lags / leads and (iterated) differences and growth rates / log-differences on (unordered) time-series and panel-data. Panel-data > (ts-)array conversions, panel- auto, partial- and cross - correlation functions \tab \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}}, \code{\link[=fgrowth]{fgrowth/G}}, \code{\link{psmat}}, \code{\link{psacf}}, \code{\link{pspacf}}, \code{\link{psccf}}  \cr \cr \cr

\link[=qsu]{Summary Statistics} \tab Extremely fast (one-pass, grouped and weighted), summary statistics for cross-sectional and panel-data (within- and between-groups) \tab \code{qsu} \cr

List Processing \tab (Recursive) row-binding / unlisting in 2-dimensions, data-apply, search and extract list-elements / list-subsetting, list-identification.  \tab \code{\link{unlist2d}}, \code{\link{rapply2d}}, \code{\link{get_elem}}, \code{\link{has_elem}}, \code{\link{atomic_elem}}, \code{\link{list_elem}}, \code{\link{reg_elem}}, \code{\link{irreg_elem}}, \code{\link{is.regular}}, \code{\link{is.unlistable}}, \code{\link{ldepth}} \cr
%Visual Data Exploration \tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab histplot, densplot, scplot, corplot \cr
Small Functions \tab set and extract variable labels, return object with dimnames, row- or col- names set, seq row col, pairwise cov cor, vector remove NA's, check all elements identical, recode values, replace NA/NAN/Inf,-Inf. \tab \code{\link{vlabels}}, \code{\link{namlab}}, \code{\link{seq_row}}, \code{\link{seq_col}}, \code{\link{setRownames}}, \code{\link{setColnames}}, \code{\link{setDimnames}}, \code{\link{pwcor}}, \code{\link{pwcov}}, \code{\link{na.rm}}, \code{\link{all.identical}}, Recode, replace_non_finite \cr
}

}

\author{
Sebastian Krantz

Maintainer: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}
}
%\references{
%  This optional section can contain literature or other references for
%  background information.
%}
\keyword{ package }
%\seealso{
%  Optional links to other man pages
%}
\examples{
## To Do: A Taste of collapse:
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.
  }
}
