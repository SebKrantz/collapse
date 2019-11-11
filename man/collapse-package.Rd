\name{collapsedev-package}
\alias{collapsedev-package}
\alias{collapsedev}
\docType{package}
\title{
  \code{collapse:} Advanced Data Transformation and Exploration in R % Advanced Data Aggregation, Transformation and Exploration in R % Fast Advanced Data Transformation and Exploration % Advanced and Fast % in R % Advanced Data Aggregation, Transformation and Exploration in R % or: Advanced and Fast Data Manipulation and Exploration. Or: Advanced and Fast Data Transformation and Exploration.. Or: Advanced Data Manipulation and Exploration. Or: Advanced and Fast Data Aggregation, Transformation and Exploration. Or: Fast Advanced Data Manipulation and Exploration. (Nah, 2 adjectives is bad, and you have to earn the title Fast !!) % Advanced and Fast Data Manipulation and Exploration (Manipulation sounds so meticulous, but your package moves in broad strokes)
  % Fazit: The title is Great !!
}
\description{
\code{collapse} is a C/C++ based package for data manipulation in R. It's aims are 1. to facilitate complex data transformation and exploration tasks and 2. to help make R code fast, flexible, parsimonious and programmer friendly. % \code{collapse} is a fast  %to facilitate (advanced) data manipulation in R   % To achieve the latter,
% collapse provides a broad set.. -> Nah, it's not a misc package
}
\details{
% Put this below bullet points ?? % \code{collapse} provides a carefully conceived
\code{collapse} provides a compact set of functions organized into several topics (see below). Most are S3 generic with methods for vectors, matrices and data.frames. Inputs are quickly passed to C/C++ functions, enabling flexible and parsimonious code at extreme speeds.

\code{collapse} consciously avoids lazy-evaluation or complex input formats and exports key methods, for maximum programmability.  % Most are S3 generic with methods for common \code{R} objects (vectors, matrices, data.frames, ...) % high computation  %(aggregation and transformations ~10x \code{data.table} on data <1 Mio obs.).
Smart attribute handling and extra methods ensure compatibility with \code{dplyr}, \code{data.table} and the \code{plm} approach to panel-data. % and \code{base/stats}. % extra methods warrant, provide ?

%Beyond speed, flexibility and parsimony in coding, a central objective of \code{collapse} is to facilitate advanced / complex operations on data.
\bold{Key functionality added by \code{collapse}:}     (More details in table below) % Key areas Key topics addressed by \code{collapse} are: % where \code{collapse} offers innovative solutions are:

\tabular{ll}{ % \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
                 (1) \tab \emph{Advanced data programming}: A large set of fast statistical functions supporting groupwise and weighted computations on vectors, matrices and data.frames. Quick (ordered) and reusable grouping, quick data-converions, and quick select and replace variables. \cr
                 (2) \tab \emph{Advanced aggregation}: Easy multi-data-type, multi-function, weighted, parallelized and highly customized data aggregation. \cr
                 (3) \tab \emph{Advanced transformations}: Efficient (grouped, weighted) sweeping out of statistics, scaling, centering, higher-dimensional centering, complex linear prediction and partialling-out. \cr
                 (4) \tab \emph{Advanced time-computations}: Efficient (sequences of) lags/leads, and (iterated) differences and growth rates on (unordered) time-series and panel data. Auto, partial and cross-correlation functions for panel data. Panel data to (ts-)matrix/3D-array expansion. \cr
                 (5) \tab \emph{List Processing}: Recursive list-extraction/subsetting, apply, and row-binding/unlisting in 2D. \cr
                 (6) \tab \emph{Advanced data exploration}: Fast (grouped, weighted) summary statistics for cross-sectional and panel data. Elegant (multivariate, grouped, weighted and panel-decomposed) histograms, density, scatter and correlation plots.  \cr
}

\code{collapse} is built using \code{Rcpp}, imports \code{C} functions from \code{data.table} and \code{lfe}, and uses \code{ggplot2} visualizations.

}
\section{Functions and Topics}{
The following table fully summarizes the content of \code{collapse}. Each topic further has it's own overview page in the documentation.
\tabular{lll}{ \emph{ Topic } \Sexpr{"\u200B"} \Sexpr{"\u200B"} \tab \emph{ Main Features / Keywords} \tab \emph{ Functions }  \cr
  Fast Statistical Functions \tab Fast (grouped and weighted) statistical functions for vector, matrix, data.frame and grouped_df (dplyr compatible).  \tab fsum, fprod, fmean, fmedian, fmode, fvar, fsd, fmin, fmax, ffirst, flast, fNobs, fNdistinct \cr
Quick Grouping \tab Fast (ordered or unordered) and reusable groupings from vectors, data.frames, lists \tab GRP \cr
Select / Replace Variables \tab select and replace variables by type, using a function or a regular expression. %, get data, variables names, variable indices
\tab get_vars, num_vars, cat_vars, char_vars, fact_vars, logi_vars, Date_vars \cr
Quick Data-Converions \tab Fast conversions, vector > factor / matrix / data.frame or table, matrix > list / data.frame or table, array > matrix / data.frame or table, list > data.frame or table \tab qDF, qDT, qM, qF, qG, mrtl, mctl \cr
Data Aggregation \tab Fast Aggregate multi-data-type, multi-function, weighted, fully custom, parallelized. Efficent (parallelized) split-apply-combine computing. \tab collap, BY \cr
Data Transformations \tab Row- column- data-apply. Fast (grouped and weighted) replace, sweep, scale / standardize, (higher-dimensional) within- and between-transformations (e.g. centering and averaging), complex lin. prediction and partialling out \tab dapply, TRA, BY, fscale / STD, W, B, HDW, HDB. \cr
Time-Series / Panel-Series \tab Efficiently compute (sequences of), lags / leads and (iterated) differences and growth rates / log-differences on (unordered) time-series and panel-data. Panel-data > (ts-)array conversions, panel- auto, partial- and cross - correlation functions \tab flag / L, flead / F, fdiff / D, fgrowth / G, psmat, psacf, pspacf, psccf  \cr
Summary Statistics \tab Extremely fast (one-pass, grouped and weighted), summary statistics for cross-sectional and panel-data (within- and between-groups) \tab qsu \cr
List Processing \tab (Recursive) row-binding / unlisting in 2-dimensions, data-apply, search and extract list-elements / list-subsetting, list-identification.  \tab unlist2d, rapply2d, get_elem, has_elem, atomic_elem, list_elem, reg_elem, irreg_elem, is.regular, is.unlistable, ldepth \cr
Visual Data Exploration \tab Pretty (grouped, weighted, and panel-decomposed) histogram-, density-, scatter- and correlation plots \tab histplot, densplot, scplot, corplot \cr
Small Functions \tab set and extract variable labels, return object with attributes set/duplicated and dimnames, row- or col- names set, seq row col, pairwise cov cor, fast vector remove NA's, faster nlevels, check all elements identical, recode values, replace NA/NAN/Inf,-Inf.  \tab vlabels, namlab, seq_row, seq_col, setRownames, setColnames, setDimnames, setAttributes, duplAttributes, pwcor, pwcov, na_rm, all.identical, Recode, replace_non_finite \cr
}

}

\author{
Sebastian Krantz

Maintainer: Sebastian Krantz <sebastian.krantz@graduateinstitute.ch>
}
\references{
  This optional section can contain literature or other references for
  background information.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
## To Do: A Taste of collapse:
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.
  }
}
