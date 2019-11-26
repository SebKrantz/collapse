\name{collapse-package}
\alias{collapse-package}
\alias{collapse}
\docType{package}
\title{ %   \code{collapse:}
Advanced and Fast Data Transformation % and Exploration % in R % and Fast  % Advanced Data Aggregation, Transformation and Exploration in R % Fast Advanced Data Transformation and Exploration % Advanced and Fast % in R % Advanced Data Aggregation, Transformation and Exploration in R % or: Advanced and Fast Data Manipulation and Exploration. Or: Advanced and Fast Data Transformation and Exploration.. Or: Advanced Data Manipulation and Exploration. Or: Advanced and Fast Data Aggregation, Transformation and Exploration. Or: Fast Advanced Data Manipulation and Exploration. (Nah, 2 adjectives is bad, and you have to earn the title Fast !!) % Advanced and Fast Data Manipulation and Exploration (Manipulation sounds so meticulous, but your package moves in broad strokes)
  % Fazit: The title is Great !!
}
\description{
\code{collapse} is a C/C++ based package for data manipulation in R. It's aims are 1. to facilitate complex data transformation and exploration tasks and 2. to help make R code fast, flexible, parsimonious and programmer friendly. % \code{collapse} is a fast  %to facilitate (advanced) data manipulation in R   % To achieve the latter,
% collapse provides a broad set.. -> Nah, it's not a misc package
}
\details{
% Put this below bullet points ?? % \code{collapse} provides a carefully conceived
\code{collapse} provides a compact set of functions organized into several topics (see \link[=collapse-documentation]{Collapse Overiew}). Most are S3 generic with methods for vectors, matrices and data.frames. Inputs are quickly passed to compiled C/C++ code, enabling flexible and parsimonious coding at extreme speeds.

The package avoids lazy-evaluation and exports core methods for maximum programmability.  % Most are S3 generic with methods for common \code{R} objects (vectors, matrices, data.frames, ...) % high computation  %(aggregation and transformations ~10x \code{data.table} on data <1 Mio obs.).
Smart attribute handling and additional methods ensure compatibility and support for \code{dplyr}, \code{data.table} and the \code{plm} approach to panel-data. % and \code{base/stats}. % extra methods warrant, provide ?

%Beyond speed, flexibility and parsimony in coding, a central objective of \code{collapse} is to facilitate advanced / complex operations on data.
\bold{Key functionality:}     (More details in \link[=collapse-documentation]{Collapse Overview}) % Key areas Key topics addressed by \code{collapse} are: % where \code{collapse} offers innovative solutions are:

\tabular{ll}{ % \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
                 (1) \tab \emph{Advanced data programming}: A set of Fast statistical functions supporting groupwise and weighted computations on vectors, matrices and data.frames. Quick (ordered) and reusable grouping, quick data-converions, and quick select and replace variables. \cr
                 (2) \tab \emph{Advanced aggregation}: Fast and easy multi-data-type, multi-function, weighted, parallelized and fully customized data aggregation. \cr
                 (3) \tab \emph{Advanced transformations}: Efficient (grouped, weighted) sweeping out of statistics, scaling, centering, higher-dimensional centering, complex linear prediction and partialling-out. \cr
                 (4) \tab \emph{Advanced time-computations}: Efficient (sequences of) lags/leads, and (iterated) differences and growth rates on (unordered) time-series and panel data. Auto, partial and cross-correlation functions for panel data. Panel data to (ts-)array conversions. \cr
                 (5) \tab \emph{List Processing}: Recursive list-identification, extraction/subsetting, apply, and row-binding/unlisting in 2D. \cr
                 (6) \tab \emph{Advanced data exploration}: Fast (grouped, weighted) summary statistics for cross-sectional and panel data. \cr % Elegant (multivariate, grouped, weighted and panel-decomposed) histograms, density, scatter and correlation plots.  \cr
}
\code{collapse} is built using \code{Rcpp} and imports \code{C} functions from \code{data.table} and \code{lfe}. %, and uses \code{ggplot2} visualizations.

}
\section{Getting Started}{
Please see \link[=collapse-documentation]{Collapse Documentation & Overview}, or one of the Vignettes.
}
\author{
Sebastian Krantz

Maintainer: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}
}
\section{Developing / Feature Requests / Bug Reporting}{
\itemize{
\item If you are interested in extending or optimizing this package, see the source code at \url{https://github.com/SebKrantz/collapse/tree/master}, fork and send pull-requests, or e-mail me.
\item Please send feature requests via e-mail.
\item Please report issues at \url{https://github.com/SebKrantz/collapse/issues} or e-mail me.
}
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
