\name{quick-conversion}
\alias{A4-quick-conversion}
\alias{quick-conversion}
\alias{qDF}
\alias{qDT}
\alias{qTBL}
\alias{qM}
\alias{mctl}
\alias{mrtl}
\alias{as_numeric_factor}
\alias{as_character_factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quick Data Conversion}
\description{
Fast, flexible and precise conversion of common data objects, without method dispatch and extensive checks:
\itemize{
\item \code{qDF}, \code{qDT} and \code{qTBL} convert vectors, matrices, higher-dimensional arrays and suitable lists to data frame, \emph{data.table} and \emph{tibble}, respectively.
\item \code{qM} converts vectors, higher-dimensional arrays, data frames and suitable lists to matrix.
\item \code{mctl} and \code{mrtl} column- or row-wise convert a matrix to list, data frame or \emph{data.table}. They are used internally by \code{qDF/qDT/qTBL}, \code{\link{dapply}}, \code{\link{BY}}, etc\dots
\item \code{\link{qF}} converts atomic vectors to factor (documented on a separate page).
\item \code{as_numeric_factor} and \code{as_character_factor} convert factors, or all factor columns in a data frame / list, to character or numeric (by converting the levels).
}
}
\usage{
# Converting between matrices, data frames / tables / tibbles

 qDF(X, row.names.col = FALSE, keep.attr = FALSE, class = "data.frame")
 qDT(X, row.names.col = FALSE, keep.attr = FALSE, class = c("data.table", "data.frame"))
qTBL(X, row.names.col = FALSE, keep.attr = FALSE, class = c("tbl_df","tbl","data.frame"))
  qM(X, row.names.col = NULL , keep.attr = FALSE, class = NULL, sep = ".")

# Programmer functions: matrix rows or columns to list / DF / DT - fully in C++

mctl(X, names = FALSE, return = "list")
mrtl(X, names = FALSE, return = "list")

# Converting factors or factor columns

  as_numeric_factor(X, keep.attr = TRUE)
as_character_factor(X, keep.attr = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, factor, matrix, higher-dimensional array, data frame or list. \code{mctl} and \code{mrtl} only accept matrices, \code{as_numeric_factor} and \code{as_character_factor} only accept factors, data frames or lists.}
  \item{row.names.col}{can be used to add an column saving names or row.names when converting objects to data frame using \code{qDF/qDT/qTBL}. \code{TRUE} will add a column \code{"row.names"}, or you can supply a name e.g. \code{row.names.col = "variable"}. With \code{qM}, the argument has the opposite meaning, and can be used to select one or more columns in a data frame/list which will be used to create the rownames of the matrix e.g. \code{qM(iris, row.names.col = "Species")}. In this case the column(s) can be specified using names, indices, a logical vector or a selector function. See Examples.}
  \item{keep.attr}{logical. \code{FALSE} (default) yields a \emph{hard} / \emph{thorough} object conversion: All unnecessary attributes are removed from the object yielding a plain matrix / data.frame / \emph{data.table}. \code{FALSE} yields a \emph{soft} / \emph{minimal} object conversion: Only the attributes 'names', 'row.names', 'dim', 'dimnames' and 'levels' are modified in the conversion. Other attributes are preserved. See also \code{class}.}
  \item{class}{if a vector of classes is passed here, the converted object will be assigned these classes. If \code{NULL} is passed, the default classes are assigned: \code{qM} assigns no class, \code{qDF} a class \code{"data.frame"}, and \code{qDT} a class \code{c("data.table", "data.frame")}. If \code{keep.attr = TRUE} and \code{class = NULL} and the object already inherits the default classes, further inherited classes are preserved. See Details and the Example. }
  \item{sep}{character. Separator used for interacting multiple variables selected through \code{row.names.col}.}


  \item{names}{logical. Should the list be named using row/column names from the matrix?}
    \item{return}{an integer or string specifying what to return. The options are:
      \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "list"   \tab\tab returns a plain list \cr
                 2 \tab\tab "data.frame" \tab\tab returns a plain data.frame \cr
                 3 \tab\tab "data.table"   \tab\tab returns a plain \emph{data.table} \cr
  }
  }
}
\details{
Object conversions using these functions are maximally efficient and involve 3 consecutive steps: (1) Converting the storage mode / dimensions / data of the object, (2) converting / modifying the attributes and (3) modifying the class of the object:

(1) is determined by the choice of function and the optional \code{row.names.col} argument. Higher-dimensional arrays are converted by expanding the second dimension (adding columns, same as \code{as.matrix, as.data.frame, as.data.table}).

(2) is determined by the \code{keep.attr} argument: \code{keep.attr = TRUE} seeks to preserve the attributes of the object. Its effect is like copying \code{attributes(converted) <- attributes(original)}, and then modifying the \code{"dim", "dimnames", "names", "row.names"} and \code{"levels"} attributes as necessitated by the conversion task. \code{keep.attr = FALSE} only converts / assigns / removes these attributes and drops all others.

(3) is determined by the \code{class} argument: Setting \code{class = "myclass"} will yield a converted object of class \code{"myclass"}, with any other / prior classes being removed by this replacement. Setting \code{class = NULL} does NOT mean that a class \code{NULL} is assigned (which would remove the class attribute), but rather that the default classes are assigned: \code{qM} assigns no class, \code{qDF} a class \code{"data.frame"}, and \code{qDT} a class \code{c("data.table", "data.frame")}. At this point there is an interaction with \code{keep.attr}: If \code{keep.attr = TRUE} and \code{class = NULL} and the object converted already inherits the respective default classes, then any other inherited classes will also be preserved (with \code{qM(x, keep.attr = TRUE, class = NULL)} any class will be preserved if \code{is.matrix(x)} evaluates to \code{TRUE}.)

The default \code{keep.attr = FALSE} ensures \emph{hard} conversions so that all unnecessary attributes are dropped. Furthermore in \code{qDF/qDT/qTBL} the default classes were explicitly assigned. This is to ensure that the default methods apply, even if the user chooses to preserve further attributes. For \code{qM} a more lenient default setup was chosen to enable the full preservation of time series matrices with \code{keep.attr = TRUE}. If the user wants to keep attributes attached to a matrix but make sure that all default methods work properly, either one of \code{qM(x, keep.attr = TRUE, class = "matrix")} or \code{unclass(qM(x, keep.attr = TRUE))} should be employed.

}
\value{
\code{qDF} - returns a data.frame\cr
\code{qDT} - returns a \emph{data.table}\cr
\code{qTBL} - returns a \emph{tibble}\cr
\code{qM} - returns a matrix\cr
\code{mctl}, \code{mrtl} - return a list, data frame or \emph{data.table} \cr
\code{qF} - returns a factor\cr
\code{as_numeric_factor} - returns X with factors converted to numeric variables\cr
\code{as_character_factor} - returns X with factors converted to character variables
}
% \note{
% \code{qTBL} works similarly to \code{qDT} assigning different classes, i.e. \code{qTBL(x)} is equivalent to \code{qDT(x, class = c("tbl_df", "tbl", "data.frame"))}. Similar converters for other data frame based classes are easily created from \code{qDF} and \code{qDT}. The principle difference between them is that \code{qDF} preserves rownames whereas \code{qDT} always assigns integer rownames.
% }
\seealso{
\code{\link{qF}}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Basic Examples
mtcarsM <- qM(mtcars)                   # Matrix from data.frame
mtcarsDT <- qDT(mtcarsM)                # data.table from matrix columns
mtcarsTBL <- qTBL(mtcarsM)              # tibble from matrix columns
head(mrtl(mtcarsM, TRUE, "data.frame")) # data.frame from matrix rows, etc..
head(qDF(mtcarsM, "cars"))              # Adding a row.names column when converting from matrix
head(qDT(mtcars, "cars"))               # Saving row.names when converting data frame to data.table
head(qM(iris, "Species"))               # Examples converting data to matrix, saving information
head(qM(GGDC10S, is.character))         # as rownames
head(qM(gv(GGDC10S, -(2:3)), 1:3, sep = "-")) # plm-style rownames

# mrtl() and mctl() are very useful for iteration over matrices
# Think of a coordninates matrix e.g. from sf::st_coordinates()
coord <- matrix(rnorm(10), ncol = 2, dimnames = list(NULL, c("X", "Y")))
# Then we can
for (d in mrtl(coord)) {
  cat("lon =", d[1], ", lat =", d[2], fill = TRUE)
  # do something complicated ...
}
rm(coord)

## Factors
cylF <- qF(mtcars$cyl)                  # Factor from atomic vector
cylF

# Factor to numeric conversions
identical(mtcars,  as_numeric_factor(dapply(mtcars, qF)))
% ## Explaining the interaction of keep.attr and class. Consider the time series   EuStockMarkets
% plot()
}
\keyword{manip}
\keyword{documentation}
