\name{fselect-get_vars-add_vars}
% \name{select-replace-vars}
% \alias{select-replace-vars}
\alias{fselect}
\alias{fselect<-}
\alias{slt}
\alias{slt<-}
\alias{get_vars}
\alias{gv}
\alias{gvr}
\alias{num_vars}
\alias{nv}
\alias{cat_vars}
\alias{char_vars}
\alias{fact_vars}
\alias{logi_vars}
\alias{date_vars}
\alias{add_vars}
\alias{av}
\alias{get_vars<-}
\alias{gv<-}
\alias{gvr<-}
\alias{num_vars<-}
\alias{nv<-}
\alias{cat_vars<-}
\alias{char_vars<-}
\alias{fact_vars<-}
\alias{logi_vars<-}
\alias{date_vars<-}
\alias{add_vars<-}
\alias{av<-}
\title{Fast Select, Replace or Add Data Frame Columns}
\description{
Efficiently select and replace (or add) a subset of columns from (to) a data frame. This can be done by data type, or using expressions, column names, indices, logical vectors, selector functions or regular expressions matching column names.
}
\usage{
## Select and replace variables, analgous to dplyr::select but significantly faster
fselect(.x, \dots, return = "data")
fselect(x, \dots) <- value
slt(.x, \dots, return = "data")   # Shorthand for fselect
slt(x, \dots) <- value            # Shorthand for fselect<-

## Select and replace columns by names, indices, logical vectors,
## regular expressions or using functions to identify columns

get_vars(x, vars, return = "data", regex = FALSE, rename = FALSE, \dots)
      gv(x, vars, return = "data", \dots)   # Shorthand for get_vars
     gvr(x, vars, return = "data", \dots)   # Shorthand for get_vars(..., regex = TRUE)

get_vars(x, vars, regex = FALSE, \dots) <- value
      gv(x, vars, \dots) <- value           # Shorthand for get_vars<-
     gvr(x, vars, \dots) <- value           # Shorthand for get_vars<-(..., regex = TRUE)

## Add columns at any position within a data.frame

add_vars(x, \dots, pos = "end")
add_vars(x, pos = "end") <- value
      av(x, \dots, pos = "end")             # Shorthand for add_vars
      av(x, pos = "end") <- value         # Shorthand for add_vars<-

## Select and replace columns by data type

num_vars(x, return = "data")
num_vars(x) <- value
      nv(x, return = "data")       # Shorthand for num_vars
      nv(x) <- value               # Shorthand for num_vars<-
cat_vars(x, return = "data")       # Categorical variables, see is_categorical
cat_vars(x) <- value
char_vars(x, return = "data")
char_vars(x) <- value
fact_vars(x, return = "data")
fact_vars(x) <- value
logi_vars(x, return = "data")
logi_vars(x) <- value
date_vars(x, return = "data")      # See is_date
date_vars(x) <- value
}

\arguments{
  \item{x, .x}{a data frame or list.}
  \item{value}{a data frame or list of columns whose dimensions exactly match those of the extracted subset of \code{x}. If only 1 variable is in the subset of \code{x}, \code{value} can also be an atomic vector or matrix, provided that \code{NROW(value) == nrow(x)}.}
  \item{vars}{a vector of column names, indices (can be negative), a suitable logical vector, or a vector of regular expressions matching column names (if \code{regex = TRUE}). It is also possible to pass a function returning \code{TRUE} or \code{FALSE} when applied to the columns of \code{x}.}
  \item{return}{an integer or string specifying what the selector function should return. The options are:
      \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "data"   \tab\tab subset of data frame (default) \cr
                 2 \tab\tab "names" \tab\tab column names \cr
                 3 \tab\tab "indices"   \tab\tab column indices \cr
                 4 \tab\tab "named_indices"  \tab\tab named column indices \cr
                 5 \tab\tab "logical"   \tab\tab logical selection vector \cr
                 6 \tab\tab "named_logical"     \tab\tab named logical vector \cr
  }
  \emph{Note}: replacement functions only replace data, however column names are replaced together with the data (if available).
  }
\item{regex}{logical. \code{TRUE} will do regular expression search on the column names of \code{x} using a (vector of) regular expression(s) passed to \code{vars}. Matching is done using \code{\link{grep}}.}
\item{rename}{logical. If \code{vars} is a named vector of column names or indices, \code{rename = TRUE} will use the (non missing) names to rename columns.}
\item{pos}{the position where columns are added in the data frame. \code{"end"} (default) will append the data frame at the end (right) side. "front" will add columns in front (left). Alternatively one can pass a vector of positions (matching \code{length(value)} if value is a list). In that case the other columns will be shifted around the new ones while maintaining their order. }
\item{\dots}{for \code{fselect}: column names and expressions e.g. \code{fselect(mtcars, newname = mpg, hp, carb:vs)}. for \code{get_vars}: further arguments passed to \code{\link{grep}}, if \code{regex = TRUE}. For \code{add_vars}: multiple lists/data frames or vectors (which should be given names e.g. \code{name = vector}). A single argument passed may also be an (unnamed) vector or matrix.}
}

\details{
\code{get_vars(<-)} is around 2x faster than \code{`[.data.frame`} and 8x faster than \code{`[<-.data.frame`}, so the common operation \code{data[cols] <- someFUN(data[cols])} can be made 10x more efficient (abstracting from computations performed by \code{someFUN}) using \code{get_vars(data, cols) <- someFUN(get_vars(data, cols))} or the shorthand \code{gv(data, cols) <- someFUN(gv(data, cols))}.

Similarly type-wise operations like \code{data[sapply(data, is.numeric)]} or \code{data[sapply(data, is.numeric)] <- value} are facilitated and more efficient using \code{num_vars(data)} and \code{num_vars(data) <- value} or the shortcuts \code{nv} and \code{nv<-} etc.

\code{fselect} provides an efficient alternative to \code{dplyr::select}, allowing the selection of variables based on expressions evaluated within the data frame, see Examples. It is about 100x faster than \code{dplyr::select} but also more simple as it does not provide special methods (except for 'sf' and 'data.table' which are handled internally) .

Finally, \code{add_vars(data1, data2, data3, \dots)} is a lot faster than \code{cbind(data1, data2, data3, \dots)}, and preserves the attributes of \code{data1} (i.e. it is like adding columns to \code{data1}). The replacement function \code{add_vars(data) <- someFUN(get_vars(data, cols))} efficiently appends \code{data} with computed columns. The \code{pos} argument allows adding columns at positions other than the end (right) of the data frame, see Examples. \emph{Note} that \code{add_vars} does not check duplicated column names or \code{NULL} columns, and does not evaluate expressions in a data environment, or replicate length 1 inputs like \code{\link{cbind}}. All of this is provided by \code{\link{ftransform}}.

All functions introduced here perform their operations class-independent. They all basically work like this: (1) save the attributes of \code{x}, (2) unclass \code{x}, (3) subset, replace or append \code{x} as a list, (4) modify the "names" component of the attributes of \code{x} accordingly and (5) efficiently attach the attributes again to the result from step (3).
Thus they can freely be applied to data.table's, grouped tibbles, panel data frames and other classes and will return an object of exactly the same class and the same attributes.

% secure w.r.t. redefinitions of \code{`[.data.frame`} or \code{`[<-.data.frame`} for other classes (i.e. data.table's, tibbles etc.) and preserve all attributes of the data
}
\note{
In many cases functions here only check the length of the first column, which is one of the reasons why they are so fast. When lists of unequal-length columns are offered as replacements this yields a malformed data frame (which will also print a warning in the console i.e. you will notice that).
}
\seealso{
\code{\link{fsubset}}, \code{\link{ftransform}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Wold Development Data
head(fselect(wlddev, Country = country, Year = year, ODA)) # Fast dplyr-like selecting
head(fselect(wlddev, -country, -year, -PCGDP))
head(fselect(wlddev, country, year, PCGDP:ODA))
head(fselect(wlddev, -(PCGDP:ODA)))
fselect(wlddev, country, year, PCGDP:ODA) <- NULL          # Efficient deleting
head(wlddev)
rm(wlddev)

head(num_vars(wlddev))                                     # Select numeric variables
head(cat_vars(wlddev))                                     # Select categorical (non-numeric) vars
head(get_vars(wlddev, is_categorical))                     # Same thing

num_vars(wlddev) <- num_vars(wlddev)                       # Replace Numeric Variables by themselves
get_vars(wlddev,is.numeric) <- get_vars(wlddev,is.numeric) # Same thing

head(get_vars(wlddev, 9:12))                               # Select columns 9 through 12, 2x faster
head(get_vars(wlddev, -(9:12)))                            # All except columns 9 through 12
head(get_vars(wlddev, c("PCGDP","LIFEEX","GINI","ODA")))   # Select using column names
head(get_vars(wlddev, "[[:upper:]]", regex = TRUE))        # Same thing: match upper-case var. names
head(gvr(wlddev, "[[:upper:]]"))                           # Same thing

get_vars(wlddev, 9:12) <- get_vars(wlddev, 9:12)           # 9x faster wlddev[9:12] <- wlddev[9:12]
add_vars(wlddev) <- STD(gv(wlddev,9:12), wlddev$iso3c)     # Add Standardized columns 9 through 12
head(wlddev)                                               # gv and av are shortcuts

get_vars(wlddev, 14:17) <- NULL                            # Efficient Deleting added columns again
av(wlddev, "front") <- STD(gv(wlddev,9:12), wlddev$iso3c)  # Again adding in Front
head(wlddev)
get_vars(wlddev, 1:4) <- NULL                              # Deleting
av(wlddev,c(10,12,14,16)) <- W(wlddev,~iso3c, cols = 9:12, # Adding next to original variables
                               keep.by = FALSE)
head(wlddev)
get_vars(wlddev, c(10,12,14,16)) <- NULL                   # Deleting

head(add_vars(wlddev, new = STD(wlddev$PCGDP)))                  # Can also add columns like this
head(add_vars(wlddev, STD(nv(wlddev)), new = W(wlddev$PCGDP)))   # etc...

head(add_vars(mtcars, mtcars, mpg = mtcars$mpg, mtcars), 2)      # add_vars does not check names!
}
\keyword{manip}
