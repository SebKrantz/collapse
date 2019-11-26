\name{select-replace-vars}
\alias{select-replace-vars}
\alias{get_vars}
\alias{num_vars}
\alias{cat_vars}
\alias{char_vars}
\alias{fact_vars}
\alias{logi_vars}
\alias{Date_vars}
\alias{get_vars<-}
\alias{num_vars<-}
\alias{cat_vars<-}
\alias{char_vars<-}
\alias{fact_vars<-}
\alias{logi_vars<-}
\alias{Date_vars<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quick Select and Replace Data Frame Columns} % by Type, Function or Regular Expression}
\description{
Efficiently select or replace a subset of columns from a data frame. This can be done by data type, or using column names, indices, logical vectors, functions or regular expressions.

The performance is generally faster and more secure than \code{`[`} (secure w.r.t. redefinitions of \code{`[.data.frame`} or \code{`[<-.data.frame`} for other classes and preventing the loss of attributes, but not in the sense of performing all kinds of costly checks in the presence of malformed data.frame's).
}
\usage{
## Select and replace columns by data type
num_vars(x, return = c("data","names","indices","named_indices"))
num_vars(x) <- value
cat_vars(x, return = c("data","names","indices","named_indices"))
cat_vars(x) <- value
char_vars(x, return = c("data","names","indices","named_indices"))
char_vars(x) <- value
fact_vars(x, return = c("data","names","indices","named_indices"))
fact_vars(x) <- value
logi_vars(x, return = c("data","names","indices","named_indices"))
logi_vars(x) <- value
Date_vars(x, return = c("data","names","indices","named_indices"))
Date_vars(x) <- value

## Select and replace columns by names, indices, logical vectors,
## regular expressions or using other functions to identify columns
get_vars(x, vars, return = c("data","names","indices","named_indices"),
         regex = FALSE, ...)
get_vars(x, vars, regex = FALSE, ...) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame. \emph{Note}: Functions als work with data.table's etc... as long as \code{x} is a list and \code{!is.null(dim(x))}.}
  \item{value}{a data.frame or list of columns whose dimensions exactly match those of the extracted subet of \code{x}.}
  \item{vars}{column names (can be negative), indices, a suitable logical vector, a vector of regular expressions matching column names if \code{regex = TRUE}, or a function returning \code{TRUE} or \code{FALSE} when applied to the columns of \code{x}.}
  \item{return}{an integer or string specifying what to return. The options are:
  \enumerate{
  \item "data"\cr
  \item "names"\cr
  \item "indices"\cr
  \item "named_indices"
  }
  \emph{Note}: replacement functions currently only replace data, not column names or indices (ordering). However column names are replaced together with the data.
  }
\item{regex}{logical. \code{TRUE} will do regular expression search on the column names of \code{x} using a (vector of) regular expression(s) passed to \code{vars}.}
\item{...}{further arguments passed to \code{\link{grep}}, if \code{regex = TRUE}.}
}
\details{
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Wold Development Data
head(num_vars(wlddev))                                     # Select numeric variables
head(get_vars(wlddev, is.numeric))                         # Same thing
head(cat_vars(wlddev))                                     # Select categorical (non-numeric) variable
head(get_vars(wlddev, is.categorical))                     # Same thing

num_vars(wlddev) <- num_vars(wlddev)                       # Replace Numeric Variables by themselves
get_vars(wlddev,is.numeric) <- get_vars(wlddev,is.numeric) # Same thing

head(get_vars(wlddev, 9:12))                               # Select columns 9 through 12. 2x faster than wlddev[9:12]
head(get_vars(wlddev, -(9:12)))                            # All except columns 9 through 12
head(get_vars(wlddev, c("PCGDP","LIFEEX","GINI","ODA")))   # Select using column names
head(get_vars(wlddev, "[[:upper:]]", regex = TRUE))        # Same thing: Selecting variable names with only upper-case letters

get_vars(wlddev, 9:12) <- get_vars(wlddev, 9:12)           # Again simple replacement, 4-5x faster than wlddev[9:12] <- wlddev[9:12]
get_vars(wlddev, 9:12) <-                                  # Standardizing columns 9 through 12
                STD(get_vars(wlddev,9:12), wlddev$iso3c)   # Note that this also change the names in the data, whereas
head(wlddev)                                               # wlddev[9:12] <- ... replaces only data, not names!
rm(wlddev)
}
