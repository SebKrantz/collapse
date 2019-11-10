\name{select_variables}
\alias{select_variables}
\alias{num.vars}
\alias{cat.vars}
\alias{char.vars}
\alias{fact.vars}
\alias{logi.vars}
\alias{get.vars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Select and Replace Data Frame Columns by Type, Function or Regular Expression}
\description{
Quickly return a subset of a data.frame, or alternatively of its variable names or indices, matching either certain data types, an identifier function applied to the columns, or a regular expression applied to the variable names.
}
\usage{
## Select columns by type
num.vars(x, return = "data")
cat.vars(x, return = "data")
char.vars(x, return = "data")
fact.vars(x, return = "data")
logi.vars(x, return = "data")

## Select columns by function or regular expression
get.vars(x, FoR, return = "data", ...)

## Replace columns by type
num.vars(x) <- y
cat.vars(x) <- y
char.vars(x) <- y
fact.vars(x) <- y
logi.vars(x) <- y

## Replace columns by function or regular expression
get.vars(x, FoR, ...) <- y
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data.frame.}
  \item{y}{a data.frame or list of columns whose dimensions exactly match those of the extracted subet of \code{x}.}
  \item{FoR}{a function returning TRUE or FALSE when applied to the columns of \code{x}, or a character vector of regular expressions.}
  \item{return}{an integer or string specifying what to return. The options are:
  \enumerate{
  \item "data"\cr
  \item "names"\cr
  \item "indices"\cr
  \item "named.indices"
  }
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
