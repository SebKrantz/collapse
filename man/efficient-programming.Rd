\name{efficient-programming}
\alias{AA2-efficient-programming}
\alias{efficient-programming}
\alias{anyv}
\alias{allv}
\alias{allNA}
\alias{whichv}
\alias{whichNA}
\alias{alloc}
\alias{copyv}
\alias{setv}
\alias{setop}
\alias{\%==\%}
\alias{\%!=\%}
\alias{\%+=\%}
\alias{\%-=\%}
\alias{\%*=\%}
\alias{\%/=\%}
\alias{cinv}
\alias{vlengths}
\alias{vtypes}
\alias{vgcd}
\alias{fnlevels}
\alias{fnrow}
\alias{fncol}
\alias{fdim}
\alias{missing_cases}
\alias{na_rm}
\alias{na_omit}
\alias{na_insert}
\alias{seq_row}
\alias{seq_col}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Small Functions to Make R Programming More Efficient
}
\description{
 A small set of functions to addresses some common inefficiencies in R, such as the creation of logical vectors to compare quantities, unnecessary copies of objects in elementary mathematical or subsetting operations, obtaining information about objects (esp. data frames), or dealing with missing values. % It makes sense to use them when dealing with > 10,000 obs. on modern computers.
}
\usage{
anyv(x, value)              # Faster than any(x == value). See also kit::panyv()
allv(x, value)              # Faster than all(x == value). See also kit::pallv()
allNA(x)                    # Faster than all(is.na(x)). See also kit::pallNA()
whichv(x, value,            # Faster than which(x == value)
       invert = FALSE)      # or which(x != value). See also Note (3)
whichNA(x, invert = FALSE)  # Faster than which((!)is.na(x))
x \%==\% value                # Infix for whichv(v, value, FALSE), use e.g. in fsubset()
x \%!=\% value                # Infix for whichv(v, value, TRUE). See also Note (3)
alloc(value, n)             # Fast rep_len(value, n) or replicate(n, value). See Details
copyv(X, v, R, \dots, invert  # Fast replace(X, v, R), replace(X, X (!/=)= v, R) or
    = FALSE, vind1 = FALSE, # replace(X, (!)v, R[(!)v]). See Details and Note (4).
    xlist = FALSE)          # For multi-replacement see also kit::vswitch()
setv(X, v, R, \dots, invert   # Same for X[v] <- r, X[x (!/=)= v] <- r or
    = FALSE, vind1 = FALSE, # x[(!)v] <- r[(!)v]. Modifies X by reference, fastest.
    xlist = FALSE)          # X/R/V can also be lists/DFs. See Details and Examples.
setop(X, op, V, \dots,        # Faster than X <- X +\-\*\/ V (modifies by reference)
      rowwise = FALSE)      # optionally can also add v to rows of a matrix or list
X \%+=\% V                    # Infix for setop(X, "+", V). See also Note (2)
X \%-=\% V                    # Infix for setop(X, "-", V). See also Note (2)
X \%*=\% V                    # Infix for setop(X, "*", V). See also Note (2)
X \%/=\% V                    # Infix for setop(X, "/", V). See also Note (2)
na_rm(x)                    # Fast: if(anyNA(x)) x[!is.na(x)] else x,
                            # also removes NULL / empty elements from list
na_omit(X, cols = NULL,     # Faster na.omit for matrices and data frames,
        na.attr = FALSE,    # can use selected columns to check, attach indices,
        prop = 0, ...)      # and remove cases with a proportion of values missing
na_insert(X, prop = 0.1,    # Insert missing values at random
          value = NA)
missing_cases(X, cols=NULL, # The oposite of complete.cases(), faster for DF's.
  prop = 0, count = FALSE)  # See also kit::panyNA(), kit::pallNA(), kit::pcountNA()
vlengths(X, use.names=TRUE) # Faster lengths() and nchar() (in C, no method dispatch)
vtypes(X, use.names = TRUE) # Get data storage types (faster vapply(X, typeof, ...))
vgcd(x)                     # Greatest common divisor of positive integers or doubles
fnlevels(x)                 # Faster version of nlevels(x) (for factors)
fnrow(X)                    # Faster nrow for data frames (not faster for matrices)
fncol(X)                    # Faster ncol for data frames (not faster for matrices)
fdim(X)                     # Faster dim for data frames (not faster for matrices)
seq_row(X)                  # Fast integer sequences along rows of X
seq_col(X)                  # Fast integer sequences along columns of X
cinv(x)                     # Choleski (fast) inverse of symmetric PD matrix, e.g. X'X
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X, V, R}{a vector, matrix or data frame.}
  \item{x, v}{a (atomic) vector or matrix (\code{na_rm} also supports lists).}
  \item{value}{a single value of any (atomic) vector type. For \code{whichv} it can also be a \code{length(x)} vector.}
  \item{invert}{logical. \code{TRUE} considers elements \code{x != value}.}
  \item{vind1}{logical. If \code{length(v) == 1L}, setting \code{vind1 = TRUE} will interpret \code{v} as an index, rather than a value to search and replace.}
  \item{xlist}{logical. If \code{X} is a list, the default is to treat it like a data frame and replace rows. Setting \code{xlist = TRUE} will treat \code{X} and its replacement \code{R} like 1-dimensional list vectors.}
  \item{op}{an integer or character string indicating the operation to perform.
    \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab \code{"+"}   \tab\tab add \code{V} \cr
                 2 \tab\tab \code{"-"} \tab\tab subtract \code{V} \cr
                 3 \tab\tab \code{"*"}   \tab\tab multiply by \code{V} \cr
                 4 \tab\tab \code{"/"}  \tab\tab divide by \code{V} \cr
    }
  }
  \item{rowwise}{logical. \code{TRUE} performs the operation between \code{V} and each row of \code{X}.}
   \item{cols}{select columns to check for missing values using column names, indices, a logical vector or a function (e.g. \code{is.numeric}). The default is to check all columns, which could be inefficient.}
  \item{n}{integer. The length of the vector to allocate with \code{value}.}
  \item{na.attr}{logical. \code{TRUE} adds an attribute containing the removed cases. For compatibility reasons this is exactly the same format as \code{na.omit} i.e. the attribute is called "na.action" and of class "omit".}
  \item{prop}{double. For \code{na_insert}: the proportion of observations to be randomly replaced with \code{NA}. For \code{missing_cases} and \code{na_omit}: the proportion of values missing for the case to be considered missing (within \code{cols} if specified). For matrices this is implemented in R as \code{rowSums(is.na(X)) >= max(as.integer(prop * ncol(X)), 1L)}. The C code for data frames works equivalently, and skips list- and raw-columns (\code{ncol(X)} is adjusted downwards).}
  \item{count}{logical. \code{TRUE} returns the row-wise missing value count (within \code{cols}). This ignores \code{prop}.}
  \item{use.names}{logical. Preserve names if \code{X} is a list. }
  \item{\dots}{for \code{na_omit}: further arguments passed to \code{[} for vectors and matrices. With indexed data it is also possible to specify the \code{drop.index.levels} argument, see \link{indexing}. For \code{copyv}, \code{setv} and \code{setop}, the argument is unused, and serves as a placeholder for possible future arguments.}
}
\details{
\code{alloc} is a fusion of \code{\link{rep_len}} and \code{\link{replicate}} that is faster in both cases. If \code{value} is a length one atomic vector (logical, integer, double, string, complex or raw), the functionality is as \code{rep_len(value, n)} i.e. the output is a length \code{n} atomic vector with the same attributes as \code{value} (apart from \code{"names"}, \code{"dim"} and \code{"dimnames"}). For all other cases the functionality is as \code{replicate(n, value, simplify = FALSE)} i.e. the output is a length-\code{n} list of the objects. If \code{value} is a list itself, the list is not deep-copied on each replication.

\code{copyv} and \code{setv} are designed to optimize operations that require replacing data in objects in the broadest sense. The only difference between them is that \code{copyv} first deep-copies \code{X} before doing replacements whereas \code{setv} modifies \code{X} in place and returns the result invisibly. There are 3 ways these functions can be used:

\enumerate{
\item To replace a single value, \code{setv(X, v, R)} is an efficient alternative to \code{X[X == v] <- R}, and  \code{copyv(X, v, R)} is more efficient than \code{replace(X, X == v, R)}. This can be inverted using \code{setv(X, v, R, invert = TRUE)}, equivalent to \code{X[X != v] <- R}.

\item To do standard replacement with integer or logical indices i.e. \code{X[v] <- R} is more efficient using \code{setv(X, v, R)}, and, if \code{v} is logical, \code{setv(X, v, R, invert = TRUE)} is efficient for \code{X[!v] <- R}. To distinguish this from use case (1) when \code{length(v) == 1}, the argument \code{vind1 = TRUE} can be set to ensure that \code{v} is always interpreted as an index.

\item To copy values from objects of equal size i.e. \code{setv(X, v, R)} is faster than \code{X[v] <- R[v]}, and \code{setv(X, v, R, invert = TRUE)} is faster than \code{X[!v] <- R[!v]}.
}

Both \code{X} and \code{R} can be atomic or data frames / lists. If \code{X} is a list, the default behavior is to interpret it like a data frame, and apply \code{setv/copyv} to each element/column of \code{X}. If \code{R} is also a list, this is done using \code{\link{mapply}}. Thus \code{setv/copyv} can also be used to replace elements or rows in data frames, or copy rows from equally sized frames. Note that for replacing subsets in data frames \code{\link[data.table]{set}} from \code{data.table} provides a more convenient interface (and there is also \code{\link[data.table]{copy}} if you just want to deep-copy an object without any modifications to it).

If \code{X} should not be interpreted like a data frame, setting \code{xlist = TRUE} will interpret it like a 1D list-vector analogous to atomic vectors, except that use case (1) is not permitted i.e. no value comparisons on list elements.
%\code{copyv} and \code{setv} perform different tasks, depending on the input. If \code{v} is a scalar, the elements of \code{X} are compared to \code{v}, and the matching ones (or non-matching ones if \code{invert = TRUE}) are replaced with \code{R}, where \code{R} can be either a scalar or an object of the same dimensions as \code{X}. If \code{X} is a data frame, \code{R} can also be a column-vector matching \code{fnrow(X)}. The second option is if \code{v} is either a logical or integer vector of indices with \code{length(v) > 1L}, indicating the elements of a vector / matrix (or rows if \code{X} is a data frame) to replace with corresponding elements from \code{R}. Thus \code{R} has to be of equal dimensions as \code{X}, but could also be a column-vector if \code{X} is a data frame. Setting \code{vind1 = TRUE} ensures that \code{v} is always interpreted as an index, even if \code{length(v) == 1L}.

% In this case \code{r} has to be a vector of the same length as \code{x}, and the corresponding elements in \code{v} are replaced with their counterparts in \code{r}. \code{copyv} does all that by first creating a copy of \code{x}, whereas \code{setv} modifies \code{x} directly and is thus more efficient.
}

\note{

\enumerate{
\item None of these functions (apart from \code{alloc}) currently support complex vectors.

\item \code{setop} and the operators \code{\%+=\%}, \code{\%-=\%}, \code{\%*=\%} and \code{\%/=\%} also work with integer data, but do not perform any integer related checks. R's integers are bounded between +-2,147,483,647 and \code{NA_integer_} is stored as the value -2,147,483,648. Thus computations resulting in values exceeding +-2,147,483,647 will result in integer overflows, and \code{NA_integer_} should not occur on either side of a \code{setop} call. These are programmers functions and meant to provide the most efficient math possible to responsible users.

\item It is possible to compare factors by the levels (e.g. \code{iris$Species \%==\% "setosa")}) or using integers (\code{iris$Species \%==\% 1L}). The latter is slightly more efficient. Nothing special is implemented for other objects apart from basic types, e.g. for dates (which are stored as doubles) you need to generate a date object i.e. \code{wlddev$date \%==\% as.Date("2019-01-01")}. Using \code{wlddev$date \%==\% "2019-01-01"} will give \code{integer(0)}.

\item \code{setv/copyv} only allow positive integer indices being passed to \code{v}, and, for efficiency reasons, they only check the first and the last index. Thus if there are indices in the middle that fall outside of the data range it will terminate R.

}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[=data-transformations]{Data Transformations}, \link[=small-helpers]{Small (Helper) Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{

## Which value
whichNA(wlddev$PCGDP)                # Same as which(is.na(wlddev$PCGDP))
whichNA(wlddev$PCGDP, invert = TRUE) # Same as which(!is.na(wlddev$PCGDP))
whichv(wlddev$country, "Chad")       # Same as which(wlddev$county == "Chad")
wlddev$country \%==\% "Chad"           # Same thing
whichv(wlddev$country, "Chad", TRUE) # Same as which(wlddev$county != "Chad")
wlddev$country \%!=\% "Chad"           # Same thing
lvec <- wlddev$country == "Chad"     # If we already have a logical vector...
whichv(lvec, FALSE)                  # is fastver than which(!lvec)
rm(lvec)

# Using the \%==\% operator can yield tangible performance gains
fsubset(wlddev, iso3c \%==\% "DEU") # 3x faster than:
fsubset(wlddev, iso3c == "DEU")

## Math by reference: permissible types of operations
x <- alloc(1.0, 1e5) # Vector
x \%+=\% 1
x \%+=\% 1:1e5
xm <- matrix(alloc(1.0, 1e5), ncol = 100) # Matrix
xm \%+=\% 1
xm \%+=\% 1:1e3
setop(xm, "+", 1:100, rowwise = TRUE)
xm \%+=\% xm
xm \%+=\% 1:1e5
xd <- qDF(replicate(100, alloc(1.0, 1e3), simplify = FALSE)) # Data Frame
xd \%+=\% 1
xd \%+=\% 1:1e3
setop(xd, "+", 1:100, rowwise = TRUE)
xd \%+=\% xd
rm(x, xm, xd)

## setv() and copyv()
x <- rnorm(100)
y <- sample.int(10, 100, replace = TRUE)
setv(y, 5, 0)            # Faster than y[y == 5] <- 0
setv(y, 4, x)            # Faster than y[y == 4] <- x[y == 4]
setv(y, 20:30, y[40:50]) # Faster than y[20:30] <- y[40:50]
setv(y, 20:30, x)        # Faster than y[20:30] <- x[20:30]
rm(x, y)

# Working with data frames, here returning copies of the frame
copyv(mtcars, 20:30, ss(mtcars, 10:20))
copyv(mtcars, 20:30, fscale(mtcars))
ftransform(mtcars, new = copyv(cyl, 4, vs))
# Column-wise:
copyv(mtcars, 2:3, fscale(mtcars), xlist = TRUE)
copyv(mtcars, 2:3, mtcars[4:5], xlist = TRUE)

## Missing values
mtc_na <- na_insert(mtcars, 0.15)    # Set 15\% of values missing at random
fnobs(mtc_na)                        # See observation count
missing_cases(mtc_na)                # Fast equivalent to !complete.cases(mtc_na)
missing_cases(mtc_na, cols = 3:4)    # Missing cases on certain columns?
missing_cases(mtc_na, count = TRUE)  # Missing case count
missing_cases(mtc_na, prop = 0.8)    # Cases with 80\% or more missing
missing_cases(mtc_na, cols = 3:4, prop = 1)     # Cases mssing columns 3 and 4
missing_cases(mtc_na, cols = 3:4, count = TRUE) # Missing case count on columns 3 and 4

na_omit(mtc_na)                      # 12x faster than na.omit(mtc_na)
na_omit(mtc_na, prop = 0.8)          # Only remove cases missing 80\% or more
na_omit(mtc_na, na.attr = TRUE)      # Adds attribute with removed cases, like na.omit
na_omit(mtc_na, cols = .c(vs, am))   # Removes only cases missing vs or am
na_omit(qM(mtc_na))                  # Also works for matrices
na_omit(mtc_na$vs, na.attr = TRUE)   # Also works with vectors
na_rm(mtc_na$vs)                     # For vectors na_rm is faster ...
rm(mtc_na)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
\keyword{utilities} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{manip}
\keyword{math}
\keyword{documentation}
