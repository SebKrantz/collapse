\name{efficient-programming}
\alias{AA2-efficient-programming}
\alias{efficient-programming}
\alias{anyv}
\alias{allv}
\alias{allNA}
\alias{whichv}
\alias{whichNA}
\alias{frange}
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
anyv(x, value)              # Faster than any(x == value)
allv(x, value)              # Faster than all(x == value)
allNA(x)                    # Faster than all(is.na(x))
whichv(x, value,            # Faster than which(x == value)
       invert = FALSE)      # or which(x != value). See also Note (3)
whichNA(x, invert = FALSE)  # Faster than which((!)is.na(x))
x \%==\% value                # Infix for whichv(v, value, FALSE), use e.g. in fsubset
x \%!=\% value                # Infix for whichv(v, value, TRUE). See also Note (3)
alloc(value, n)             # Faster than rep_len(value, n)
copyv(X, v, R, \dots, invert  # Faster than replace(x, x == v, r) or replace(x, v, r[v])
  = FALSE, vind1 = FALSE)   # or replace(x, x != v, r) or replace(x, !v, r[!v])
setv(X, v, R, \dots, invert   # Same for x[x (!/=)= v] <- r or x[(!)v] <- r[(!)v]
  = FALSE, vind1 = FALSE)   # modifies x by reference, fastest
setop(X, op, V, \dots,        # Faster than X <- X +\-\*\/ V (modifies by reference)
      rowwise = FALSE)      # optionally can also add v to rows of a matrix or list
X \%+=\% V                    # Infix for setop(X, "+", V). See also Note (2)
X \%-=\% V                    # Infix for setop(X, "-", V). See also Note (2)
X \%*=\% V                    # Infix for setop(X, "*", V). See also Note (2)
X \%/=\% V                    # Infix for setop(X, "/", V). See also Note (2)
na_rm(x)                    # Fast: if(anyNA(x)) x[!is.na(x)] else x,
                            # also removes NULL / empty elements from list
na_omit(X, cols = NULL,     # Faster na.omit for matrices and data frames,
      na.attr = FALSE, ...) # can use selected columns and attach indices
na_insert(X, prop = 0.1,    # Insert missing values at random
          value = NA)
missing_cases(X,            # The oposite of complete.cases(), faster for
              cols = NULL)  # data frames
vlengths(X, use.names=TRUE) # Faster version of lengths() (in C, no method dispatch)
vtypes(X, use.names = TRUE) # Get data storage types (faster vapply(X, typeof, ...))
vgcd(x)                     # Greatest common divisor of positive integers or doubles
frange(x, na.rm = TRUE)     # Much faster base::range, for integer and double objects
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
  \item{vind1}{logical. If \code{length(v) == 1L}, setting \code{vind1 = TRUE} will interpret \code{v} as an index of \code{X} and \code{R}, rather than a value to search and replace.}
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
  \item{prop}{double. Specify the proportion of observations randomly replaced with \code{NA}.}
  \item{use.names}{logical. Preserve names if \code{X} is a list. }
  \item{na.rm}{logical. \code{TRUE} omits missing values by skipping them in the computation. \code{FALSE} terminates the computation once a missing value is encountered and returns 2 missing values.}
  \item{\dots}{for \code{na_omit}: further arguments passed to \code{[} for vectors and matrices. With indexed data it is also possible to specify the \code{drop.index.levels} argument, see \link{indexing}. For \code{copyv}, \code{setv} and \code{setop}, the argument is unused, and serves as a placeholder for possible future arguments.}
}
\details{
\code{copyv} and \code{setv} are designed to optimize operations that require replacing a single value in an object e.g. \code{X[X == value] <- r} or \code{X[X == value] <- R[R == value]} or simply copying parts of an existing object into another object e.g. \code{X[v] <- R[v]}. Thus they only cover cases where base R is inefficient by either creating a logical vector or materializing a subset to do some replacement. No alternative is provided in cases where base R is efficient e.g. \code{x[v] <- r} or cases provided by \code{\link[data.table]{set}} and \code{\link[data.table]{copy}} from the \emph{data.table} package. Both functions work equivalently, with the difference that \code{copyv} creates a deep copy of the data before making the replacements and returns the copy, whereas \code{setv} modifies the data directly without creating a copy and returns the modified object invisibly. Thus \code{setv} is considerably more efficient.

\code{copyv} and \code{setv} perform different tasks, depending on the input. If \code{v} is a scalar, the elements of \code{X} are compared to \code{v}, and the matching ones (or non-matching ones if \code{invert = TRUE}) are replaced with \code{R}, where \code{R} can be either a scalar or an object of the same dimensions as \code{X}. If \code{X} is a data frame, \code{R} can also be a column-vector matching \code{fnrow(X)}. The second option is if \code{v} is either a logical or integer vector of indices with \code{length(v) > 1L}, indicating the elements of a vector / matrix (or rows if \code{X} is a data frame) to replace with corresponding elements from \code{R}. Thus \code{R} has to be of equal dimensions as \code{X}, but could also be a column-vector if \code{X} is a data frame. Setting \code{vind1 = TRUE} ensures that \code{v} is always interpreted as an index, even if \code{length(v) == 1L}.

% In this case \code{r} has to be a vector of the same length as \code{x}, and the corresponding elements in \code{v} are replaced with their counterparts in \code{r}. \code{copyv} does all that by first creating a copy of \code{x}, whereas \code{setv} modifies \code{x} directly and is thus more efficient.
}

\note{
(1) None of these functions currently support complex vectors.

(2) \code{setop} and the operators \code{\%+=\%}, \code{\%-=\%}, \code{\%*=\%} and \code{\%/=\%} also work with integer data, but do not perform any integer related checks. R's integers are bounded between +-2,147,483,647 and \code{NA_integer_} is stored as the value -2,147,483,648. Thus computations resulting in values exceeding +-2,147,483,647 will result in integer overflows, and \code{NA_integer_} should not occur on either side of a \code{setop} call. These are programmers functions and meant to provide the most efficient math possible to responsible users.

(3) It is possible to compare factors by the levels (e.g. \code{iris$Species \%==\% "setosa")}) or using integers (\code{iris$Species \%==\% 1L}). The latter is slightly more efficient. Nothing special is implemented for other objects apart from basic types, e.g. for dates (which are stored as doubles) you need to generate a date object i.e. \code{wlddev$date \%==\% as.Date("2019-01-01")}. Using \code{wlddev$date \%==\% "2019-01-01"} will give \code{integer(0)}.
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

## Missing values
mtc_na <- na_insert(mtcars, 0.15)    # Set 15\% of values missing at random
fnobs(mtc_na)                        # See observation count
na_omit(mtc_na)                      # 12x faster than na.omit(mtc_na)
na_omit(mtc_na, na.attr = TRUE)      # Adds attribute with removed cases, like na.omit
na_omit(mtc_na, cols = c("vs","am")) # Removes only cases missing vs or am
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
