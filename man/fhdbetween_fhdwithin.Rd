\name{fhdbetween-fhdwithin}
\alias{fhdbetween}
\alias{fhdbetween.default}
\alias{fhdbetween.matrix}
\alias{fhdbetween.data.frame}
\alias{fhdbetween.pseries}
\alias{fhdbetween.pdata.frame}
\alias{fhdwithin}
\alias{fhdwithin.default}
\alias{fhdwithin.matrix}
\alias{fhdwithin.data.frame}
\alias{fhdwithin.pseries}
\alias{fhdwithin.pdata.frame}
\alias{HDW}
\alias{HDW.default}
\alias{HDW.matrix}
\alias{HDW.data.frame}
\alias{HDW.pseries}
\alias{HDW.pdata.frame}
\alias{HDB}
\alias{HDB.default}
\alias{HDB.matrix}
\alias{HDB.data.frame}
\alias{HDB.pseries}
\alias{HDB.pdata.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Higher-Dimensional Centering and Linear Prediction
}
\description{
\code{fhdbetween} is a generalization of \code{fbetween} to efficiently predict with multiple factors and linear models (i.e. predict with vectors/factors, matrices, or data frames/lists where the latter may contain multiple factor variables). Similarly, \code{fhdwithin} is a generalization of \code{fwithin} to center on multiple factors and partial-out linear models.

The corresponding operators \code{HDB} and \code{HDW} additionally allow to predict / partial out full \code{lm()} formulas with interactions between variables.

}
\usage{
fhdbetween(x, \dots)
 fhdwithin(x, \dots)
       HDB(x, \dots)
       HDW(x, \dots)

\method{fhdbetween}{default}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", \dots)
\method{fhdwithin}{default}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", \dots)
\method{HDB}{default}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", \dots)
\method{HDW}{default}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", \dots)

\method{fhdbetween}{matrix}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", \dots)
\method{fhdwithin}{matrix}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, lm.method = "qr", \dots)
\method{HDB}{matrix}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, stub = "HDB.",
    lm.method = "qr", \dots)
\method{HDW}{matrix}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE, stub = "HDW.",
    lm.method = "qr", \dots)

\method{fhdbetween}{data.frame}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE,
           variable.wise = FALSE, lm.method = "qr", \dots)
\method{fhdwithin}{data.frame}(x, fl, w = NULL, na.rm = .op[["na.rm"]], fill = FALSE,
          variable.wise = FALSE, lm.method = "qr", \dots)
\method{HDB}{data.frame}(x, fl, w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]], fill = FALSE,
    variable.wise = FALSE, stub = "HDB.", lm.method = "qr", \dots)
\method{HDW}{data.frame}(x, fl, w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]], fill = FALSE,
    variable.wise = FALSE, stub = "HDW.", lm.method = "qr", \dots)

# Methods for indexed data / compatibility with plm:

\method{fhdbetween}{pseries}(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, \dots)
\method{fhdwithin}{pseries}(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, \dots)
\method{HDB}{pseries}(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, \dots)
\method{HDW}{pseries}(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE, \dots)

\method{fhdbetween}{pdata.frame}(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE,
           variable.wise = TRUE, \dots)
\method{fhdwithin}{pdata.frame}(x, effect = "all", w = NULL, na.rm = .op[["na.rm"]], fill = TRUE,
          variable.wise = TRUE, \dots)
\method{HDB}{pdata.frame}(x, effect = "all", w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]],
    fill = TRUE, variable.wise = TRUE, stub = "HDB.", \dots)
\method{HDW}{pdata.frame}(x, effect = "all", w = NULL, cols = is.numeric, na.rm = .op[["na.rm"]],
    fill = TRUE, variable.wise = TRUE, stub = "HDW.", \dots)
}

\arguments{
    \item{x}{a numeric vector, matrix, data frame, 'indexed_series' ('pseries') or 'indexed_frame' ('pdata.frame').}
  \item{fl}{a numeric vector, factor, matrix, data frame or list (which may or may not contain factors). In the \code{HDW/HDB} data frame method \code{fl} can also be a one-or two sided \code{lm()} formula with variables contained in \code{x}. Interactions \code{(:)} and full interactions \code{(*)} are supported. See Examples and the Note.}
  \item{w}{a vector of (non-negative) weights.}
  \item{cols}{\emph{data.frame methods}: Select columns to center (partial-out) or predict using column names, indices, a logical vector or a function. Unless specified otherwise all numeric columns are selected. If \code{NULL}, all columns are selected.}
  \item{na.rm}{remove missing values from both \code{x} and \code{fl}. by default rows with missing values in \code{x} or \code{fl} are removed. In that case an attribute "na.rm" is attached containing the rows removed.}
  \item{fill}{If \code{na.rm = TRUE}, \code{fill = TRUE} will not remove rows with missing values in \code{x} or \code{fl}, but fill them with \code{NA}'s.}
  \item{variable.wise}{\emph{(p)data.frame methods}: Setting \code{variable.wise = TRUE} will process each column individually i.e. use all non-missing cases in each column and in \code{fl} (\code{fl} is only checked for missing values if \code{na.rm = TRUE}). This is a lot less efficient but uses all data available in each column. }
\item{effect}{\emph{plm} methods: Select which panel identifiers should be used for centering. 1L takes the first variable in the \link[=indexing]{index}, 2L the second etc.. Index variables can also be called by name using a character vector. The keyword \code{"all"} uses all identifiers. }
\item{stub}{a prefix or stub to rename all transformed columns. \code{FALSE} will not rename columns.}
\item{lm.method}{character. The linear fitting method. Supported are \code{"chol"} and \code{"qr"}. See \code{\link{flm}}.}
  \item{\dots}{further arguments passed to \code{fixest::demean} (other than \code{notes} and \code{im_confident}) and \code{\link{chol}} / \code{\link{qr}}. Possible choices are \code{tol} to set a uniform numerical tolerance for the entire fitting process, or \code{nthreads} and \code{iter} to govern the higher-order centering process.}

}
\details{
\code{fhdbetween/HDB} and \code{fhdwithin/HDW} are powerful functions for high-dimensional linear prediction problems involving large factors and datasets, but can just as well handle ordinary regression problems. They are implemented as efficient wrappers around \code{\link[=fwithin]{fbetween / fwithin}}, \code{\link{flm}} and some C++ code from the \code{fixest} package that is imported for higher-order centering tasks (thus \code{fixest} needs to be installed for problems involving more than one factor).

Intended areas of use are to efficiently obtain residuals and predicted values from data, and to prepare data for complex linear models involving multiple levels of fixed effects. Such models can now be fitted using \code{(g)lm()} on data prepared with \code{fhdwithin / HDW} (relying on bootstrapped SE's for inference, or implementing the appropriate corrections). See Examples.

If \code{fl} is a vector or matrix, the result are identical to \code{lm} i.e. \code{fhdbetween / HDB} returns \code{fitted(lm(x ~ fl))} and \code{fhdwithin / HDW} \code{residuals(lm(x ~ fl))}. If \code{fl} is a list containing factors, all variables in \code{x} and non-factor variables in \code{fl} are centered on these factors using either \code{\link[=fwithin]{fbetween / fwithin}} for a single factor or \code{fixest} C++ code for multiple factors. Afterwards the centered data is regressed on the centered predictors. If \code{fl} is just a list of factors, \code{fhdwithin/HDW} returns the centered data and \code{fhdbetween/HDB} the corresponding means. Take as a most general example a list \code{fl = list(fct1, fct2, ..., var1, var2, ...)} where \code{fcti} are factors and \code{vari} are continuous variables. The output of \code{fhdwithin/HDW | fhdbetween/HDB} will then be identical to calling \code{resid | fitted} on \code{lm(x ~ fct1 + fct2 + ... + var1 + var2 + ...)}. The computations performed by \code{fhdwithin/HDW} and \code{fhdbetween/HDB} are however much faster and more memory efficient than \code{lm} because factors are not passed to \code{\link{model.matrix}} and expanded to matrices of dummies but projected out beforehand.

The formula interface to the data.frame method (only supported by the operators \code{HDW | HDB}) provides ease of use and allows for additional modeling complexity. For example it is possible to project out formulas like \code{HDW(data, ~ fct1*var1  + fct2:fct3 + var2:fct2:fct3 + var2:var3 + poly(var5,3)*fct5)} containing simple \code{(:)} or full \code{(*)} interactions of factors with continuous variables or polynomials of continuous variables, and two-or three-way interactions of factors and continuous variables. If the formula is one-sided as in the example above (the space left of \code{(~)} is left empty), the formula is applied to all variables selected through \code{cols}. The specification provided in \code{cols} (default: all numeric variables not used in the formula) can be overridden by supplying one-or more dependent variables. For example \code{HDW(data, var1 + var2 ~ fct1 + fct2)} will return a data.frame with \code{var1} and \code{var2} centered on \code{fct1} and \code{fct2}.

The special methods for 'indexed_series' (\code{plm::pseries}) and 'indexed_frame's (\code{plm::pdata.frame}) center a panel series or variables in a panel data frame on all panel-identifiers. By default in these methods \code{fill = TRUE} and \code{variable.wise = TRUE}, so missing values are kept. This change in the default arguments was done to ensure a coherent framework of functions and operators applied to \emph{plm} panel data classes.
}
\note{

% \subsection{Caution with full (*) and factor-continuous variable interactions:}{
% In general full interactions specified with \code{(*)} can be very slow on large data, and \code{lfe::demeanlist} is also not very speedy on interaction between factors and continuous variables, so these structures should be used with caution (don't just specify an interaction like that on a large dataset, start with smaller data and see how long computations take. Upon further updates of \code{lfe::demeanlist}, performance might improve).
% }

\subsection{On the differences between \code{fhdwithin/HDW}\dots and \code{fwithin/W}\dots:}{
\itemize{
\item \code{fhdwithin/HDW} can center data on multiple factors and also partial out continuous variables and factor-continuous interactions while \code{fwithin/W} only centers on one factor or the interaction of a set of factors, and does that very efficiently.
\item \code{HDW(data, ~ qF(group1) + qF(group2))} simultaneously centers numeric variables in data on \code{group1} and \code{group2}, while \code{W(data, ~ group1 + group2)} centers data on the interaction of \code{group1} and \code{group2}. The equivalent operation in \code{HDW} would be: \code{HDW(data, ~ qF(group1):qF(group2))}.
\item \code{W} always does computations on the variable-wise complete observations (in both matrices and data frames), whereas by default \code{HDW} removes all cases missing in either \code{x} or \code{fl}. In short, \code{W(data, ~ group1 + group2)} is actually equivalent to \code{HDW(data, ~ qF(group1):qF(group2), variable.wise = TRUE)}. \code{HDW(data, ~ qF(group1):qF(group2))} would remove any missing cases.
\item \code{fbetween/B} and \code{fwithin/W} have options to fill missing cases using group-averages and to add the overall mean back to group-demeaned data. These options are not available in \code{fhdbetween/HDB} and \code{fhdwithin/HDW}. Since \code{HDB} and \code{HDW} by default remove missing cases, they also don't have options to keep grouping-columns as in \code{B} and \code{W}.
}
}
}
\value{
\code{HDB} returns fitted values of regressing \code{x} on \code{fl}. \code{HDW} returns residuals. See Details and Examples.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link[=fbetween]{fbetween, fwithin}}, \code{\link{fscale}}, \code{\link{TRA}}, \code{\link{flm}}, \code{\link{fFtest}}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
HDW(mtcars$mpg, mtcars$carb)                   # Simple regression problems
HDW(mtcars$mpg, mtcars[-1])
HDW(mtcars$mpg, qM(mtcars[-1]))
head(HDW(qM(mtcars[3:4]), mtcars[1:2]))
head(HDW(iris[1:2], iris[3:4]))                # Partialling columns 3 and 4 out of colums 1 and 2
head(HDW(iris[1:2], iris[3:5]))                # Adding the Species factor -> fixed effect
\donttest{ % We don't test because this code depends on suggested package.
head(HDW(wlddev, PCGDP + LIFEEX ~ iso3c + qF(year))) # Partialling out 2 fixed effects
head(HDW(wlddev, PCGDP + LIFEEX ~ iso3c + qF(year), variable.wise = TRUE)) # Variable-wise
head(HDW(wlddev, PCGDP + LIFEEX ~ iso3c + qF(year) + ODA)) # Adding ODA as a continuous regressor
head(HDW(wlddev, PCGDP + LIFEEX ~ iso3c:qF(decade) + qF(year) + ODA)) # Country-decade and year FE's

head(HDW(wlddev, PCGDP + LIFEEX ~ iso3c*year))          # Country specific time trends
head(HDW(wlddev, PCGDP + LIFEEX ~ iso3c*poly(year, 3))) # Country specific cubic trends

# More complex examples
lm(HDW.mpg ~ HDW.hp, data = HDW(mtcars, ~ factor(cyl)*carb + vs + wt:gear + wt:gear:carb))
lm(mpg ~ hp + factor(cyl)*carb + vs + wt:gear + wt:gear:carb, data = mtcars)

lm(HDW.mpg ~ HDW.hp, data = HDW(mtcars, ~ factor(cyl)*carb + vs + wt:gear))
lm(mpg ~ hp + factor(cyl)*carb + vs + wt:gear, data = mtcars)

lm(HDW.mpg ~ HDW.hp, data = HDW(mtcars, ~ cyl*carb + vs + wt:gear))
lm(mpg ~ hp + cyl*carb + vs + wt:gear, data = mtcars)

lm(HDW.mpg ~ HDW.hp, data = HDW(mtcars, mpg + hp ~ cyl*carb + factor(cyl)*poly(drat,2)))
lm(mpg ~ hp + cyl*carb + factor(cyl)*poly(drat,2), data = mtcars)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line
\keyword{multivariate}
